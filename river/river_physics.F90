module river_physics_mod 

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! For the full text of the GNU General Public License,               
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="klf@gfdl.noaa.gov"> Kirsten Findell </CONTACT> 
! <CONTACT EMAIL="z1l@gfdl.noaa.gov"> Zhi Liang </CONTACT> 

  use mpp_mod,         only : mpp_sync_self, mpp_send, mpp_recv, EVENT_RECV, EVENT_SEND
  use mpp_mod,         only : mpp_npes, mpp_error, FATAL, mpp_get_current_pelist
  use mpp_mod,         only : mpp_root_pe, mpp_pe, mpp_max
  use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : ZERO, NINETY, MINUS_NINETY, mpp_update_domains 
  use mpp_domains_mod, only : mpp_get_compute_domains
  use mpp_domains_mod, only : mpp_get_num_overlap, mpp_get_overlap
  use fms_mod,         only : stdlog, open_namelist_file, write_version_number
  use fms_mod,         only : close_file, check_nml_error
  use river_type_mod,  only : river_type, Leo_Mad_trios
  use constants_mod,   only : tfreeze, hlf

  implicit none
  private


!--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river_physics.F90,v 16.0 2008/07/30 22:13:03 fms Exp $'
  character(len=128) :: tagname = '$Name: perth_2008_10 $'


! ---- public interfaces -----------------------------------------------------

  public :: river_physics_init, river_physics_step

!----------------------------------------------------------------------
  real               :: clw = 4218.
  real               :: csw = 2106.
  real,    parameter :: sec_in_day = 86400.

! ---- namelist interface
  character*6 :: algor = 'linear'

  namelist /river_physics_nml/ algor

  integer, parameter, dimension(8) :: di=(/1,1,0,-1,-1,-1,0,1/)
  integer, parameter, dimension(8) :: dj=(/0,-1,-1,-1,0,1,1,1/)
  integer                          :: isc, iec, jsc, jec  ! compute domain
  integer                          :: isd, ied, jsd, jed  ! data domain
  integer                          :: maxtravel
  integer                          :: npes
  integer                          :: num_species
 
  type comm_type
     integer          :: count
     integer          :: pe
     integer, pointer :: i(:) => NULL()
     integer, pointer :: j(:) => NULL()
     integer, pointer :: k(:) => NULL()
  end type comm_type

  type halo_update_type
     type(comm_type), pointer :: send(:) => NULL();
     type(comm_type), pointer :: recv(:) => NULL();
  end type halo_update_type

  type(halo_update_type),  allocatable :: halo_update(:)
  real, dimension(:),      allocatable :: send_buffer, recv_buffer
  logical, dimension(:,:), allocatable :: in_domain
  integer, dimension(:,:), allocatable :: nlev

contains

!#######################################################################

  subroutine river_physics_init(River, domain )
    type(river_type), intent(inout) :: River
    type(domain2d),   intent(inout) :: domain
    integer                         :: unit, io_status, ierr
    integer                         :: i, j


!--- read namelist -------------------------------------------------
    unit = open_namelist_file()
    ierr = 1;
    do while ( ierr/=0 )
       read  (unit, river_physics_nml, iostat=io_status, end=10)
       ierr = check_nml_error(io_status,'river_physics_nml')
    enddo
10  call close_file (unit)

!--- write version and namelist info to logfile --------------------
    call write_version_number(version,tagname)
    write (stdlog(), river_physics_nml)  

    npes     = mpp_npes()

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)

    num_species = size(River%outflow_c,3)
    maxtravel = maxval(River%travel)
    call mpp_max(maxtravel)

!--- set up the halo update 
    call setup_halo_update(River, domain)

    do j = jsc, jec
       do i = isc, iec
          if(River%tocell(i,j) > 0) then
             River%i_tocell(i,j) = i + di(River%tocell(i,j))
             River%j_tocell(i,j) = j + dj(River%tocell(i,j))
          end if
       end do
    end do
 

  end subroutine river_physics_init

!#####################################################################

  subroutine river_physics_step(River, cur_travel)

    type(river_type),     intent(inout) :: River
    integer,                 intent(in) :: cur_travel

! ---- local vars ----------------------------------------------------------
    integer   :: i, j
    real      :: Q0, dQ_dV, avail, out_frac, qmelt
    real      :: v_r_d(River%num_species-River%num_c+1:River%num_species)
    real      :: conc(1:River%num_species)


    ! do for all cells at current number of steps from river mouth
    do j = jsc, jec 
       do i = isc, iec
          if (River%travel(i,j)==cur_travel) then

! avail is amount of water to be partitioned between outflow and new storage
              avail = River%storage(i,j)  &
                   +(River%inflow(i,j)+River%infloc(i,j))*River%dt_slow
! determine total water storage at end of step
              if (algor.eq.'linear') then
! linear algorithm assumes outflow = Q0+dQ_dV*dS
                  if (River%storage(i,j) .le. 0.) then
                      Q0 = 0.; dQ_dV = 0.
                  else
                      Q0=River%o_coef(i,j)*River%storage(i,j)**River%o_exp
                      dQ_dV=River%o_exp*Q0/River%storage(i,j)
                  endif
                  River%storage(i,j) = River%storage(i,j) + River%dt_slow *   &
                       (River%inflow(i,j)+River%infloc(i,j)-Q0) &
                       /(1.+River%dt_slow*dQ_dV)
              else if (algor.eq.'nonlin') then
! nonlin algorithm assumes all inflow at start of step 
! and then integrates analytically
                  if (avail .gt. 0.) then
                      River%storage(i,j) = (avail**(1.-River%o_exp) &
                           + River%o_coef(i,j)*(River%o_exp-1.)*River%dt_slow) &
                           **(1./(1.-River%o_exp))
                  else
                      River%storage(i,j) = avail
                  endif
              endif
! determine total water outflow during step
              River%outflow(i,j) = (avail - River%storage(i,j)) / River%dt_slow
! given outflow, determine flow width, depth, velocity
              if (River%outflow(i,j) .le. 0.) then
                  River%depth(i,j) = 0.
                  River%width(i,j) = 0.
                  River%vel(i,j)   = 0.
              else
                  River%depth(i,j) = River%d_coef(i,j) &
                       * River%outflow(i,j)**River%d_exp
                  River%width(i,j) = River%w_coef(i,j) &
                       * River%outflow(i,j)**River%w_exp
                  River%vel(i,j) = River%outflow(i,j) /                   &
                                        (River%width(i,j) * River%depth(i,j))
              endif
! given water outflow and storage, apportion other tracked stuff in same ratio
              out_frac = 0.
              if (avail .gt. 0.) out_frac = River%outflow(i,j)/avail
              River%outflow_c(i,j,:) = out_frac * (River%storage_c(i,j,:) &
                   +(River%inflow_c(i,j,:)+River%infloc_c(i,j,:))*River%dt_slow)
! add outflows to the inputs for downstream cells

              River%storage_c(i,j,:) = River%storage_c(i,j,:)       &
                      + (River%inflow_c(i,j,:) - River%outflow_c(i,j,:)   &
                      + River%infloc_c(i,j,:) ) * River%dt_slow
! define intensive variables for diagnostics and for use in transformations.
! along the way, melt swept snow as necessary. freeze will be a separate
! process, added later; it will be different in that frozen river water will
! be stationary, thus a different species
              if (River%storage(i,j) .gt. 0.) then
                  conc(1) = River%storage_c(i,j,1)/River%storage(i,j)
                  conc(2) = tfreeze + River%storage_c(i,j,2) /  &
                       ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
                  if (River%storage_c(i,j,1).gt.0. .and. conc(2).gt.tfreeze) then
                      qmelt = min(hlf*River%storage_c(i,j,1), River%storage_c(i,j,2))
                      River%storage_c(i,j,1) = River%storage_c(i,j,1) - qmelt/hlf
                      River%storage_c(i,j,2) = River%storage_c(i,j,2) - qmelt
                      conc(2) = tfreeze + River%storage_c(i,j,2) /  &
                           ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
                   endif
              else
                  conc(1) = 0.
                  conc(2) = tfreeze   ! TEMPORARY
              endif

              if (River%do_age) then
                  River%removal_c(i,j,River%num_phys+1) = -River%storage(i,j)/sec_in_day
                  River%storage_c(i,j,River%num_phys+1) = River%storage_c(i,j,River%num_phys+1) &
                       - River%removal_c(i,j,River%num_phys+1)*River%dt_slow
              endif

              if (River%storage(i,j) .gt. 0.) then
                  conc(River%num_phys+1:River%num_species) = &
                       River%storage_c(i,j,River%num_phys+1:River%num_species)/River%storage(i,j)
              else
                  conc(River%num_phys+1:River%num_species) = 0.
              endif

              if(River%num_c.gt.0) then
                  if (River%depth(i,j).gt.0.) then
                      v_r_d = River%vf_ref * River%Q10**((conc(2)-River%t_ref)/10.)&
                           / ((1+River%kinv*conc(River%num_species-River%num_c+1:River%num_species)) &
                           *River%depth(i,j))
! next should not be necessary if storage_c is positive, but maybe it's not.
                      v_r_d = River%vf_ref * River%Q10**((conc(2)-River%t_ref)/10.)&
                           / ((1+River%kinv*max(0.,conc(River%num_species-River%num_c+1:River%num_species)))*River%depth(i,j))
                  else
                      v_r_d = 0.
                  endif
                  River%removal_c(i,j,River%num_species-River%num_c+1:River%num_species) = &
                       River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) &
                       * (1-exp( -v_r_d * River%dt_slow)) &
                       / River%dt_slow
                  River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) = &
                       River%storage_c(i,j,River%num_species-River%num_c+1:River%num_species) &
                       - River%removal_c(i,j,River%num_species-River%num_c+1:River%num_species)* River%dt_slow
              endif

          endif
       enddo
    enddo

    call do_halo_update(River, halo_update(cur_travel))

  end subroutine river_physics_step

!#####################################################################
  subroutine setup_halo_update(River, domain)
    type(river_type),         intent(in) :: River
    type(domain2d),        intent(inout) :: domain

    integer, parameter                   :: MAXCOMM      = 8  ! should be no larger than 8.
    integer                              :: travelnow, nsend2, p, toc
    integer                              :: spos, rpos, n, m, l
    integer                              :: buffer_pos, pos, msgsize
    integer                              :: i, j, i1, j1, i2, j2, i3, j3, i4, j4, k, kk
    integer                              :: send_size, recv_size, siz, i_dest, j_dest
    integer, allocatable, dimension(:,:) :: tocell
    integer, allocatable, dimension(:,:) :: is_recv, ie_recv, js_recv, je_recv
    integer, allocatable, dimension(:,:) :: is1_send, ie1_send, js1_send, je1_send
    integer, allocatable, dimension(:,:) :: is2_send, ie2_send, js2_send, je2_send
    integer, allocatable, dimension(:,:) :: rot_send, rot_recv, dir_send, dir_recv
    integer, allocatable, dimension(:)   :: nsend, nrecv, pelist
    integer, allocatable, dimension(:)   :: isl, iel, jsl, jel  
    integer, allocatable, dimension(:)   :: sbuf, rbuf
    type(comm_type), pointer             :: send => NULL()
    integer, allocatable, dimension(:,:,:) :: i_send, j_send, t_send, p_send, n_send

    call mpp_get_data_domain   (domain, isd, ied, jsd, jed)
    allocate(pelist  (0:npes-1        )                             )
    call mpp_get_current_pelist(pelist)
    
    allocate(halo_update(maxtravel) )
    do travelnow = 1, maxtravel
       allocate(halo_update(travelnow)%send(0:npes-1))
       allocate(halo_update(travelnow)%recv(0:npes-1))
       halo_update(travelnow)%send(:)%count = 0
       halo_update(travelnow)%recv(:)%count = 0
       do p = 0, npes-1
          halo_update(travelnow)%send(p)%count = 0
          halo_update(travelnow)%recv(p)%count = 0
          halo_update(travelnow)%send(p)%pe    = pelist(p)
          halo_update(travelnow)%recv(p)%pe    = pelist(p)
       end do
    end do

    !--- first get the travel and tocell information onto data domain
    allocate(tocell(isd:ied,jsd:jed))
    allocate(isl(0:npes-1), iel(0:npes-1), jsl(0:npes-1), jel(0:npes-1) )
    tocell(isc:iec,jsc:jec) = River%tocell(isc:iec,jsc:jec)
    call mpp_update_domains(tocell, domain)
    call mpp_get_compute_domains(domain, xbegin=isl, xend=iel, ybegin=jsl, yend=jel)

    !--- first get the halo update information for send and recv.
    allocate(is_recv (0:npes-1,MAXCOMM), ie_recv (0:npes-1,MAXCOMM) )
    allocate(js_recv (0:npes-1,MAXCOMM), je_recv (0:npes-1,MAXCOMM) )
    allocate(is1_send(0:npes-1,MAXCOMM), ie1_send(0:npes-1,MAXCOMM) )
    allocate(js1_send(0:npes-1,MAXCOMM), je1_send(0:npes-1,MAXCOMM) )
    allocate(is2_send(0:npes-1,MAXCOMM), ie2_send(0:npes-1,MAXCOMM) )
    allocate(js2_send(0:npes-1,MAXCOMM), je2_send(0:npes-1,MAXCOMM) )
    allocate(rot_recv(0:npes-1,MAXCOMM), rot_send(0:npes-1,MAXCOMM) )
    allocate(dir_recv(0:npes-1,MAXCOMM), dir_send(0:npes-1,MAXCOMM) )
    allocate(rbuf    (npes*MAXCOMM*4  ), sbuf    (npes*MAXCOMM*4  ) )
    allocate(nrecv   (0:npes-1        ), nsend   (0:npes-1        ) )

    spos = 0
    do p = 0, npes-1
       nrecv(p) = mpp_get_num_overlap(domain, EVENT_RECV, p)
       call mpp_send(nrecv(p),     plen = 1,       to_pe = pelist(p))
       if(nrecv(p) >0) then
          if(nrecv(p) > MAXCOMM) call mpp_error(FATAL, &
          "river_mod: number of overlapping for recving is larger than MAXCOMM, increase MAXCOMM")
          call mpp_get_overlap(domain, EVENT_RECV, p, is_recv(p,1:nrecv(p)), ie_recv(p,1:nrecv(p)), &
                               js_recv(p,1:nrecv(p)), je_recv(p,1:nrecv(p)), dir_recv(p,1:nrecv(p)), rot_recv(p,1:nrecv(p)) )
          !--- send the information to the process that send data.
          do n = 1, nrecv(p)
             sbuf(spos+(n-1)*4+1) = is_recv(p,n)
             sbuf(spos+(n-1)*4+2) = ie_recv(p,n)
             sbuf(spos+(n-1)*4+3) = js_recv(p,n)
             sbuf(spos+(n-1)*4+4) = je_recv(p,n)
          end do

          call mpp_send(sbuf(spos+1), plen = 4*nrecv(p), to_pe = pelist(p))
          spos = spos + 4*nrecv(p)
       end if
    end do

    rpos = 0
    do p = 0, npes-1
       nsend(p) = mpp_get_num_overlap(domain, EVENT_SEND, p)
       call mpp_recv(nsend2,     glen = 1,       from_pe = pelist(p))
       if(nsend(p) .NE. nsend2) call mpp_error(FATAL, &
              "river_mod: number of send and recv between two processors are not equal")
       if(nsend(p) >0) then
          if(nsend(p) > MAXCOMM) call mpp_error(FATAL, &
          "river_mod: number of overlapping for sending is larger than MAXCOMM, increase MAXCOMM")
          call mpp_get_overlap(domain, EVENT_SEND, p, is1_send(p,1:nsend(p)), ie1_send(p,1:nsend(p)), &
                               js1_send(p,1:nsend(p)), je1_send(p,1:nsend(p)), dir_send(p,1:nsend(p)), rot_send(p,1:nsend(p)) )
          call mpp_recv(rbuf(rpos+1), glen=4*nsend(p), from_pe=pelist(p))
          do n = 1, nsend(p)
             is2_send(p,n) = rbuf(rpos+(n-1)*4+1)
             ie2_send(p,n) = rbuf(rpos+(n-1)*4+2)
             js2_send(p,n) = rbuf(rpos+(n-1)*4+3)
             je2_send(p,n) = rbuf(rpos+(n-1)*4+4)
          end do
          rpos = rpos + 4*nsend(p)           
       end if
    end do

    call mpp_sync_self()

    do p = 0, npes-1
       !--- configure points need to receive from other pe.
       !--- (i,j) --- halo index on the pe sent to, one neighbor pe data domain
       !--- (i1,j1) --- neighbor index of (i,j) on the pe sent to, on neighbor pe compute domain
       !--- (i2,j2) --- my index corresponding to (i,j), on my compute domain
       !--- (i3,j3) --- neighbor index of (i2,j2), on my data domain
       !--- (i4,j4) --- index of (i1,j1) tocell. 
       do n = 1, nsend(p)
          select case ( dir_send(p,n) )
          case(1)  ! east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             i1 = i - 1
             do j = js2_send(p,n), je2_send(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsl(p) .OR. j1 > jel(p) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! w->e
                      i2 = is1_send(p,n)
                      i3 = i2 -1
                      j2 = js1_send(p,n) + j  - js2_send(p,n)
                      j3 = js1_send(p,n) + j1 - js2_send(p,n)
                   case (NINETY) ! s->e
                      i2 = is1_send(p,n) + (je2_send(p,n) - j )
                      i3 = is1_send(p,n) + (je2_send(p,n) - j1)
                      j2 = js1_send(p,n)
                      j3 = j2 - 1
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                      end if
                   end if
                end do
             end do
          case(2)  ! south east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i - 1
             j1 = j + 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! nw->se 
                i3 = i2 - 1
                j3 = j2 + 1
             case (NINETY)
                i3 = i2 - 1
                j3 = j2 - 1
             case (MINUS_NINETY)
                i3 = i2 + 1
                j3 = j2 + 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                end if
             end if
          case(3)  ! south
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             j1 = j + 1
             do i = is2_send(p,n), ie2_send(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isl(p) .OR. i1 > iel(p) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! n->s
                      i2 = is1_send(p,n) + i  - is2_send(p,n)
                      i3 = is1_send(p,n) + i1 - is2_send(p,n)
                      j2 = js1_send(p,n)
                      j3 = j2 + 1
                   case (MINUS_NINETY) ! e->s
                      i2 = is1_send(p,n)
                      i3 = i2 + 1
                      j2 = js1_send(p,n) + (ie2_send(p,n) - i )
                      j3 = js1_send(p,n) + (ie2_send(p,n) - i1)
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                      end if
                   end if
                end do
             end do
          case(4)  ! south west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i + 1
             j1 = j + 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! ne->sw 
                i3 = i2 + 1
                j3 = j2 + 1
             case (NINETY) !  
                i3 = i2 - 1
                j3 = j2 + 1
             case (MINUS_NINETY)
                i3 = i2 + 1
                j3 = j2 - 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                end if
             end if
          case(5)  ! west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             i1 = i + 1
             do j = js2_send(p,n), je2_send(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsl(p) .OR. j1 > jel(p) ) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! e->w
                      i2 = is1_send(p,n)
                      i3 = i2 + 1
                      j2 = js1_send(p,n) + j  - js2_send(p,n)
                      j3 = js1_send(p,n) + j1 - js2_send(p,n)
                   case (NINETY) ! n->w
                      i2 = is1_send(p,n) + (je2_send(p,n) - j )
                      i3 = is1_send(p,n) + (je2_send(p,n) - j1)
                      j2 = js1_send(p,n)
                      j3 = j2 + 1
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                      end if
                   end if
                end do
             end do
          case(6)  ! north west
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i + 1
             j1 = j - 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! se->nw 
                i3 = i2 + 1
                j3 = j2 - 1
             case (NINETY) !  
                i3 = i2 + 1
                j3 = j2 + 1
             case (MINUS_NINETY)
                i3 = i2 - 1
                j3 = j2 - 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                end if
             end if
          case(7)  ! north
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             j1 = j - 1
             do i = is2_send(p,n), ie2_send(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isl(p) .OR. i1 > iel(p)) cycle
                   select case(rot_send(p,n))
                   case (ZERO) ! s->n
                      i2 = is1_send(p,n) + i  - is2_send(p,n)
                      i3 = is1_send(p,n) + i1 - is2_send(p,n)
                      j2 = js1_send(p,n)
                      j3 = j2 - 1
                   case (MINUS_NINETY) ! w->n
                      i2 = is1_send(p,n)
                      i3 = i2 - 1
                      j2 = js1_send(p,n) + (ie2_send(p,n) - i )
                      j3 = js1_send(p,n) + (ie2_send(p,n) - i1)
                   end select
                   if(River%travel(i3,j3) >0) then
                      toc = tocell(i3,j3)
                      i4 = i1 + di(toc)
                      j4 = j1 + dj(toc) 
                      if(i4 == i .AND. j4 == j) then
                         call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                      end if
                   end if
                end do
             end do
          case(8)  ! north east
             i = is2_send(p,n)  ! is2_send(p,n) == ie2_send(p,n)
             j = js2_send(p,n)  ! js2_send(p,n) == je2_send(p,n)
             i1 = i - 1
             j1 = j - 1
             i2 = is1_send(p,n)  ! is1_send(p,n) == ie1_send(p,n)
             j2 = js1_send(p,n)  ! js1_send(p,n) == je1_send(p,n)
             select case(rot_send(p,n))
             case (ZERO) ! sw->ne 
                i3 = i2 - 1
                j3 = j2 - 1
             case (NINETY) !  
                i3 = i2 + 1
                j3 = j2 - 1
             case (MINUS_NINETY)
                i3 = i2 - 1
                j3 = j2 + 1
             end select  
             if(River%travel(i3,j3) >0) then
                toc = tocell(i3,j3)
                i4 = i1 + di(toc)
                j4 = j1 + dj(toc) 
                if(i4 == i .AND. j4 == j) then
                   call add_single_overlap(halo_update(River%travel(i3,j3))%recv(p), i2, j2)
                end if
             end if
          end select
       end do       

       !--- configure points need to send to other pe.
       !--- (i,j) ---  index on my data domain
       !--- (i1,j1) --- index on my compute domain corresponding to (i,j)
       !--- (i2,j2) --- index of (i1,j1) tocell
       do n = 1, nrecv(p)
          select case ( dir_recv(p,n) )
          case(1)  ! east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             i1 = i - 1
             do j = js_recv(p,n), je_recv(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsc .OR. j1 > jec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                   end if
                end do
             end do
          case(2)  ! south east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i - 1
             j1 = j + 1
             if(River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                end if
             end if
          case(3)  ! south
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             j1 = j + 1
             do i = is_recv(p,n), ie_recv(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isc .OR. i1 > iec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                   end if
                end do
             end do
          case(4)  ! south west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i + 1
             j1 = j + 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                end if
             end if
          case(5)  ! west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             i1 = i + 1
             do j = js_recv(p,n), je_recv(p,n)
                do l = -1,1
                   j1 = j + l
                   if(j1<jsc .OR. j1 > jec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                   end if
                end do
             end do
          case(6)  ! north west
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i + 1
             j1 = j - 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                end if
             end if
          case(7)  ! north
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             j1 = j - 1
             do i = is_recv(p,n), ie_recv(p,n)
                do l = -1,1
                   i1 = i + l
                   if(i1<isc .OR. i1 > iec .OR. River%travel(i1,j1) < 1) cycle
                   i2 = i1 + di(tocell(i1,j1))
                   j2 = j1 + dj(tocell(i1,j1)) 
                   if(i2 == i .AND. j2 == j) then
                      call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                   end if
                end do
             end do
          case(8)  ! north east
             i = is_recv(p,n)  ! is_recv(p,n) == ie_recv(p,n)
             j = js_recv(p,n)  ! js_recv(p,n) == je_recv(p,n)
             i1 = i - 1
             j1 = j - 1
             if( River%travel(i1,j1) > 0) then
                i2 = i1 + di(tocell(i1,j1))
                j2 = j1 + dj(tocell(i1,j1)) 
                if(i2 == i .AND. j2 == j) then
                   call add_single_overlap(halo_update(River%travel(i1,j1))%send(p), i1, j1)
                end if
             end if
          end select
       end do       
    end do ! end do p = 0, npes-1

    allocate(in_domain(isc:iec,jsc:jec))
    in_domain = .true.
    do j = jsc, jec
       do i = isc, iec
          if(River%tocell(i,j) > 0) then
             i_dest = i + di(River%tocell(i,j))
             j_dest = j + dj(River%tocell(i,j))
             if(i_dest < isc .OR. i_dest > iec .OR. j_dest < jsc .OR. j_dest > jec) then 
                LOOP_TRAVEL: do travelnow = 1, maxtravel
                   do p = 0, npes-1
                      send => halo_update(travelnow)%send(p)
                      do n = 1, send%count
                         if(send%i(n) == i .AND. send%j(n) == j) then
                            in_domain(i,j) = .false.
                            exit LOOP_TRAVEL 
                         end if
                      end do
                   end do
                end do LOOP_TRAVEL
                if(in_domain(i,j)) then
                   i_dest = i
                   j_dest = j
                end if
             end if
             River%i_tocell(i,j) = i_dest
             River%j_tocell(i,j) = j_dest
          end if
       end do
    end do

    !--- add points that sent to self.
    p = mpp_pe() - mpp_root_pe()
    do j = jsc, jec
       do i = isc, iec
          m = River%travel(i,j)
          if(m >0 .and. in_domain(i,j) ) then
             call add_single_overlap(halo_update(m)%send(p), i, j)
             call add_single_overlap(halo_update(m)%recv(p), River%i_tocell(i, j), River%j_tocell(i, j))
          end if
       end do
    end do

    !--- the following is for the purpose of bitwise reproduce between processor count
    send_size = 0
    do p = 0, npes-1
       do m = 1, maxtravel
          send_size = send_size + halo_update(m)%send(p)%count
       end do
    end do
    send_size = send_size*2
    allocate(send_buffer(send_size))
    pos = 0
    do p=0, npes-1
       buffer_pos = pos
       do m = 1, maxtravel
          do n = 1, halo_update(m)%send(p)%count
             send_buffer(pos+1) = halo_update(m)%send(p)%i(n)
             send_buffer(pos+2) = halo_update(m)%send(p)%j(n)
             pos = pos + 2
          end do
       end do
       msgsize = pos - buffer_pos
       call mpp_send(msgsize, plen=1, to_pe = pelist(p))
       if(msgsize >0) then
          call mpp_send(send_buffer(buffer_pos+1), plen = msgsize, to_pe = pelist(p))
       end if      
    end do

    allocate(i_send(isc:iec,jsc:jec,8), j_send(isc:iec,jsc:jec,8) )
    allocate(p_send(isc:iec,jsc:jec,8), t_send(isc:iec,jsc:jec,8) )
    allocate(n_send(isc:iec,jsc:jec,8), nlev(isc:iec,jsc:jec))
    nlev = 0
    do p=0, npes-1
       recv_size = 0
       do m = 1, maxtravel
          recv_size = recv_size + halo_update(m)%recv(p)%count
       end do
       recv_size = recv_size*2
       call mpp_recv(msgsize, glen = 1, from_pe = pelist(p) )
       if(msgsize .NE. recv_size) call mpp_error(FATAL, "river_physics_mod: mismatch at send size and recv size")

       if(recv_size >0) then
          allocate(recv_buffer(recv_size))
          call mpp_recv(recv_buffer(1), glen = recv_size, from_pe = pelist(p))

          pos = 0
          do m = 1, maxtravel
             do n = 1, halo_update(m)%recv(p)%count
                i = halo_update(m)%recv(p)%i(n)
                j = halo_update(m)%recv(p)%j(n)
                i1 = recv_buffer(pos+1)
                j1 = recv_buffer(pos+2)
                pos = pos + 2
                do k = 1, nlev(i,j)
                   if( j1 < j_send(i,j,k) .OR. (j1 == j_send(i,j,k) .AND. i1 < i_send(i,j,k) ) ) then
                      do kk = nlev(i,j)+1, k+1, -1
                         i_send(i,j,kk) = i_send(i,j,kk-1)
                         j_send(i,j,kk) = j_send(i,j,kk-1)
                         p_send(i,j,kk) = p_send(i,j,kk-1)
                         t_send(i,j,kk) = t_send(i,j,kk-1)
                         n_send(i,j,kk) = n_send(i,j,kk-1)
                         halo_update(t_send(i,j,kk))%recv(p_send(i,j,kk))%k(n_send(i,j,kk)) = &
                             halo_update(t_send(i,j,kk))%recv(p_send(i,j,kk))%k(n_send(i,j,kk)) + 1
                      end do
                      exit
                   end if
                end do
                nlev(i,j) = nlev(i,j) + 1
                i_send(i,j,k) = i1
                j_send(i,j,k) = j1
                p_send(i,j,k) = p
                t_send(i,j,k) = m
                n_send(i,j,k) = n         
                halo_update(m)%recv(p)%k(n) = k               
             end do
          end do
          deallocate(recv_buffer)
       end if
    end do

    call mpp_sync_self()    
    deallocate(send_buffer)

    !--- set up buffer for send and recv.
    send_size = 0
    do m = 1, maxtravel
       siz = 0
       do p = 0, npes-1
          send_size = send_size + halo_update(m)%send(p)%count
       end do
    end do
    send_size = send_size*(num_species+1)
    if(send_size > 0) allocate(send_buffer(send_size))

    recv_size = 0
    do m = 1, maxtravel
       siz = 0
       do p = 0, npes-1
          siz = siz + halo_update(m)%recv(p)%count
       end do
       recv_size = max(recv_size, siz)
    end do
    recv_size = recv_size*(num_species+1)
    if(recv_size > 0) allocate(recv_buffer(recv_size))

    deallocate(tocell)
    deallocate(isl, iel, jsl, jel )
    deallocate(sbuf, rbuf)
    deallocate(is_recv, ie_recv, js_recv, je_recv)
    deallocate(is1_send, ie1_send, js1_send, je1_send)
    deallocate(is2_send, ie2_send, js2_send, je2_send)
    deallocate(rot_send, rot_recv, nsend, nrecv, pelist)
    return

  end subroutine setup_halo_update

!###############################################################################
  subroutine do_halo_update(River, update)
     type(river_type),    intent(inout) :: River
     type(halo_update_type), intent(in) :: update
     type(comm_type), pointer           :: send=>NULL()
     type(comm_type), pointer           :: recv=>NULL()
     integer                            :: buffer_pos, pos
     integer                            :: p, n, i, j, count, l, k
     real                               :: wrk_c(isc:iec,jsc:jec,num_species, 8)
     real                               :: wrk  (isc:iec,jsc:jec, 8)

     !--- send the data
     pos = 0
     do p = 0, npes-1
        send => update%send(p)
        count = send%count
        if(count == 0) cycle
        buffer_pos = pos
        do n = 1, count
           i = send%i(n)
           j = send%j(n)
           pos = pos + 1
           send_buffer(pos)   = River%outflow(i,j)
           do l = 1, num_species
              pos = pos + 1
              send_buffer(pos) = River%outflow_c(i,j,l)
           end do
        end do
        call mpp_send(send_buffer(buffer_pos+1), plen=count*(num_species+1), to_pe = send%pe ) 
     end do

     !--- receive the data, must receive the data in the order 0, 1, 2, ... npes-1 to bitwise reproducing.
     nlev = 0
     do p = 0, npes-1
        recv => update%recv(p)
        count = recv%count
        if(count == 0) cycle
        call mpp_recv(recv_buffer(1), glen=count*(num_species+1), from_pe=recv%pe ) 
        pos = 0
        do n = 1, count
           i = recv%i(n)
           j = recv%j(n)
           k = recv%k(n)
           pos = pos + 1
           wrk(i,j,k) = recv_buffer(pos)
           nlev(i,j) = nlev(i,j)+1
           do l = 1, num_species
              pos = pos + 1
              wrk_c(i,j,l,k) = recv_buffer(pos)
           end do
        end do
     end do

     do j = jsc, jec
        do i = isc, iec
           do k = 1, nlev(i,j)
              River%inflow(i,j)   = River%inflow(i,j) + wrk(i,j,k)
              River%inflow_c(i,j,:) = River%inflow_c(i,j,:) + wrk_c(i,j,:,k)
           end do
        end do
     end do

     call mpp_sync_self()


     return

  end subroutine do_halo_update


!###############################################################################
!  This routine will add one point for send/recv into the data type, allocate 
!  memory or expand memory if needed 

  subroutine add_single_overlap(comm, i, j)
    type(comm_type), intent(inout) :: comm
    integer,         intent(in)    :: i, j
    integer, parameter             :: DEFAULT_SIZE = 40 ! too big or too small
    integer                        :: count, maxcount
    integer, allocatable           :: tmp(:)

    count = comm%count 
    count = count + 1
    if(count == 1) then ! assign default space to hold index
       allocate(comm%i(DEFAULT_SIZE))
       allocate(comm%j(DEFAULT_SIZE))
       allocate(comm%k(DEFAULT_SIZE))
    end if
    maxcount = size(comm%i)

    if(count > maxcount) then ! need to expend the size to hold the index
       allocate(tmp(maxcount))
       tmp = comm%i
       deallocate(comm%i)
       allocate(comm%i(2*maxcount))
       comm%i(1:maxcount) = tmp
       tmp = comm%j
       deallocate(comm%j)
       allocate(comm%j(2*maxcount))
       comm%j(1:maxcount) = tmp
       deallocate(tmp)
       deallocate(comm%k)
       allocate(comm%k(2*maxcount))
    end if
    comm%i(count) = i
    comm%j(count) = j    
    comm%count    = count

    return    

  end subroutine add_single_overlap

!#####################################################################

end module river_physics_mod
