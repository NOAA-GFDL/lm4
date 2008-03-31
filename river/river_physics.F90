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

  use fms_mod,         only : stdlog, open_namelist_file, write_version_number
  use fms_mod,         only : close_file, check_nml_error
  use river_type_mod,  only : river_type, Leo_Mad_trios
  use constants_mod,   only : DENS_H2O, tfreeze, hlf

  implicit none
  private


  !--- version information ---------------------------------------------
  character(len=128) :: version = '$Id: river_physics.F90,v 15.0.2.4 2007/12/05 19:41:35 slm Exp $'
  character(len=128) :: tagname = '$Name: omsk_2008_03 $'


  ! ---- public interfaces -----------------------------------------------------

  public :: river_physics_init, river_physics_step

  !----------------------------------------------------------------------
  real               :: clw = 4218.
  real               :: csw = 2106.
  real,    parameter :: sec_in_day = 86400.
  real               :: dens_inv
  integer, parameter, dimension(8) :: di=(/1,1,0,-1,-1,-1,0,1/)
  integer, parameter, dimension(8) :: dj=(/0,-1,-1,-1,0,1,1,1/)
  integer idest,jdest

  ! ---- namelist interface
  character*6 :: algor = 'linear'

  namelist /river_physics_nml/ algor


contains

  !#######################################################################

  subroutine river_physics_init( )

    integer            :: unit, io_status, ierr, siz(4)
    integer            :: isc, iec, jsc, jec       ! domain decomposition of land grid
    character(len=128) :: filename

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

    dens_inv = 1.0/DENS_H2O

  end subroutine river_physics_init

  !#####################################################################

  subroutine river_physics_step(River, petravel,cur_travel)

    type(river_type),     intent(inout) :: River
    integer, dimension(:,:), intent(in) :: petravel
    integer,                 intent(in) :: cur_travel

    ! ---- local vars ----------------------------------------------------------
    integer   :: i, j, i_species
    real      :: Q0, dQ_dV, avail, out_frac, qmelt
    real      :: v_r_d(River%num_species-River%num_c+1:River%num_species), &
                 conc(1:River%num_species)
    
! do for all cells at current number of steps from river mouth
    do j = River%js(cur_travel),River%je(cur_travel) 
       do i = River%is(cur_travel),River%ie(cur_travel)
          if (petravel(i,j)==cur_travel) then

! avail is amount of water to be partitioned between outflow and new storage
              avail = River%storage(i,j)  &
                 +(River%inflow(i,j)+River%infloc(i,j))*River%dt_slow
! determine total water storage at end of step
              if (algor.eq.'linear') then
                  ! linear algorithm assumes outflow = Q0+dQ_dV*dS
                  if (River%storage(i,j) .le. 0.) then
                      Q0 = 0.; dQ_dV = 0.
                    else
                      Q0=River%o_coef(i,j)*River%storage(i,j)**River%o_exp(i,j)
                      dQ_dV=River%o_exp(i,j)*Q0/River%storage(i,j)
                    endif
                  River%storage(i,j) = River%storage(i,j) + River%dt_slow *   &
                    (River%inflow(i,j)+River%infloc(i,j)-Q0) &
                            /(1.+River%dt_slow*dQ_dV)
                else if (algor.eq.'nonlin') then
                  ! nonlin algorithm assumes all inflow at start of step 
                  ! and then integrates analytically
                  if (avail .gt. 0.) then
                      River%storage(i,j) = (avail**(1.-River%o_exp(i,j)) &
                        + River%o_coef(i,j)*(River%o_exp(i,j)-1.)*River%dt_slow) &
                                                  **(1./(1.-River%o_exp(i,j)))
                    else
                      River%storage(i,j) = avail
                    endif
                endif
! determine total water outflow during step
              River%outflow(i,j) = (avail - River%storage(i,j)) / River%dt_slow
! given outflow, determine flow width, depth, velocity
              if (River%outflow(i,j) .le. 0.) then
                      River%depth(i,j) = 0.
!                      River%width(i,j) = 0.
!                      River%vel(i,j)   = 0.
                else
                      River%depth(i,j) = River%d_coef(i,j) &
                       * River%outflow(i,j)**River%d_exp(i,j)
!                      River%width(i,j) = DHG_coef%on_w                        &
!                       *(River%outflowmean(i,j)**(DHG_exp%on_w-AAS_exp%on_w)) &
!                       *(River%outflow(i,j)**AAS_exp%on_w)
!                      River%vel(i,j) = River%outflow(i,j) /                   &
!                                        (River%width(i,j) * River%depth(i,j))
                endif
! given water outflow and storage, apportion other tracked stuff in same ratio
              out_frac = 0.
              if (avail .gt. 0.) out_frac = River%outflow(i,j)/avail
              River%outflow_c(i,j,:) = out_frac * (River%storage_c(i,j,:) &
                 +(River%inflow_c(i,j,:)+River%infloc_c(i,j,:))*River%dt_slow)
! add outflows to the inputs for downstream cells
              idest=i+di(River%tocell(i,j))
              jdest=j+dj(River%tocell(i,j))
              if (idest.eq.0           ) idest=River%nlon
              if (idest.eq.River%nlon+1) idest=1
              if (jdest.eq.0 .or. jdest.eq.River%nlat+1) then
                  write (*,*) 'WARNING: JDEST=', jdest
                  endif
              River%inflow(idest,jdest) = &
                  River%inflow(idest,jdest) + River%outflow(i,j) 
              River%inflow_c(idest,jdest,:) = &
                  River%inflow_c(idest,jdest,:) + River%outflow_c(i,j,:) 
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
!                  if (River%storage_c(i,j,1).gt.0. .and. conc(2).gt.tfreeze) then
!                      qmelt = min(DENS_H2O*hlf*River%storage_c(i,j,1), River%storage_c(i,j,2))
!                      conc(2) = conc(2) - qmelt / ( clw*River%storage(i,j) &
!                                                       + (csw-clw)*River%storage_c(i,j,1) )
!                      River%storage_c(i,j,1) = River%storage_c(i,j,1) - qmelt/hlf
!                      River%storage_c(i,j,2) = (DENS_H2O * &
!                               ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1) )) &
!                                * (conc(2) - tfreeze)
!                    endif
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

  end subroutine river_physics_step

  !#####################################################################

end module river_physics_mod
