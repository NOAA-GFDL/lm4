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
  character(len=128) :: version = '$Id: river_physics.F90,v 15.0 2007/08/14 03:59:36 fms Exp $'
  character(len=128) :: tagname = '$Name: omsk $'


  ! ---- public interfaces -----------------------------------------------------

  public :: river_physics_init, river_physics_step

  ! ---- namelist interface
  real, dimension(3) :: ave_DHG_exp = (/0.49,0.33,0.18/)  ! (/B, F, M for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_AAS_exp = (/0.19,0.39,0.42/)  ! (/b, f, m for avg of many rivers, 15Nov05/)
  real, dimension(3) :: ave_DHG_coef = (/4.62,0.26,0.82/) ! (/A, C, K for avg of many rivers, 15Nov05/)
  !     Coefficients are based on equations in m3/s, m/s, and m. 
  character*6 :: algor = 'linear'

  namelist /river_physics_nml/ ave_DHG_exp, ave_AAS_exp, ave_DHG_coef, algor

  !----------------------------------------------------------------------
  real,   parameter  :: NaN_real = -1.e8
  real               :: clw = 4218.
  real               :: csw = 2106.
  real,    parameter :: sec_in_day = 86400.
  real,    parameter :: v_r = 0.1 / sec_in_day

  real :: dens_inv

contains

  !#######################################################################

  subroutine river_physics_init( )

    integer            :: unit, io_status, ierr, siz(4)
    integer            :: isc, iec, jsc, jec       ! domain decomposition of land grid
    character(len=128) :: filename

    !--- read namelist -------------------------------------------------
    unit = open_namelist_file()
    read  (unit, river_physics_nml,iostat=io_status)
    ierr = check_nml_error(io_status,'river_physics_nml')
    call close_file (unit)

    !--- write version and namelist info to logfile --------------------
    call write_version_number(version,tagname)
    write (stdlog(), river_physics_nml)  

    dens_inv = 1.0/DENS_H2O

  end subroutine river_physics_init

  !#####################################################################

  subroutine river_physics_step(River, runoff, runoff_s, runoff_h, &
                                        runoff_c,    petravel,cur_travel)

    type(river_type),     intent(inout) :: River
    real, dimension(:,:),    intent(in) :: runoff, runoff_s, runoff_h
    real, dimension(:,:,:),    intent(in) :: runoff_c
    integer, dimension(:,:), intent(in) :: petravel
    integer,                 intent(in) :: cur_travel

    ! ---- local vars ----------------------------------------------------------
    integer   :: i, j, i_species, num_species
    real      :: Q0, dQ_dV, avail, out_frac, aaa, bbb, qmelt
    real, dimension (size(runoff,1),size(runoff,2)) :: inflow, inflow_s, inflow_h
    real, dimension (size(runoff,1),size(runoff,2),size(runoff_c,3)) :: inflow_c
    type(Leo_Mad_trios)   :: DHG_exp            ! downstream equation exponents
    type(Leo_Mad_trios)   :: DHG_coef           ! downstream equation coefficients
    type(Leo_Mad_trios)   :: AAS_exp            ! at-a-station equation exponents 

    call get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)
    bbb = 1./ (AAS_exp%on_w + AAS_exp%on_d)
    num_species = size(runoff_c,3)
    
! gather neighbors' outflows to form lateral inflow  
    call gather(River, River%outflow,   inflow,   petravel, cur_travel)
    call gather(River, River%outflow_s, inflow_s, petravel, cur_travel)
    call gather(River, River%outflow_h, inflow_h, petravel, cur_travel)
    do i_species = 1, num_species
      call gather(River, River%outflow_c(:,:,i_species), &
                        inflow_c(:,:,i_species), petravel, cur_travel)
      enddo

! do for all cells at current number of steps from river mouth
    do j = River%js(cur_travel),River%je(cur_travel) 
       do i = River%is(cur_travel),River%ie(cur_travel)
          if (petravel(i,j)==cur_travel) then

              aaa = River%outflowmean(i,j) / &
               (River%celllength(i,j)*DHG_coef%on_w*DHG_coef%on_d &
                 *River%outflowmean(i,j)**(DHG_exp%on_w+DHG_exp%on_d))**bbb
! form River copy of local source
              River%runoff  (i,j) = runoff  (i,j)*River%cellarea(i,j)*dens_inv
              River%runoff_s(i,j) = runoff_s(i,j)*River%cellarea(i,j)*dens_inv
              River%runoff_h(i,j) = runoff_h(i,j)*River%cellarea(i,j)*dens_inv
              do i_species = 1, num_species
                River%runoff_c(i,j,i_species) = &
                          runoff_c(i,j,i_species)*River%cellarea(i,j)*dens_inv
                enddo
! add local source to inflow
              inflow  (i,j) = inflow  (i,j) + River%runoff  (i,j)
              inflow_s(i,j) = inflow_s(i,j) + River%runoff_s(i,j)
              inflow_h(i,j) = inflow_h(i,j) + River%runoff_h(i,j)
              do i_species = 1, num_species
                inflow_c(i,j,i_species) = &
                              inflow_c(i,j,i_species) + River%runoff_c(i,j,i_species)
                enddo
! avail is amount of water to be partitioned between outflow and new storage
              avail = River%storage(i,j)+inflow(i,j)*River%dt_slow
! determine total water storage at end of step
              if (algor.eq.'linear') then
                  ! linear algorithm assumes outflow = Q0+dQ_dV*dS
                  if (River%storage(i,j) .le. 0.) then
                      Q0 = 0.; dQ_dV = 0.
                    else
                      Q0=aaa*River%storage(i,j)**bbb
                      dQ_dV=bbb*Q0/River%storage(i,j)
                    endif
                  River%storage(i,j) = River%storage(i,j) + River%dt_slow *   &
                                       (inflow(i,j)-Q0)/(1.+River%dt_slow*dQ_dV)
                else if (algor.eq.'nonlin') then
                  ! nonlin algorithm assumes all inflow at start of step 
                  ! and then integrates analytically
                  if (avail .gt. 0.) then
                      River%storage(i,j) =  &
                        (avail**(1.-bbb) + aaa*(bbb-1.)*River%dt_slow) &
                                                            **(1./(1.-bbb))
                    else
                      River%storage(i,j) = avail
                    endif
                endif
! determine total water outflow during step
              River%outflow(i,j) = (avail - River%storage(i,j)) / River%dt_slow
! given outflow, determine flow width, depth, velocity
              if (River%outflow(i,j) .le. 0.) then
                      River%width(i,j) = 0.
                      River%depth(i,j) = 0.
                      River%vel(i,j)   = 0.
                else
                      River%width(i,j) = DHG_coef%on_w                        &
                       *(River%outflowmean(i,j)**(DHG_exp%on_w-AAS_exp%on_w)) &
                       *(River%outflow(i,j)**AAS_exp%on_w)
                      River%depth(i,j) = DHG_coef%on_d                        &
                       *(River%outflowmean(i,j)**(DHG_exp%on_d-AAS_exp%on_d)) &
                       *(River%outflow(i,j)**AAS_exp%on_d)
                      River%vel(i,j) = River%outflow(i,j) /                   &
                                        (River%width(i,j) * River%depth(i,j))
                endif
! given water outflow and storage, apportion other tracked stuff in same ratio
              out_frac = 0.
              if (avail .gt. 0.) out_frac = River%outflow(i,j)/avail
              River%outflow_s(i,j) = out_frac * &
                          (River%storage_s(i,j)+inflow_s(i,j)*River%dt_slow)
              River%outflow_h(i,j) = out_frac * &
                          (River%storage_h(i,j)+inflow_h(i,j)*River%dt_slow)
              do i_species = 1, num_species
                River%outflow_c(i,j,i_species) = out_frac * &
                          (River%storage_c(i,j,i_species) &
                                     +inflow_c(i,j,i_species)*River%dt_slow)
                enddo
              River%storage_s(i,j) = River%storage_s(i,j)    &
                   + (inflow_s(i,j) - River%outflow_s(i,j)) * River%dt_slow
              River%storage_h(i,j) = River%storage_h(i,j)    &
                   + (inflow_h(i,j) - River%outflow_h(i,j)) * River%dt_slow
              do i_species = 1, num_species
                River%storage_c(i,j,i_species) = River%storage_c(i,j,i_species)&
                   + (inflow_c(i,j,i_species) - River%outflow_c(i,j,i_species))&
                                                            * River%dt_slow
                enddo
! define intensive variables for diagnostics and for use in transformations.
! along the way, melt swept snow as necessary. freeze will be a separate
! process, added later; it will be different in that frozen river water will
! be stationary, thus a different species
              if (River%storage(i,j) .gt. 0.) then
                  River%wtr_temp(i,j) = tfreeze + River%storage_h(i,j) /  &
                     ( clw*River%storage(i,j) + (csw-clw)*River%storage_s(i,j) )
!                  if (River%storage_s(i,j).gt.0. .and. River%wtr_temp(i,j).gt.tfreeze) then
!                      qmelt = min(hlf*River%storage_s(i,j), River%storage_h(i,j))
!                      River%wtr_temp(i,j) = River%wtr_temp(i,j) - qmelt / ( clw*River%storage(i,j) &
!                                + (csw-clw)*River%storage_s(i,j) )
!                      River%storage_s(i,j) = River%storage_s(i,j) - qmelt/hlf
!                      River%storage_h(i,j) = &
!                               ( clw*River%storage(i,j) + (csw-clw)*River%storage_s(i,j) ) &
!                                * (River%wtr_temp(i,j) - tfreeze)
!                    endif
                  River%snow_frac(i,j) = River%storage_s(i,j)/River%storage(i,j)
                else
                  River%wtr_temp(i,j)  = NaN_real
                  River%snow_frac(i,j) = NaN_real
                endif
! define intensive variables for diagnostics and for use in transformations.
! along the way, melt swept snow as necessary. freeze will be a separate
! process, added later; it will be different in that frozen river water will
! be stationary, thus a different species
              if (River%storage(i,j) .gt. 0.) then
                  River%conc(i,j,2) = tfreeze + River%storage_c(i,j,2) /  &
                    ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1))
!                  if (River%storage_c(i,j,1).gt.0. .and. River%conc(i,j,2).gt.tfreeze) then
!                      qmelt = min(DENS_H2O*hlf*River%storage_c(i,j,1), River%storage_c(i,j,2))
!                      River%conc(i,j,2) = River%conc(i,j,2) - qmelt / ( clw*River%storage(i,j) &
!                                                       + (csw-clw)*River%storage_c(i,j,1) )
!                      River%storage_c(i,j,1) = River%storage_c(i,j,1) - qmelt/hlf
!                      River%storage_c(i,j,2) = (DENS_H2O * &
!                               ( clw*River%storage(i,j) + (csw-clw)*River%storage_c(i,j,1) )) &
!                                * (River%conc(i,j,2) - tfreeze)
!                    endif
                  River%conc(i,j,1) = River%storage_c(i,j,1)/River%storage(i,j)
                else
                  River%conc(i,j,2) = NaN_real
                  River%conc(i,j,1) = NaN_real
                endif
              River%storage_c(i,j,3) = River%storage_c(i,j,3) &
                                         + River%storage(i,j)*River%dt_slow/sec_in_day
              River%storage_c(i,j,4) = River%storage_c(i,j,4) &
                                         * exp( -v_r * River%dt_slow/River%depth(i,j))
              if (River%storage(i,j) .gt. 0.) then
                  River%conc(i,j,3:5) = River%storage_c(i,j,3:5)/River%storage(i,j)
                else
                  River%conc(i,j,3:5) = NaN_real
                endif

          endif
       enddo
    enddo

  end subroutine river_physics_step

  !#####################################################################

  subroutine gather(River, source, dest, petravel, cur_travel)
    type(river_type),      intent(inout) :: River
    real,    dimension(0:River%nlon+1,0:River%nlat+1), intent(in)  :: source
    real,    dimension(:,:), intent(out) :: dest
    integer, dimension(:,:), intent(in)  :: petravel
    integer,                 intent(in)  :: cur_travel
    real                                 :: totinflow
    integer                              :: i, j
    do j = River%js(cur_travel),River%je(cur_travel)
       do i = River%is(cur_travel),River%ie(cur_travel)
          if (petravel(i,j)==cur_travel) then
             totinflow = 0.0
             if(River%fromcell_coef(i,j,1) == 1) totinflow = totinflow + source(i+1,j  )
             if(River%fromcell_coef(i,j,2) == 1) totinflow = totinflow + source(i+1,j-1  )
             if(River%fromcell_coef(i,j,3) == 1) totinflow = totinflow + source(i,  j-1  )
             if(River%fromcell_coef(i,j,4) == 1) totinflow = totinflow + source(i-1,j-1  )
             if(River%fromcell_coef(i,j,5) == 1) totinflow = totinflow + source(i-1,j  )
             if(River%fromcell_coef(i,j,6) == 1) totinflow = totinflow + source(i-1,j+1  )
             if(River%fromcell_coef(i,j,7) == 1) totinflow = totinflow + source(i,  j+1  )
             if(River%fromcell_coef(i,j,8) == 1) totinflow = totinflow + source(i+1,j+1  )
             dest(i,j) = totinflow
          endif
       enddo
    enddo
  end subroutine gather
  !#####################################################################

  subroutine get_Leo_Mad_params(DHG_exp, DHG_coef, AAS_exp)

    type(Leo_Mad_trios), intent(inout) :: DHG_exp  ! Exponents for downstream equations
    type(Leo_Mad_trios), intent(inout) :: DHG_coef ! Coefficients for downstream equations
    type(Leo_Mad_trios), intent(inout) :: AAS_exp  ! Exponents for at-a-station equations

!!! Exponents for the downstream hydraulic geometry equations
    DHG_exp%on_w = ave_DHG_exp(1) 
    DHG_exp%on_d = ave_DHG_exp(2)
    DHG_exp%on_V = ave_DHG_exp(3)

!!! Coefficients for the downstream hydraulic geometry equations
    DHG_coef%on_w = ave_DHG_coef(1)
    DHG_coef%on_d = ave_DHG_coef(2)
    DHG_coef%on_V = ave_DHG_coef(3)

!!! Exponents for the at-a-station hydraulic geometry equations
    AAS_exp%on_w = ave_AAS_exp(1)
    AAS_exp%on_d = ave_AAS_exp(2)
    AAS_exp%on_V = ave_AAS_exp(3)

  end subroutine get_Leo_Mad_params

  !#####################################################################

end module river_physics_mod
!             C = River%cel(i,j) * River%dt_slow / River%celllength(i,j)
!             if (use_fixes) then
!                if (C .GE. 2.0) then
!                    C = 2.0
!                    River%cel(i,j) = C*River%celllength(i,j)/River%dt_slow
!                  endif
!             endif
!
!             D = 2*River%diff(i,j)/(River%cel(i,j)*River%celllength(i,j)) 
!             if (use_fixes) then
!                if (D .LT. 1.0 - C) then
!                    D = 1.-C
!                    River%diff(i,j) = 0.5*D*River%cel(i,j)*River%celllength(i,j)
!                  else if (D .gt. 1+C) then
!                    D = 1.+C
!                    River%diff(i,j) = 0.5*D*River%cel(i,j)*River%celllength(i,j)
!                  endif
!               endif
!
!             totalinv = 1.0/(1+C+D)
!             C0 = (-1+C+D)*totalinv
!             C1 = (1+C-D)*totalinv
!             C2 = (1-C+D)*totalinv
!             River%outflow(i,j) = C0*River%inflow(i,j)+     &
!                                  C1*River%inflowlast(i,j)+ &
!                                  C2*River%outflowlast(i,j)
