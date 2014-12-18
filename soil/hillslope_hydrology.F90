! ============================================================================
! hillslope_hydrology_module
! ============================================================================
module hillslope_hydrology_mod

#include "../shared/debug.inc"

use fms_mod, only : write_version_number, error_mesg, FATAL, NOTE
use soil_tile_mod, only : &
     soil_tile_type, clw, initval, soil_data_hydraulic_properties
use land_tile_mod, only : land_tile_type, land_tile_enum_type, &
     first_elmt, tail_elmt, next_elmt, current_tile, operator(/=), nitems
use fms_mod, only : write_version_number
use land_data_mod,      only : land_state_type, lnd
use land_debug_mod, only : is_watch_point, set_current_point, get_current_point, &
                           check_conservation, is_watch_cell
use hillslope_mod, only : do_hillslope_model, strm_depth_penetration, use_hlsp_aspect_in_gwflow, &
                          use_geohydrodata, stiff_do_explicit, dammed_strm_bc, simple_inundation, &
                          surf_flow_velocity, limit_intertile_flow, flow_ratio_limit, exp_inundation
use constants_mod,      only: tfreeze, &
                              dens_h2o, epsln
      ! Use global tfreeze in energy flux calculations, not local freezing-point-depression temperature.
use fms_mod, only: error_mesg, FATAL
use time_manager_mod, only : time_type, time_type_to_real
use land_tile_diag_mod, only : diag_buff_type, register_tiled_diag_field, &
     send_tile_data

implicit none
private

! ==== public interfaces =====================================================
public :: hlsp_hydrology_1
public :: hlsp_hydrology_2
public :: hlsp_hydro_lev_init
public :: hlsp_hydro_init
public :: stiff_explicit_gwupdate
! =====end of public interfaces ==============================================


! ==== module constants ======================================================
character(len=*), parameter, private   :: &
    module_name = 'hillslope_hydrology',&
    version     = '$Id: hillslope_hydrology.F90,v 21.0 2014/12/15 21:51:13 fms Exp $',&
    tagname     = '$Name: ulm $'

! ---- diagnostic field IDs
!integer :: id_gtos_hlsp, & ! ground to stream runoff for hillslope area
!           id_gtosh_hlsp   ! ground to stream runoff heat for hillslope
!                           ! normalized to hillslope area (mm/s), (W/m^2), respectively
integer :: id_gdiv, & !  groundwater divergence (excl. to stream) (mm/s)
           id_ghdiv   !  heat flux associated with groundwater divergence (excl. to stream) (W/m^2)
integer :: id_gtos, & !  groundwater divergence from tile to stream (mm/s)
           id_gtosh   !  heat flux associated with groundwater divergence to stream (W/m^2)

! ==== module variables ======================================================
logical :: module_is_initialized =.FALSE.
integer :: num_l = -1 ! # of water layers
real, allocatable    :: dz    (:)   ! thicknesses of layers
real, allocatable    :: zfull (:)

contains

! ============================================================================
subroutine hlsp_hydro_lev_init(num_l_in, dz_in, zfull_in)
  integer, intent(in) :: num_l_in ! # of layers
  real   , intent(in) :: &
       dz_in(:), &  ! layer thickness
       zfull_in(:)  ! layer centers

  call write_version_number(version, tagname)
  module_is_initialized =.TRUE.

  num_l = num_l_in
  allocate(dz(num_l), zfull(num_l))
  dz(1:num_l)    = dz_in(1:num_l)
  zfull(1:num_l) = zfull_in(1:num_l)
end subroutine hlsp_hydro_lev_init

! Called from soil_step_2
! Retrieves inter-tile fluxes calculated in hlsp_hydrology_1, and calculates
! saturated fraction, storage, water table, and surface water diagnostics.
! Currently, saturated fraction is 0 or 1, based on water table depth previously
! calculated in soil_step_2. Storage_frac is either kept as the integrated-hillslope
! value if use_geohydrodata, or set based on zwt.  Surface_water is likewise the
! amount that the top soil layer exceeds saturation.  These will be more sophisticated when
! microtopography is implemented.
! ============================================================================
subroutine hlsp_hydrology_2(soil, psi, vlc, vsc, div_it, hdiv_it, &
              sat_frac, inun_frac, storage_frac, zwt, surface_water, runoff_frac)

   type(soil_tile_type), intent(inout) :: soil ! soil tile pointer
   real, intent(in)    :: psi(num_l) ! soil moisture potential wrt local elevation [m]
   real, intent(in)    :: vlc(num_l) ! volumetric liquid water content [-]
   real, intent(in)    :: vsc(num_l) ! volumetric solid water content [-]
   real, intent(out)   :: sat_frac   ! saturated area fraction used for sat runoff [-]
   real, intent(out)   :: inun_frac  ! diagnostic inundated area fraction [-]
   real, intent(inout) :: storage_frac ! diagnostic fraction of groundwater storage above stream elev [-]
   real, intent(in)    :: zwt ! depth to water table [m]; first saturated layer from top
   real, intent(out)   :: surface_water ! [m] surface water storage
   real, intent(out)   :: div_it(num_l) ! [mm/s] divergence of water due to inter-tile flow
                                        ! (incl. to stream)
   real, intent(out)   :: hdiv_it(num_l) ! divergence of heat due to inter-tile water flow [W/m^2]
   real, intent(out)   :: runoff_frac ! [-] effective saturated fraction used for calculating surface runoff, used
                        ! when simple inundation scheme is applied
   real :: delta_time ! [s]
   integer :: l


   delta_time = time_type_to_real(lnd%dt_fast)


   ! Retrieve water and energy divergences.
   div_it(1:num_l) = soil%div_hlsp(1:num_l)
   hdiv_it(1:num_l) = soil%div_hlsp_heat(1:num_l)


   ! Diagnostics
   ! Trivial for now
   surface_water = max(0., (vlc(1) + vsc(1) - soil%pars%vwc_sat)*dz(1))

   ! Set saturated fraction
   if (surface_water > 0.) then
      sat_frac = 1.
   else
      sat_frac = 0.
   end if
   
   ! Inundation and runoff_frac
   if (.not. simple_inundation .and. .not. exp_inundation) then
      if (psi(1) >= -soil%pars%microtopo) then
         inun_frac = min((soil%pars%microtopo + psi(1))/soil%pars%microtopo, 1.)
      else
         inun_frac = 0.
      end if
      runoff_frac = sat_frac
   else ! simple inundation scheme
        ! if exp_inundation then only diagnostic
      inun_frac = min(exp(-zwt / soil%pars%microtopo),1.)
      runoff_frac = inun_frac*inun_frac * soil%pars%tile_hlsp_slope * &
                    surf_flow_velocity / soil%pars%tile_hlsp_length * delta_time
                    ! extra inun_frac factor for channel connectivity
                    !       [-]               *         [-]         *
                    !       [m/s]      /                [m]         *    [s]
      runoff_frac = 1. - exp(-runoff_frac) ! To smooth behavior approaching 1.
      ! ZMS Could set surface_water here to soil%pars%microtopo*inun_frac,
      ! the simple integral, unless the above diagnostic proves useful.
   end if

   ! Storage above stream elevation
   if (.not. use_geohydrodata) then ! else leave storage_frac as is for comparison
      storage_frac = max((soil%pars%tile_hlsp_elev - zwt) / soil%pars%tile_hlsp_elev, 0.)
   end if

   if (is_watch_point()) then
      write(*,*)'hlsp_hydrology2,sat_frac:', sat_frac
      write(*,*)'psi(1:num_l):', psi(1:num_l)
      write(*,*)'vlc(1:num_l):', vlc(1:num_l)
      write(*,*)'storage_frac:', storage_frac
      write(*,*)'depth_to_saturation:', zwt
      write(*,*)'surface_water:', surface_water
!      write(*,*)'div(1:num_l):', div(1:num_l)
   end if

end subroutine hlsp_hydrology_2


! Calculate fluxes of water and associated heat between tiles in each gridcell.
! Called from update_land_model_fast. Occurs outside main tile loop. 
! ============================================================================
subroutine hlsp_hydrology_1(ground_to_stream, ground_to_stream_heat, ground_to_stream_tracers, &
                        num_species)
   ! Arguments
   integer, intent(in)  :: num_species ! number of tracer species

   real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je), intent(out) :: &
       ground_to_stream,  &  ! groundwater runoff directly to stream (mm/s)
       ground_to_stream_heat ! groundwater runoff heat directly to stream (W/m^2)

   real, dimension(lnd%is:lnd%ie,lnd%js:lnd%je,num_species), intent(out) :: &
       ground_to_stream_tracers ! groundwater runoff tracers directly to stream (1/m/s)

   integer ::     j,i,l,k,s
   type(land_tile_type), pointer :: tile, tile2 ! pointers to tile list elements
   type(land_tile_enum_type)     :: te, ce, ce2  ! tail and current tile list elements
   type(soil_tile_type), pointer :: soil, soil2 ! pointers to soil tiles

   ! Tile fractional area sums, needed for normalizing flux into tile.
   real    ::     area_above  ! running sum of contributing tile area (fraction)
   real    ::     area_level  ! running sum of tile area at same hillslope level (fraction)
   real    ::     area_below  ! running sum of tile area downslope (fraction)
   ! Gridcell fractional area sums
!   real    ::     area_stream ! running sum of tile area contributing to stream (fraction)

   ! Fluxes: the following fluxes are all + OUT of tile.
   real    ::     div_above(1:num_l)   ! running sum of contributing tile water flux (mm/s)
   real    ::     div_level(1:num_l)   ! running sum of tile water flux from same hillslope level (mm/s)
   real    ::     div_below(1:num_l)   ! running sum of tile water flux downslope (mm/s)
   real    ::     hdiv_above(1:num_l)  ! running sum of contributing tile energy flux (W/m^2)
   real    ::     hdiv_level(1:num_l)  ! running sum of tile energy flux from same hillslope level (W/m^2)
   real    ::     hdiv_below(1:num_l)  ! running sum of tile energy flux downslope  (W/m^2) (+ out of tile)
   real    ::     tdiv_above(1:num_l,1:num_species) ! running sum of contributing tile tracer flux (1/m^2/s)
   real    ::     tdiv_level(1:num_l,1:num_species) ! running sum of tile tracer flux from same hillslope
                                                    ! level (1/m^2/s)
   real    ::     tdiv_below(1:num_l,1:num_species) ! running sum of tile tracer flux downslope (1/m^2/s)
   real    ::     wbal        ! total water for balance check (mm/s)
   real    ::     ebal        ! total energy for balance check (W/m^2)
   real    ::     tbal(1:num_species)  ! total tracer for balance check (1/m^2/s)
   real    ::     k_hat       ! effective hydraulic conductivity between tiles for layer (mm/s)
   real    ::     ks          ! factor for limiting penetration of stream (-)
   real    ::     L_hat       ! mean hillslope length (m), possibly corrected for aspect
   real    ::     L1, L2      ! hillslope lengths 1 and 2 (m)
   real    ::     w1, w2      ! hillslope normalized widths 1 and 2 (-)
   real    ::     A1, A2      ! tile area fractions 1 and 2
   real    ::     w_hat       ! mean hillslope width (m)
   real    ::     deltapsi    ! Absolute difference in hydraulic head for tile 1 - tile 2 (m)
   real    ::     wflux       ! water flux, temporary (mm/s)
   real    ::     eflux       ! energy flux, temproary (W/m^2)
   real    ::     delta_h     ! elevation difference between tile 1 and 2 (m)
   real    ::     y           ! disturbance lengthscale (m)
   real, parameter :: wthresh = 1.e-14 ! water balance error threshold (mm/s)
   real, parameter :: ethresh = 1.e-8  ! energy balance error threshold (W/m^2)
   real    ::     frl         ! flow ratio limit: maximum ratio of head difference to length
                              ! allowed when limiting intertile flows

   ! For calling soil_data_hydraulic_properties
   real, dimension(num_l)  :: vlc, vsc, & ! volumetric fractions of water and ice in the layer (-)
                  DThDP, K_x, K_z, DKDP ! soil hydraulic parameters (not used)
   real :: DPsi_min, DPsi_max ! soil hydraulic parameters (not used)

!   ! For diagnostics
!   integer ::     max_hidx_k ! maximum hillslope parent index
!   real, allocatable, dimension(:)  :: gtos, gtosh ! ground to stream, ground to stream heat
!           ! for each hillslope, normalized per m^2 of hillslope area
!   real, allocatable, dimension(:)  :: hlsp_area ! total area of hillslope (-)
   real, allocatable :: gtos_bytile(:,:) ! [mm/s] gwater divergence to stream
                                                     ! dimension (tiles, num_l)
   real, allocatable :: gtosh_bytile(:,:) ! [W/m^2] gwater heat to stream
                                                     ! dimension (tiles, num_l)
   integer :: numtiles ! number of tiles in gridcell

   if (.not. do_hillslope_model) return

   ! Begin calculations
   
   frl = flow_ratio_limit

   ! Loop over gridcells

   do j=lnd%js,lnd%je
      do i=lnd%is,lnd%ie

         te = tail_elmt (lnd%tile_map(i,j))

         ! Initial loop over tile list
         ! ZMS for now this is an extra loop to calculate soil hydraulic props.
         ! This will need to be consolidated later.
         ce = first_elmt(lnd%tile_map(i,j))
         do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ce = next_elmt(ce)
            if (.not.associated(tile%soil)) cycle
            soil => tile%soil

            do l = 1,num_l
               vlc(l) = max(0., soil%wl(l) / (dens_h2o*dz(l)))
               vsc(l) = max(0., soil%ws(l) / (dens_h2o*dz(l)))
            end do

            call soil_data_hydraulic_properties ( soil, vlc, vsc, &
                 soil%psi, DThDP, K_x, K_z, DKDP, DPsi_min, DPsi_max)

            ! ZMS Find max hidx_k?
         
         end do


         ! Initialize for gridcell
         ground_to_stream(i,j) = 0.
         ground_to_stream_heat(i,j) = 0.
         ground_to_stream_tracers(i,j,:) = 0.
         ! ZMS Add diagnostic for hillslope-specific streamflow?
!         area_stream = 0.

         ! Get first pointer to tile list.
         ce = first_elmt(lnd%tile_map(i,j))
         k = 1

         ! Allocate and initialize gtos_bytile
         numtiles = nitems(lnd%tile_map(i,j))
         allocate(gtos_bytile(numtiles, num_l), gtosh_bytile(numtiles, num_l))
         gtos_bytile(:,:) = 0.
         gtosh_bytile(:,:) = 0.
         
         do while(ce /= te)
            tile=>current_tile(ce)  ! get pointer to current tile
            ! advance position to the next tile at bottom of loop, as ce needed in comparison
            if (.not.associated(tile%soil)) then
               ce = next_elmt(ce)
               k=k+1
               cycle
            end if

            soil => tile%soil

            ! Debug
            call set_current_point(i,j,k)
            if (is_watch_cell()) then
               write(*,*)'In hlsp_hydrology_1. In watch cell. At hidx_k, hidx_j:', soil%hidx_k, soil%hidx_j
               do l = 1,num_l
                  write(*,'(i3.3,99(x,a,g23.16))') l, 'soil%psi', soil%psi(l), 'soil%hyd_cond_horz', soil%hyd_cond_horz(l)
               enddo
            end if
         
            ! Initialize sums
            area_above = 0.
            area_level = 0.
            area_below = 0.
            div_above(:)  = 0.
            div_level(:)  = 0.
            div_below(:)  = 0.
            hdiv_above(:) = 0.
            hdiv_level(:) = 0.
            hdiv_below(:) = 0.
            tdiv_above(:,:) = 0.
            tdiv_level(:,:) = 0.
            tdiv_below(:,:) = 0.

            ! Get pointer to second tile list
            ce2 = first_elmt(lnd%tile_map(i,j))

            ! Loop over second tile list, and calculate fluxes 
            do while (ce2 /= te)
               tile2=>current_tile(ce2)
               if (.not.associated(tile2%soil)) then
                  ce2=next_elmt(ce2)
                  cycle
               end if

               soil2 => tile2%soil

               ! Check to see if in same hillslope
               if (soil%hidx_k == soil2%hidx_k) then
                  ! See if tile2 neighbors tile, above, same level, or below:

                  if (abs(soil%hidx_j - soil2%hidx_j) == 1) then ! tile2 is a vertical neighbor
                                                                 ! of tile

                     if (soil%hidx_j == soil2%hidx_j - 1) then ! tile2 above tile
                        area_above = area_above + tile2%frac
                     else ! tile2 below tile
                        area_below = area_below + tile2%frac
                     end if
                     L1 = soil%pars%tile_hlsp_length
                     L2 = soil2%pars%tile_hlsp_length
                     w1 = soil%pars%tile_hlsp_width
                     w2 = soil2%pars%tile_hlsp_width
                     L_hat = (L1 + L2)/2.
                     w_hat = (L2*w1 + L1*w2)/(L1 + L2)
                     delta_h = soil%pars%tile_hlsp_elev - soil2%pars%tile_hlsp_elev
                     if (use_hlsp_aspect_in_gwflow) then
                        L_hat = sqrt(L_hat*L_hat + delta_h*delta_h)
                     end if

                     ! Debug
                     if (is_watch_cell()) then
                        write(*,*)'Water & energy fluxes to tile with area, hk, hj:', &
                                  tile2%frac, soil2%hidx_k, soil2%hidx_j, '.'
                     end if

                     ! Loop over vertical layers
                     do l=1,num_l
                        ! Hydraulic conductivity is harmonic mean
                        ! Conductivity will be set to epsln if layer frozen
                        if (soil%hyd_cond_horz(l) == epsln .or. soil2%hyd_cond_horz(l) == epsln) then
                           k_hat = 0.
                        else
                           k_hat = (soil%hyd_cond_horz(l)*soil2%hyd_cond_horz(l)*(L1+L2)) / &
                                   (soil%hyd_cond_horz(l)*L2 + soil2%hyd_cond_horz(l)*L1)
                        end if
                        ! Total Soil moisture potential
                        deltapsi = soil%psi(l) - soil2%psi(l) + delta_h
                        if (limit_intertile_flow .and. abs(deltapsi) > frl*L_hat) then
                           if (deltapsi > 0.) then
                              deltapsi = frl*L_hat
                           else
                              deltapsi = -frl*L_hat
                           end if
                        end if

                        ! Flux from tile --> tile2
                        ! Will later be normalized by total area of tiles above or below
                        wflux = k_hat * deltapsi/L_hat * dz(l) * tile2%frac / L1 * w_hat / w1
                       ! mm/s =  mm/s *   m     / m    *  m          -      / m   * -    / -

                        ! Energy flux
                        if (wflux < 0.) then ! water flowing into tile: heat advected in
                           eflux = wflux * (soil2%T(l)- tfreeze)* clw
                        else                 ! water flowing out of tile: heat advected out
                           eflux = wflux * (soil%T(l)- tfreeze)* clw
                        end if

                        ! Update fluxes
                        if (soil%hidx_j == soil2%hidx_j - 1) then ! tile2 above tile
                           div_above(l) = div_above(l) + wflux
                           hdiv_above(l) = hdiv_above(l) + eflux
                        else ! tile2 below tile
                           div_below(l) = div_below(l) + wflux
                           hdiv_below(l) = hdiv_below(l) + eflux
                        end if

                        do s=1,num_species
                        !! Add code here
                        end do

                        ! Debug
                        if (is_watch_cell()) then
                           write(*,*)'at level l=',l,', wflux, eflux out of tile = ', &
                                 wflux, eflux, '.'
                        end if

                     end do ! l


                  else if (soil%hidx_j == soil2%hidx_j .and. ce /= ce2) then
                     ! tile2 is in same level
                     A1 = tile%frac
                     A2 = tile2%frac
                     area_level = area_level + A2
                     ! A1 will be added to area_level at bottom.
                     y = soil%pars%disturb_scale
                     
                     ! Debug
                     if (is_watch_cell()) then
                        write(*,*)'Water & energy fluxes to tile with area, hk, hj:', &
                                  tile2%frac, soil2%hidx_k, soil2%hidx_j, '.'
                     end if
                     
                     ! Loop over vertical layers
                     do l=1,num_l
                        ! Hydraulic conductivity is harmonic mean
                        ! Conductivity will be set to epsln if layer frozen
                        if (soil%hyd_cond_horz(l) == epsln .or. soil2%hyd_cond_horz(l) == epsln) then
                           k_hat = 0.
                        else
                           k_hat = (soil%hyd_cond_horz(l)*soil2%hyd_cond_horz(l)*(A1+A2)) / &
                                   (soil%hyd_cond_horz(l)*A2 + soil2%hyd_cond_horz(l)*A1)
                        end if
                        ! Total Soil moisture potential
                        deltapsi = soil%psi(l) - soil2%psi(l) ! No delta_h
                        if (limit_intertile_flow .and. abs(deltapsi) > frl*y) then
                           if (deltapsi > 0.) then
                              deltapsi = frl*y
                           else
                              deltapsi = -frl*y
                           end if
                        end if

                        ! Flux from tile --> tile2
                        ! Will be normalized below by area_level.
                        wflux = k_hat * deltapsi / (y*y)/2. * dz(l) * A2 ! ZMS double-check this
                       ! mm/s =  mm/s *   m      /   m^2      *m      --
                        div_level(l) = div_level(l) + wflux

                        ! Energy flux
                        if (wflux < 0.) then ! water flowing into tile: heat advected in
                           eflux = wflux * (soil2%T(l)- tfreeze) * clw
                        else                 ! water flowing out of tile: heat advected out
                           eflux = wflux * (soil%T(l)- tfreeze) * clw
                        end if
                        hdiv_level(l) = hdiv_level(l) + eflux

                        do s=1,num_species
                        !! Add code here
                        end do

                        ! Debug
                        if (is_watch_cell()) then
                           write(*,*)'at level l=',l,', wflux, eflux out of tile = ', &
                                 wflux, eflux, '.'
                        end if

                     end do ! l

                  end if ! hidx_j

               end if ! hidx_k

               ce2=next_elmt(ce2)
            end do ! loop over second tile list

            ! Initialize outputs
            soil%div_hlsp(:) = 0.
            soil%div_hlsp_heat(:) = 0.

            ! Add to outputs
            if (area_above > 0.) then
               soil%div_hlsp(:) = soil%div_hlsp(:) + div_above(:) / area_above
               soil%div_hlsp_heat(:) = soil%div_hlsp_heat(:) + hdiv_above(:) / area_above
               ! Add tracer code
            end if
            
            if (area_level > 0.) then
               ! Add current tile to area_level
               area_level = area_level + tile%frac
               soil%div_hlsp(:) = soil%div_hlsp(:) + div_level(:) / area_level
               soil%div_hlsp_heat(:) = soil%div_hlsp_heat(:) + hdiv_level(:) / area_level
               ! Add tracer code
            end if
            
            if (area_below > 0.) then
               soil%div_hlsp(:) = soil%div_hlsp(:) + div_below(:) / area_below
               soil%div_hlsp_heat(:) = soil%div_hlsp_heat(:) + hdiv_below(:) / area_below
               ! Add tracer code                 
            end if


            ! Debug
            if (is_watch_cell()) then
               write(*,*)'Total divergence of water, heat (1:num_l) = ', soil%div_hlsp(1:num_l), &
                           ', ', soil%div_hlsp_heat(1:num_l), '.'
            end if

            ! Calculate flows directly to stream
            if (soil%hidx_j == 1) then
            ! ZMS for now assume psi(stream) == 0.
               A1 = tile%frac
!               area_stream = area_stream + A1
               L1 = soil%pars%tile_hlsp_length
               L_hat = L1/2.
               delta_h = soil%pars%tile_hlsp_elev
               if (use_hlsp_aspect_in_gwflow) then
                  L_hat = sqrt(L_hat*L_hat + delta_h*delta_h)
               end if

               ! Loop over vertical layers
               do l=1,num_l
                  ! Hydraulic conductivity is harmonic mean
                  ! Conductivity will be set to epsln if layer frozen
                  if (soil%hyd_cond_horz(l) == epsln) then
                     k_hat = 0.
                  else
                     k_hat = soil%hyd_cond_horz(l)
                  end if
                  ! Total Soil moisture potential
                  deltapsi = soil%psi(l) + delta_h
                  if (limit_intertile_flow .and. abs(deltapsi) > frl*L_hat) then
                     if (deltapsi > 0.) then
                        deltapsi = frl*L_hat
                     else
                        deltapsi = -frl*L_hat
                     end if
                  end if

                  if (dammed_strm_bc) deltapsi = deltapsi - zfull(l)
                  ! ZMS: Subtracting zfull here makes this a "wet valley or dammed stream"
                  ! BC rather than a "steel wall" or open stream BC.

                  ks = exp(-zfull(l)/strm_depth_penetration) ! Assume smooth relationship with
                                                             ! depth

                  ! Flux from tile --> stream, per unit area of tile
                  wflux = k_hat*ks * deltapsi/L_hat * dz(l)/L1
                                                                        ! ZMS double-check this
                  ! mm/s =  mm/s *   m     / m    *  m  /m

                  ! Energy flux
                  if (wflux < 0.) then ! water flowing into tile: heat advected in
                     ! Need to access stream state: eflux = wflux * __ * clw
                     ! ZMS For now don't allow this.
                  else  ! water flowing into stream: heat advected out
                     eflux = wflux * (soil%T(l)- tfreeze)* clw
                  end if

                  do s=1,num_species
                  !! Add code here
                  end do


                  if (wflux < 0.) then
                     ! ZMS for now only allow flow into stream
                  else
                     soil%div_hlsp(l) = soil%div_hlsp(l) + wflux
                     soil%div_hlsp_heat(l) = soil%div_hlsp_heat(l) + eflux
                     ! For diagnostics
                     gtos_bytile(k,l) = wflux
                     gtosh_bytile(k,l) = eflux

                     ! In units flowing into stream
                     ground_to_stream(i,j) = ground_to_stream(i,j) + wflux * A1
                     ground_to_stream_heat(i,j) = ground_to_stream_heat(i,j) + eflux * A1

                     do s=1,num_species
                     !! Add code here
                     end do
                  end if

               end do ! l

               ! Debug
               if (is_watch_cell()) then
                  write(*,*)'Flows to stream for tile 1: wflux(1:num_l):', gtos_bytile(k,1:num_l)
                  write(*,*)'eflux(1:num_l):', gtosh_bytile(k,1:num_l)
               end if

            end if ! Stream


            ce=next_elmt(ce)
            k=k+1
         end do  ! loop over first tile list

         ! Normalize flows to stream
         ! ZMS this is not appropriate.  Replace by hillslope areas for diagnostics.
         !ground_to_stream(i,j) = ground_to_stream(i,j) / area_stream
         !ground_to_stream_heat(i,j) = ground_to_stream_heat(i,j) / area_stream
         !ground_to_stream_tracers(i,j,:) = ground_to_stream_tracers(i,j,:) / area_stream

         ! End of gridcell calculations
         
         ! Check conservation of water, energy, and tracers, 
         ! and send tile diagnostics.
         ! Repeat single loop over tile list
         ce = first_elmt(lnd%tile_map(i,j))
         
         wbal = 0.
         ebal = 0.
         tbal(:) = 0.

         k = 0
         do while(ce /= te)
            k = k + 1
            tile=>current_tile(ce)  ! get pointer to current tile
            ce = next_elmt(ce)
            if (.not.associated(tile%soil)) cycle

            do l=1,num_l
               wbal = wbal + tile%soil%div_hlsp(l) * tile%frac
               ebal = ebal + tile%soil%div_hlsp_heat(l) * tile%frac
               do s=1,num_species
               !! Add code here
               end do
            end do

            ! Send tile diagnostics
            call send_tile_data(id_gdiv, tile%soil%div_hlsp - gtos_bytile(k,:), &
                                         tile%diag)
            call send_tile_data(id_ghdiv, tile%soil%div_hlsp_heat - gtosh_bytile(k,:), &
                                         tile%diag)
            call send_tile_data(id_gtos, gtos_bytile(k,:), tile%diag)
            call send_tile_data(id_gtosh, gtosh_bytile(k,:), tile%diag)

         end do

         wbal = wbal - ground_to_stream(i,j)       ! + is into stream
         ebal = ebal - ground_to_stream_heat(i,j)  ! + is into stream

         !! Tracers add code

         call check_conservation('hlsp_hydrology_1, between-tile fluxes (k value following is not valid)', &
                                 'Water', wbal, 0., wthresh, lnd%time, FATAL)

         call check_conservation('hlsp_hydrology_1, between-tile energy fluxes (k value following is not valid)', &
                                 'Energy', ebal, 0., ethresh, lnd%time, FATAL)

         deallocate(gtos_bytile, gtosh_bytile)
 
      end do ! i
   end do ! j
         

end subroutine hlsp_hydrology_1

! Initialize diagnostic fields. 
! ============================================================================
subroutine hlsp_hydro_init (id_lon, id_lat, id_zfull)
   integer, intent(in) :: id_lon  ! ID of land longitude (X) axis  
   integer, intent(in) :: id_lat  ! ID of land latitude (Y) axis
   integer, intent(in) :: id_zfull ! ID of vertical soil axis

   ! ---- local vars
   integer :: axes(3)

   if (.not. do_hillslope_model) return

   ! define array of axis indices
   axes = (/ id_lon, id_lat, id_zfull /)

   ! define diagnostic fields

!   id_gtos_hlsp = register_tiled_diag_field ( module_name, 'grnd_to_stream_water', axes(1:2), &
!       Time, 'groundwater flux to stream at hillslope bottom, normalized by hillslope area', 'mm/s',  missing_value=initval )
!   id_gtosh_hlsp = register_tiled_diag_field ( module_name, 'grnd_to_stream_heat', axes(1:2), &
!       Time, 'advected groundwater heat flux to stream at hillslope bottom, normalized by hillslope area', 'W/m^2',  missing_value=initval )

   id_gdiv = register_tiled_diag_field ( module_name, 'groundwater_divergence', axes, &
       lnd%time, 'groundwater divergence out of tiles, excluding flow to stream (i.e., baseflow)', 'mm/s', &
       missing_value=initval )
   id_ghdiv = register_tiled_diag_field ( module_name, 'groundwater_heat_div', axes, &
       lnd%time, 'heat flux associated with groundwater divergence (excl. to stream)', 'W/m^2', missing_value=initval )
   id_gtos = register_tiled_diag_field ( module_name, 'grounddiv_to_stream', axes, &
       lnd%time, 'groundwater divergence out of tiles directly to stream', 'mm/s', &
       missing_value=initval )
   id_gtosh = register_tiled_diag_field ( module_name, 'groundheatdiv_to_stream', axes, &
       lnd%time, 'heat flux associated with groundwater divergence to stream', 'W/m^2', missing_value=initval )

end subroutine hlsp_hydro_init

! In the case where soil profile becomes "stiff", i.e. completely frozen, during timestep,
! must apply inter-tile water & energy tendencies explicitly to maintain conservation. 
! ============================================================================
subroutine stiff_explicit_gwupdate (soil, div_it, hdiv_it, div, lrunf_bf)
   type(soil_tile_type), intent(inout) :: soil
   real, intent(inout)   :: div_it(num_l) ! [mm/s] divergence of water due to inter-tile flow
                                        ! (incl. to stream)
   real, intent(inout)   :: hdiv_it(num_l) ! divergence of heat due to inter-tile water flow [W/m^2]
   real, intent(inout)   :: div(num_l) ! total divergence of soil water [mm/s]
   real, intent(inout)   :: lrunf_bf ! baseflow [mm/s] (including inter-tile flow)
   ! Debug
   character(len=512) :: mesg
   real    :: time
   integer :: i,j,k,t

   call get_current_point(i,j,k,t)
   time = time_type_to_real(lnd%time)
   write(mesg,*)'Tile at point (i,j,k,face): ', i,',', j, ',', k, ',', t, ' has become stiff', &
                ' (i.e., completely frozen) during the timestep at time ', time, '.', &
                'Explicitly updating for inter-tile ', &
                'exchanges of water and energy in order to maintain conservation.'

   ! For now always write this.

!   if (is_watch_point()) then
      call error_mesg(module_name, mesg, NOTE)
!   end if

   if (.not. stiff_do_explicit) then
      ! Zero out all runoff, and let energy and water be balanced in the stream.
      div_it(:) = 0.
      hdiv_it(:) = 0.
      div(:) = 0.
      lrunf_bf = 0.
!      if (is_watch_point()) then
         write(mesg,*)'Namelist option "stiff_do_explicit" selected, so runoff is being zeroed out, ', &
                   'and water & energy will be balanced in the stream.'
         call error_mesg(module_name, mesg, NOTE)
!      end if
      return
   end if

   ! ZMS Need to fill in function. For now stiff_do_explicit is false.
   ! Logic: check to see if total water sink exceed ice + water mass.
   ! If so, need to use logic above to avoid negative values of wl.
   ! Else, save initial water & ice content.
   ! apply water / ice sources & sinks; remove water preferentially and then ice from each layer.
   ! Keep track of implied energy sources and sinks.
   ! If water goes negative in a layer, move water / ice from bottom-most layer with some available,
   ! proceeding upwards until negatives are eliminated.
   ! Compare implied energy sources and sinks by comparing to initial water & ice content, with
   ! prescribed energy sources and sinks from inter-tile exchange.  Correct temperature and/or phase.
   

end subroutine stiff_explicit_gwupdate

end module hillslope_hydrology_mod
