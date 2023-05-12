module soilc_type_mod

use land_data_mod, only : lnd ! only for deplete_pool
use tile_diag_buff_mod, only : diag_buff_type
use soil_tile_mod, only: soil_tile_type
use vegn_tile_mod, only: vegn_tile_type

implicit none; private

public :: soilc_t
public :: deplete_pool

! abstract type representing soil carbon model
type, abstract :: soilc_t
contains
  procedure (merge),         deferred, pass :: merge   ! merge another soil carbon tile into current one
  procedure (get_real_func), deferred, pass :: total_C ! returns total C [kgC/m2]
  procedure (get_real_func), deferred, pass :: total_N ! returns total N [kgN/m2]
  procedure (get_real_3),    deferred, pass :: rav_C   ! returns amounts of C [kgC/m2]
                                                       ! for legacy surface resistance calculations
  procedure (get_real_2D),   deferred, pass :: get_DOC ! returns DOC, by type and by layer
  procedure (get_real_2D),   deferred, pass :: get_DON ! returns DON, by type and by layer
  procedure (get_real_1D),   deferred, pass :: get_nit ! returns nitrate by layer, kgN/m2
  procedure (get_real_1D),   deferred, pass :: get_amm ! returns ammonium by layer, kgN/m2

  procedure (add_soil_carbon),   deferred, pass :: add_soil_carbon ! add new surface and sub-surface litter to soil carbon and nitrogen
  procedure (add_root_litter),   deferred, pass :: add_root_litter ! add new root litter to soil carbon and nitrogen
  procedure (add_root_exudates), deferred, pass :: add_root_exudates ! add root exudates to soil carbon

  procedure (update_soil_pools), deferred, pass :: update_soil_pools
  procedure (dsdt),              deferred, pass :: dsdt
  procedure (step3),             deferred, pass :: step3
end type

! ---- abstract interfaces for methods
abstract interface
   ! merge sc1 into sc2 with given weights
   subroutine merge(s2,w2,s1,w1)
      import :: soilc_t
      class(soilc_t), intent(inout) :: s2
      class(soilc_t), intent(in)    :: s1
      real          , intent(in)    :: w2,w1 ! merging weights
   end subroutine merge

   ! given soil carbon data, returns real number
   function get_real_func(soilC)
      import :: soilc_t ! soil carbon data structure
      class(soilc_t), intent(in) :: soilC
   end function

   ! given soil carbon data, returns three kinds of carbon
   subroutine get_real_3(soilC, fast_C, slow_C, dmic_C)
      import :: soilc_t
      class(soilc_t), intent(in)  :: soilC ! soil carbon data structure
      real, intent(out) :: &
         fast_C,    & ! fast litter carbon, [kgC/m2]
         slow_C,    & ! slow litter carbon, [kgC/m2]
         dmic_C       ! mass of dead microbes in litter, [kgC/m2]
   end subroutine

   ! given soil carbon data, returns 2D data
   subroutine get_real_2D(soilC, values)
      import :: soilc_t ! soil carbon data structure
      class(soilc_t), intent(in) :: soilC
      real,           intent(out):: values(:,:) ! in many cases (N_C_TYPES, num_l)
   end subroutine

   ! given soil carbon data, returns 2D data
   subroutine get_real_1D(soilC, values)
      import :: soilc_t ! soil carbon data structure
      class(soilc_t), intent(in) :: soilC
      real,           intent(out):: values(:) ! (num_l)
   end subroutine

   subroutine add_soil_carbon(soilC, vegn, &
          leaf_litter_C, wood_litter_C, root_litter_C, &
          leaf_litter_N, wood_litter_N, root_litter_N  )
      import :: soilc_t,vegn_tile_type

      class(soilc_t),       intent(inout) :: soilC
      type(vegn_tile_type), intent(inout) :: vegn

      real, intent(in), optional :: leaf_litter_C(:)   ! (N_C_TYPES)
      real, intent(in), optional :: wood_litter_C(:)   ! (N_C_TYPES)
      real, intent(in), optional :: root_litter_C(:,:) ! (num_l,N_C_TYPES)
      real, intent(in), optional :: leaf_litter_N(:)   ! (N_C_TYPES)
      real, intent(in), optional :: wood_litter_N(:)   ! (N_C_TYPES)
      real, intent(in), optional :: root_litter_N(:,:) ! (num_l,N_C_TYPES)
   end subroutine

   subroutine add_root_litter(soilC, vegn, litterC, litterN)
      import :: soilc_t,vegn_tile_type
      class(soilc_t),       intent(inout)  :: soilC ! soil carbon data structure
      type(vegn_tile_type), intent(in)     :: vegn ! vegetation data structure, for rhizosphere caculations
      real, intent(in) :: litterC(:, :) ! (num_l, N_C_TYPES) kgC/m2 of soil layer
      real, intent(in) :: litterN(:, :) ! (num_l, N_C_TYPES) kgN/m2 of soil layer
   end subroutine

   ! add root exudates to soil carbon
   subroutine add_root_exudates(soilC, exudateC, exudateN, ammonium, nitrate)
      import :: soilc_t
      class(soilc_t), intent(inout)  :: soilC ! soil carbon data structure
      real,intent(in)           :: exudateC(:) ! (num_l) amount of C in exudate, kgC/m2 per layer
      real,intent(in), optional :: exudateN(:) ! (num_l) amount of N in exudate, kgN/m2 per layer
      real,intent(in), optional :: ammonium(:) ! (num_l) amount of ammonium in exudate, kgN/m2(?) per layer
      real,intent(in), optional :: nitrate (:) ! (num_l) amount of  nitrate in exudate, kgN/m2(?) per layer
   end subroutine

   subroutine update_soil_pools(soilc, vegn)
      import :: soilc_t, vegn_tile_type
      class(soilc_t),       intent(inout) :: soilc
      type(vegn_tile_type), intent(inout) :: vegn
   end subroutine

   subroutine dsdt(soilc, soil, vegn, diag, soilt, theta)
      import soilc_t, soil_tile_type, vegn_tile_type, diag_buff_type
      class(soilc_t), intent(inout)       :: soilc
      type(soil_tile_type), intent(inout) :: soil
      type(vegn_tile_type), intent(inout) :: vegn
      type(diag_buff_type), intent(inout) :: diag
      real                , intent(in)    :: soilt ! average soil temperature, deg K
      real                , intent(in)    :: theta ! average soil moisture
   end subroutine

   subroutine step3(soilc, diag)
      import soilc_t, diag_buff_type
      class(soilc_t),       intent(inout) :: soilc
      type(diag_buff_type), intent(inout) :: diag
   end subroutine
end interface

contains

!> @brief Move substance from one pool to another
!!
!! given an intermediate pool of C or N, and its spending rate, move the amount
!! of mass corresponding to one fast time step from the pool to the destination.
!! The spending rate is adjusted so that intermediate pool is never depleted below zero.
!!
!! @note
!! Spending rate argument is also updated, to be correctly reported to diagnostics
subroutine deplete_pool(pool, rate, dest, accum)
   real, intent(inout) :: pool !< C or N intermediate pool, kg
   real, intent(inout) :: rate !< C or N spending rate, kg/yr
   real, intent(inout) :: dest !< C or N destination pool, kg
   real, intent(inout), optional :: accum !< accumulator for soil carbon equilibration, e.g. fs_in or ssc_in

   real :: delta ! change in pool over time step, kg

   rate  = MAX( 0.0, MIN(rate, pool/lnd%dt_fast_yr) ) ! adjust rate
   delta = rate * lnd%dt_fast_yr
   dest  = dest + delta
   pool  = pool - delta
   if (present(accum)) accum = accum + delta ! increment accumulator
end subroutine deplete_pool

end module
