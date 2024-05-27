module soil_BGC_type_mod

use land_data_mod, only : lnd ! only for deplete_pool
use land_numerics_mod, only : tridiag
use tile_diag_buff_mod, only : diag_buff_type
use soil_tile_mod, only: soil_tile_type
use vegn_tile_mod, only: vegn_tile_type

implicit none; private

public :: soil_BGC_t
public :: deplete_pool
public :: tracer_advection

! abstract type representing soil carbon model
type, abstract :: soil_BGC_t
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
  procedure (get_real_1D),   deferred, pass :: get_littC ! returns litter carbon, by litter pool, kgC/m2

  procedure (add_soil_matter),   deferred, pass :: add_soil_matter   ! add new surface and sub-surface litter to soil carbon and nitrogen
  procedure (add_root_exudates), deferred, pass :: add_root_exudates ! add root exudates to soil carbon
  procedure (burn_litter_frac),  deferred, pass :: burn_litter_frac  ! burn a fraction of sfc litter and retuen amounts of burned carbon and nitrogen
  procedure (tracer_leaching),   deferred, pass :: tracer_leaching

  procedure (deposit_N),              deferred, pass :: deposit_N
  procedure (active_root_N_uptake),   deferred, pass :: active_root_N_uptake
  procedure (myc_miner_N_uptake),     deferred, pass :: myc_miner_N_uptake
  procedure (myc_scavenger_N_uptake), deferred, pass :: myc_scavenger_N_uptake

  procedure (spend_intermediate_pools), deferred, pass :: spend_intermediate_pools
  procedure (dsdt),              deferred, pass :: dsdt
  procedure (step3),             deferred, pass :: step3
  procedure (redistribute_peat_carbon), deferred, pass :: redistribute_peat_carbon
end type

! ---- abstract interfaces for methods
abstract interface
   ! merge sc1 into sc2 with given weights
   subroutine merge(s2,w2,s1,w1)
      import :: soil_BGC_t
      class(soil_BGC_t), intent(inout) :: s2
      class(soil_BGC_t), intent(in)    :: s1
      real          , intent(in)    :: w2,w1 ! merging weights
   end subroutine merge

   ! given soil carbon data, returns real number
   function get_real_func(soilC)
      import :: soil_BGC_t ! soil carbon data structure
      class(soil_BGC_t), intent(in) :: soilC
   end function

   ! given soil carbon data, returns three kinds of carbon
   subroutine get_real_3(soilC, fast_C, slow_C, dmic_C)
      import :: soil_BGC_t
      class(soil_BGC_t), intent(in)  :: soilC ! soil carbon data structure
      real, intent(out) :: &
         fast_C,    & ! fast litter carbon, [kgC/m2]
         slow_C,    & ! slow litter carbon, [kgC/m2]
         dmic_C       ! mass of dead microbes in litter, [kgC/m2]
   end subroutine

   ! given soil carbon data, returns 2D data
   subroutine get_real_2D(soilC, values)
      import :: soil_BGC_t ! soil carbon data structure
      class(soil_BGC_t), intent(in) :: soilC
      real,           intent(out):: values(:,:) ! in many cases (N_C_TYPES, num_l)
   end subroutine

   ! given soil carbon data, returns 2D data
   subroutine get_real_1D(soilC, values)
      import :: soil_BGC_t ! soil carbon data structure
      class(soil_BGC_t), intent(in) :: soilC
      real,           intent(out):: values(:) ! (num_l)
   end subroutine

   subroutine add_soil_matter(soilC, vegn, &
          leaf_litter_C, wood_litter_C, root_litter_C, &
          leaf_litter_N, wood_litter_N, root_litter_N  )
      import :: soil_BGC_t,vegn_tile_type

      class(soil_BGC_t),       intent(inout) :: soilC
      type(vegn_tile_type), intent(inout) :: vegn

      real, intent(in), optional :: leaf_litter_C(:)   ! (N_C_TYPES)
      real, intent(in), optional :: wood_litter_C(:)   ! (N_C_TYPES)
      real, intent(in), optional :: root_litter_C(:,:) ! (num_l,N_C_TYPES)
      real, intent(in), optional :: leaf_litter_N(:)   ! (N_C_TYPES)
      real, intent(in), optional :: wood_litter_N(:)   ! (N_C_TYPES)
      real, intent(in), optional :: root_litter_N(:,:) ! (num_l,N_C_TYPES)
   end subroutine

   ! add root exudates to soil carbon
   subroutine add_root_exudates(soilC, exudateC, exudateN, ammonium, nitrate)
      import :: soil_BGC_t
      class(soil_BGC_t), intent(inout)  :: soilC ! soil carbon data structure
      real,intent(in)           :: exudateC(:) ! (num_l) amount of C in exudate, kgC/m2 per layer
      real,intent(in), optional :: exudateN(:) ! (num_l) amount of N in exudate, kgN/m2 per layer
      real,intent(in), optional :: ammonium(:) ! (num_l) amount of ammonium in exudate, kgN/m2(?) per layer
      real,intent(in), optional :: nitrate (:) ! (num_l) amount of  nitrate in exudate, kgN/m2(?) per layer
   end subroutine

   subroutine burn_litter_frac(soilc, frac, burned_C, burned_N)
      import :: soil_BGC_t
      class(soil_BGC_t), intent(inout) :: soilc
      real, intent(in)  :: frac(:) ! (N_LITTER_POOLS) fraction of litter to burn [0,1], per litter pool
      real, intent(out) :: burned_C, burned_N ! amounts of burned carbon and nitrogen
   end subroutine

   subroutine tracer_leaching(soilC, diag, &
         wl, flow, div, &
         div_hlsp_DOC, div_hlsp_DON, div_hlsp_NO3, div_hlsp_NH4, &
         ! output
         total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div )
      import :: soil_BGC_t, diag_buff_type
      type(diag_buff_type), intent(inout) :: diag
      class(soil_BGC_t), intent(inout) :: soilC
      real, intent(in) :: flow(:), div(:), wl(:) ! flow (into layer) and wl in units of mm, downward is >0  !!!xz check the unit of dz (should be m in this subroutine), flow (shoul be mm)
      real, intent(in) :: div_hlsp_DOC(:,:) ! (N_C_TYPES, num_l) [kg C/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
      real, intent(in) :: div_hlsp_DON(:,:) ! (N_C_TYPES, num_l) [kg N/m^2/s] net divergence
      real, intent(in) :: div_hlsp_NO3(:),div_hlsp_NH4(:) ! (num_l) [kg N/m^2/s] net divergence loss from tile calculated in hlsp_hydrology
      real, intent(out) :: total_DOC_div, total_DON_div, total_NO3_div, total_NH4_div
   end subroutine tracer_leaching

   subroutine deposit_N(soilc, NH4, NO3, N_org)
      import :: soil_BGC_t
      class(soil_BGC_t), intent(inout) :: soilc
      real, intent(in) :: NH4, NO3, N_org ! amounts of NH4, NO3, and organic nitrogen to deposit, kg N/m2
   end subroutine deposit_N

   subroutine active_root_N_uptake(soilc, vegn, N_uptake, dt, update_pools)
      import :: soil_BGC_t, vegn_tile_type
      class(soil_BGC_t), intent(inout) :: soilc
      type(vegn_tile_type), intent(in)    :: vegn
      real,    intent(out) :: N_uptake(:) ! Nitrogen uptake, kg N per individual
      real,    intent(in)  :: dt ! in years
      logical, intent(in)  :: update_pools
   end subroutine active_root_N_uptake

   subroutine myc_scavenger_N_uptake(soilc, vegn, N_uptake_cohorts, myc_efficiency, dt, update_pools)
      import :: soil_BGC_t, vegn_tile_type
      class(soil_BGC_t),  intent(inout) :: soilc
      type(vegn_tile_type), intent(in) :: vegn
      real,intent(out) :: N_uptake_cohorts(:) ! Units: kgN/m2 per individual
      real, intent(in) :: dt  ! dt in years
      logical, intent(in) :: update_pools
      real, intent(out) :: myc_efficiency ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero
   end subroutine

   subroutine myc_miner_N_uptake(soilc,soil,vegn,N_uptake_cohorts,C_uptake_cohorts,total_CO2prod,myc_efficiency,dt,update_pools)
      import :: soil_BGC_t, vegn_tile_type, soil_tile_type
      class(soil_BGC_t),       intent(inout) :: soilc
      type(soil_tile_type), intent(in)    :: soil
      type(vegn_tile_type), intent(in)    :: vegn
      real,    intent(out) :: N_uptake_cohorts(:), C_uptake_cohorts(:)  ! Units kg/m2 of per individual
      real,    intent(out) :: total_CO2prod ! Units of kgC/m2 (not per individual)
      real,    intent(in)  :: dt  ! dt in years
      logical, intent(in)  :: update_pools
      real,    intent(out) :: myc_efficiency  ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero
   end subroutine

   ! transfer fraction of intermediate pools (used for smoothing out contributions
   ! of spiky processes, e.g. harvesting) defined in vegetation tile data structure
   ! to soil BGC pools.
   subroutine spend_intermediate_pools(soilc, vegn)
      import :: soil_BGC_t, vegn_tile_type
      class(soil_BGC_t),    intent(inout) :: soilc
      type(vegn_tile_type), intent(inout) :: vegn
   end subroutine

   subroutine dsdt(soilc, soil, vegn, diag, soilt, theta)
      import soil_BGC_t, soil_tile_type, vegn_tile_type, diag_buff_type
      class(soil_BGC_t), intent(inout)       :: soilc
      type(soil_tile_type), intent(inout) :: soil
      type(vegn_tile_type), intent(inout) :: vegn
      type(diag_buff_type), intent(inout) :: diag
      real                , intent(in)    :: soilt ! average soil temperature, deg K
      real                , intent(in)    :: theta ! average soil moisture
   end subroutine

   subroutine step3(soilc, diag)
      import soil_BGC_t, diag_buff_type
      class(soil_BGC_t),       intent(inout) :: soilc
      type(diag_buff_type), intent(inout) :: diag
   end subroutine

   subroutine redistribute_peat_carbon(soilc)
      import :: soil_BGC_t, vegn_tile_type, soil_tile_type
      class(soil_BGC_t), intent(inout) :: soilc
   end subroutine

end interface

contains

! ============================================================================
!> @brief Move substance from one pool to another
!!
!! This is a utility subroutine that, given an intermediate pool of matte
!! (e.g. C or N), and its spending rate, move the amount of mass corresponding
!! to one fast time step from the pool to the destination. The spending rate
!! is adjusted so that intermediate pool is never depleted below zero.
!!
!! @note
!! Spending rate argument is also updated, to be correctly reported to
!! diagnostics from calling subroutine.
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

! ============================================================================
!!!xz this following subroutine is adopted from CH's code using concentration over water; please note that the units of some input variables are different. I kept tracer_advection_ORI following this subroutine
subroutine tracer_advection(tracer_mass,flow,div,dz,del_tracer,divergence_loss,wl)  ! wl was added here compared to the old version
    real,intent(inout),dimension(:):: tracer_mass  ! Per layer (not per unit water)
    real,intent(in),dimension(:)   :: flow  ! Total flow, not flow rate [mm]
    real,intent(in),dimension(:)   :: div   ! Horizontal divergence (layer total, not rate) [mm]
    real,intent(in),dimension(:)   :: dz    ! Layer thickness
    real,intent(in),dimension(:)   :: wl    ! water content [kg/m^2] by layer before Richards (1:num_l)
    real,intent(out),dimension(:)  :: del_tracer,divergence_loss ! Change in tracer mass, and divergence part

    real,dimension(size(tracer_mass)) :: aaa,bbb,ccc,ddd,wl_litter    ! Matrix coefficients for aaa*dx[i-1] + bbb*dx[i] + ccc*dx[i+1] = ddd
    real,dimension(size(tracer_mass)) :: u_minus,u_plus   ! For weighting of flow upstream/downstream
    real,dimension(size(tracer_mass)) :: tracer_concentration ! [kg C/m^3 soil]
    integer::ll,nlayers
    ! real,dimension(size(tracer_mass)) ::flow_eff ! flow adjusted to be units of [m], weighted by 1/wl  ZACK'S CODE
    real, parameter :: minwl = 0.1 ! [mm] minimum allowed wl
    real, parameter :: dens_h2o=1000.   ! kg/m3
!    real*8,parameter::porosity=0.3  !CH valore inventato  !xz volumn of water over volumn of soil; need to consider to change!!
    !real,intent(in)::theta
    !real*8::dt=1.0/(48.0*365.0)

    nlayers=size(tracer_mass)
    wl_litter(1)=dz(1)     ! m
    do ll=2,nlayers
       wl_litter(ll)=max(wl(ll-1), minwl)/dens_h2o  ! m
    enddo

    tracer_concentration=tracer_mass/wl_litter   ! kg/m3  ! concentration computed over the volume of water

    u_minus = 1.
    where (flow.lt.0.) u_minus = 0.
    do ll = 1, nlayers-1
        u_plus(ll) = 1. - u_minus(ll+1)
    enddo

    ! Top layer, uses upper bound concentration

    ll=1
    aaa(ll)= 0.0 ! flow(ll)*u_minus(ll)
    bbb(ll)= flow(ll)*(1-u_minus(ll)) - flow(ll+1)*(1-u_plus(ll)) - wl_litter(ll)
  !   m           m                          m                           m

    ! divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)   BEN CODE ORIGINAL
    divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)
    ! [kg/m^2]         =       [m]        *    [kg/m^3]

    ccc(ll)= -flow(ll+1)*u_plus(ll)
    ! m

    ddd(ll)= - tracer_concentration(ll)*(bbb(ll)+wl_litter(ll)) - tracer_concentration(ll+1)*ccc(ll)
 !    kg/m2          kg/m3                       m                    kg/m3                 m

    do ll=2,nlayers-1
        !aaa(ll)=flow(ll)*u_minus(ll)     !BEN ORIGINAL
        !bbb(ll)=flow(ll)*(1-u_minus(ll)) - flow(ll+1)*(1-u_plus(ll)) - dz(ll)
        !divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)
        !ccc(ll)=-flow(ll+1)*u_plus(ll)
        !ddd(ll)=-tracer_concentration(ll-1)*aaa(ll) - tracer_concentration(ll)*(bbb(ll)+dz(ll)) - tracer_concentration(ll+1)*ccc(ll)

!Adapted from ZACK's CODE
     aaa(ll)=flow(ll)*u_minus(ll)
        bbb(ll)=flow(ll)*(1-u_minus(ll)) - flow(ll+1)*(1-u_plus(ll)) - wl_litter(ll)
        divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)     ! [kg/m^3]
!         kg/m2            =    m           *    kg/m3
        ccc(ll)=-flow(ll+1)*u_plus(ll)   !m
        ddd(ll)=-tracer_concentration(ll-1)*aaa(ll) - tracer_concentration(ll)*(bbb(ll)+wl_litter(ll)) - tracer_concentration(ll+1)*ccc(ll)
   !    kg/m2
    enddo


    !bottom layer, flow out is zero
   ! ll=nlayers
   ! aaa(ll)=flow(ll)*u_minus(ll)
   ! bbb(ll)= flow(ll)*(1-u_minus(ll)) - dz(ll)
   ! divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)
   ! ccc(ll)= 0.0
   ! ddd(ll)=-tracer_concentration(ll-1)*aaa(ll) - tracer_concentration(ll)*(bbb(ll)+dz(ll))

!Adapted from ZACK's CODE
    ll=nlayers
    aaa(ll)=flow(ll)*u_minus(ll)
    bbb(ll)= flow(ll)*(1-u_minus(ll)) - wl_litter(ll)
    divergence_loss(ll)=max(div(ll),0.0)*tracer_concentration(ll)
    ccc(ll)= 0.0
    ddd(ll)=-tracer_concentration(ll-1)*aaa(ll) - tracer_concentration(ll)*(bbb(ll)+wl_litter(ll))


    !Solve the linear algebra problem
    if(nlayers.gt.1) then
        call tridiag(aaa,bbb,ccc,ddd,del_tracer)  !kg/m3
    else
        del_tracer=0.0
    endif

    del_tracer=del_tracer*wl_litter   !kg/m2   !variazione del tracer
    tracer_mass=tracer_mass+del_tracer    !kg/m2
    divergence_loss=divergence_loss

    where(divergence_loss>tracer_mass) divergence_loss=tracer_mass
end subroutine tracer_advection

end module
