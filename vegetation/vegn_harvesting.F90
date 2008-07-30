module vegn_harvesting_mod

use fms_mod, only : write_version_number, string, error_mesg, FATAL, NOTE, &
     mpp_pe, write_version_number, file_exist, open_namelist_file, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : axistype, mpp_get_atts, mpp_get_axis_data, &
     mpp_open, mpp_close, MPP_RDONLY, MPP_WRONLY, MPP_ASCII
use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_NTRL, HARV_POOL_SCND, &
     agf_bs, fsc_liv, fsc_wood
use vegn_tile_mod, only : &
     vegn_tile_type
use vegn_cohort_mod, only : &
     vegn_cohort_type, update_biomass_pools

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting

public :: vegn_graze_pasture
public :: vegn_harvest_cropland
public :: vegn_cut_forest
! ==== end of public interface ===============================================

! ==== module constants =====================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_harvesting.F90,v 16.0 2008/07/30 22:30:16 fms Exp $', &
     tagname = '$Name: perth $', &
     module_name = 'vegn_harvesting_mod'

! ==== module data ==========================================================

! ---- namelist variables ---------------------------------------------------
logical :: do_harvesting     = .TRUE.  ! if true, then harvesting of crops and pastures is done
real :: grazing_intensity    = 0.25    ! fraction of biomass removed each time by grazing
real :: grazing_residue      = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: fraction_wood_wasted = 0.25    ! fraction of wood wasted while cutting it
real :: crop_seed_density    = 0.1     ! biomass of seeds left after crop harvesting, kg/m2
namelist/harvesting_nml/ do_harvesting, grazing_intensity, grazing_residue, &
     fraction_wood_wasted, crop_seed_density

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init
  integer :: unit, ierr, io

  call write_version_number(version, tagname)

  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=harvesting_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'harvesting_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  
  if (mpp_pe() == mpp_root_pe()) then
     write (stdlog(), nml=harvesting_nml)
  endif
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (.not.do_harvesting) &
       return ! do nothing if no harvesting requested

  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     call vegn_graze_pasture    (vegn)
  case(LU_CROP)  ! crop
     call vegn_harvest_cropland (vegn)
  end select
end subroutine


! ============================================================================
subroutine vegn_graze_pasture(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars 
  real ::  bdead0, balive0, btotal0; ! initial combined biomass pools
  real ::  bdead1, balive1, btotal1; ! updated combined biomass pools
  type(vegn_cohort_type), pointer :: cc ! shorthand for the current cohort
  integer :: i

  balive0 = 0 ; bdead0 = 0 ;
  balive1 = 0 ; bdead1 = 0 ;

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive0 = balive0 + cc%bl + cc%blv + cc%br
     bdead0  = bdead0  + cc%bwood + cc%bsw
     ! only potential leaves are consumed
     vegn%harv_pool(HARV_POOL_PAST) = vegn%harv_pool(HARV_POOL_PAST) + &
          cc%bliving*cc%Pl*grazing_intensity*(1-grazing_residue) ;
     cc%bliving = cc%bliving - cc%bliving*cc%Pl*grazing_intensity;

     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
 
     ! calculate new combined vegetation biomass pools
     balive1 = balive1 + cc%bl + cc%blv + cc%br
     bdead1  = bdead1  + cc%bwood + cc%bsw
  enddo
  btotal0 = balive0 + bdead0
  btotal1 = balive1 + bdead1

  ! update intermediate soil carbon pools
  vegn%fsc_pool = vegn%fsc_pool + &
       (fsc_liv*(balive0-balive1)+fsc_wood*(bdead0-bdead1))*grazing_residue;
  vegn%ssc_pool = vegn%ssc_pool + &
       ((1-fsc_liv)*(balive0-balive1)+ (1-fsc_wood)*(bdead0-bdead1))*grazing_residue;
end subroutine vegn_graze_pasture


! ================================================================================
subroutine vegn_harvest_cropland(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i
  
  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! use for harvest only aboveg round living biomass and waste the correspondent below living and wood
     vegn%harv_pool(HARV_POOL_CROP) = vegn%harv_pool(HARV_POOL_CROP) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*fraction_harvested;
     vegn%fsc_pool = vegn%fsc_pool + fraction_harvested*(fsc_liv*cc%bliving*cc%Pr + &
          fsc_wood*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));
     vegn%ssc_pool = vegn%ssc_pool + fraction_harvested*((1-fsc_liv)*cc%bliving*cc%Pr + &
          (1-fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));

     cc%bliving = cc%bliving * (1-fraction_harvested);
     cc%bwood   = cc%bwood   * (1-fraction_harvested);
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
  enddo
end subroutine vegn_harvest_cropland


! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: fraction_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: idx; ! index of harvested pool, depends on type of landuse
                  ! of the current patch
  real :: delta
  integer :: i
  

  ! calculate index of the harvesting pool
  select case(vegn%landuse)
  case (LU_NTRL) 
     idx=HARV_POOL_NTRL;
  case (LU_SCND) 
     idx=HARV_POOL_SCND; 
  case default
     call error_mesg('vegn_cut_forest',&
          'landuse type is '//string(vegn%landuse)//': it must be LU_NTRL('//&
          string(LU_NTRL)//') or LU_SCND('//string(LU_SCND)//')', &
          FATAL)     
  end select

  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  fraction_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     vegn%harv_pool(idx) = vegn%harv_pool(idx) + &
          (cc%bwood+cc%bsw) * fraction_harvested*(1-fraction_wood_wasted);
     ! distribute wood and living biomass between fast and slow intermediate 
     ! soil carbon pools according to fractions specified thorough the namelists
     delta = (cc%bwood+cc%bsw)*fraction_harvested*fraction_wood_wasted;
     if(delta<0) call error_mesg('vegn_cut_forest', &
          'harvested amount of dead biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)
     vegn%ssc_pool = vegn%ssc_pool + delta*(1-fsc_wood);
     vegn%fsc_pool = vegn%fsc_pool + delta*   fsc_wood ;

     delta = balive * fraction_harvested;
     if(delta<0) call error_mesg('vegn_cut_forest', &
          'harvested amount of live biomass ('//string(delta)//' kgC/m2) is below zero', &
          FATAL)
     vegn%ssc_pool = vegn%ssc_pool + delta*(1-fsc_liv) ;
     vegn%fsc_pool = vegn%fsc_pool + delta*   fsc_liv  ;

     cc%bliving = cc%bliving*(1-fraction_harvested);
     cc%bwood   = cc%bwood*(1-fraction_harvested);
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(cc);
  enddo
end subroutine vegn_cut_forest

end module 
