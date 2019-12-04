module vegn_accessors_mod

use land_tile_mod, only: land_tile_type
use vegn_cohort_mod, only: vegn_cohort_type
use vegn_data_mod, only: spdata, FORM_WOODY, FORM_GRASS, PT_C3, PT_C4

implicit none
public

contains

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or not
logical function vegn_tile_exists(tile)
   type(land_tile_type), pointer :: tile
   vegn_tile_exists = associated(tile%vegn)
end function vegn_tile_exists

function any_vegn(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = .TRUE.
end function

function is_tree(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_WOODY)
end function

function is_grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS)
end function

function is_c3(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%pt == PT_C3)
end function

function is_c4(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%pt == PT_C4)
end function

function is_c3grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS.and.spdata(cc%species)%pt == PT_C3)
end function

function is_c4grass(cc) result(answer); logical answer
   type(vegn_cohort_type), intent(in) :: cc
   answer = (spdata(cc%species)%lifeform == FORM_GRASS.and.spdata(cc%species)%pt == PT_C4)
end function

! ============================================================================
! cohort accessor functions: given a pointer to cohort, return a pointer to a
! specific member of the cohort structure
#define DEFINE_VEGN_ACCESSOR_0D(xtype,x) subroutine vegn_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x;endif;end subroutine

#define DEFINE_VEGN_ACCESSOR_1D(xtype,x) subroutine vegn_ ## x ## _ptr(t,i,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x(i);endif;end subroutine

#define DEFINE_VEGN_ACCESSOR_2D(xtype,x) subroutine vegn_ ## x ## _ptr(t,i,j,p);\
type(land_tile_type),pointer::t;integer,intent(in)::i,j;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%vegn))p=>t%vegn%x(i,j);endif;end subroutine

#define DEFINE_COHORT_ACCESSOR(xtype,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_COHORT_COMPONENT_ACCESSOR(xtype,component,x) subroutine cohort_ ## x ## _ptr(c,p);\
type(vegn_cohort_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%component%x;end subroutine

DEFINE_VEGN_ACCESSOR_0D(integer,landuse)
DEFINE_VEGN_ACCESSOR_0D(real,age_since_disturbance)
DEFINE_VEGN_ACCESSOR_0D(real,age_since_landuse)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_ag)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_ag)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,fsc_rate_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssc_rate_bg)

DEFINE_VEGN_ACCESSOR_0D(real,fsn_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,fsn_rate_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssn_pool_bg)
DEFINE_VEGN_ACCESSOR_0D(real,ssn_rate_bg)

DEFINE_VEGN_ACCESSOR_2D(real,litter_buff_C)
DEFINE_VEGN_ACCESSOR_2D(real,litter_buff_N)
DEFINE_VEGN_ACCESSOR_2D(real,litter_rate_C)
DEFINE_VEGN_ACCESSOR_2D(real,litter_rate_N)

DEFINE_VEGN_ACCESSOR_0D(real,tc_av)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_phen)
DEFINE_VEGN_ACCESSOR_0D(real,theta_av_fire)
DEFINE_VEGN_ACCESSOR_0D(real,psist_av)
DEFINE_VEGN_ACCESSOR_0D(real,tsoil_av)
DEFINE_VEGN_ACCESSOR_0D(real,precip_av)
DEFINE_VEGN_ACCESSOR_0D(real,fuel)
DEFINE_VEGN_ACCESSOR_0D(real,lambda)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold)
DEFINE_VEGN_ACCESSOR_0D(real,ncm)
DEFINE_VEGN_ACCESSOR_0D(real,t_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,p_ann_acm)
DEFINE_VEGN_ACCESSOR_0D(real,t_cold_acm)
DEFINE_VEGN_ACCESSOR_0D(real,ncm_acm)
DEFINE_VEGN_ACCESSOR_0D(real,treeline_T_accum)
DEFINE_VEGN_ACCESSOR_0D(real,treeline_N_accum)

DEFINE_VEGN_ACCESSOR_0D(real,tc_pheno)
DEFINE_VEGN_ACCESSOR_0D(real,tc_dorm)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_pool)
DEFINE_VEGN_ACCESSOR_0D(real,csmoke_rate)
DEFINE_VEGN_ACCESSOR_0D(real,drop_wl)
DEFINE_VEGN_ACCESSOR_0D(real,drop_ws)
DEFINE_VEGN_ACCESSOR_0D(real,drop_hl)
DEFINE_VEGN_ACCESSOR_0D(real,drop_hs)
DEFINE_VEGN_ACCESSOR_0D(real,nsmoke_pool)

DEFINE_VEGN_ACCESSOR_1D(real,harv_pool_C)
DEFINE_VEGN_ACCESSOR_1D(real,harv_pool_N)
DEFINE_VEGN_ACCESSOR_1D(real,harv_rate_C)
DEFINE_VEGN_ACCESSOR_1D(real,drop_seed_C)
DEFINE_VEGN_ACCESSOR_1D(real,drop_seed_N)

DEFINE_COHORT_ACCESSOR(real,Tv)
DEFINE_COHORT_ACCESSOR(real,Wl)
DEFINE_COHORT_ACCESSOR(real,Ws)

DEFINE_COHORT_ACCESSOR(integer,species)
DEFINE_COHORT_ACCESSOR(real,bl)
DEFINE_COHORT_ACCESSOR(real,br)
DEFINE_COHORT_ACCESSOR(real,blv)
DEFINE_COHORT_ACCESSOR(real,bsw)
DEFINE_COHORT_ACCESSOR(real,bwood)
DEFINE_COHORT_ACCESSOR(real,nsc)
DEFINE_COHORT_ACCESSOR(real,bseed)
DEFINE_COHORT_ACCESSOR(real,bsw_max)
DEFINE_COHORT_ACCESSOR(real,bl_max)
DEFINE_COHORT_ACCESSOR(real,br_max)
DEFINE_COHORT_ACCESSOR(real,dbh)
DEFINE_COHORT_ACCESSOR(real,crownarea)
DEFINE_COHORT_ACCESSOR(real,bliving)
DEFINE_COHORT_ACCESSOR(real,nindivs)
DEFINE_COHORT_ACCESSOR(integer,layer)
DEFINE_COHORT_ACCESSOR(integer,firstlayer)
DEFINE_COHORT_ACCESSOR(integer,status)
DEFINE_COHORT_ACCESSOR(real,leaf_age)
DEFINE_COHORT_ACCESSOR(real,age)
DEFINE_COHORT_ACCESSOR(real,npp_previous_day)
DEFINE_COHORT_ACCESSOR(real,growth_previous_day)
DEFINE_COHORT_ACCESSOR(real,growth_previous_day_tmp)
DEFINE_COHORT_ACCESSOR(real,BM_ys)
DEFINE_COHORT_ACCESSOR(real,DBH_ys)
DEFINE_COHORT_ACCESSOR(real,topyear)
DEFINE_COHORT_ACCESSOR(real,gdd)
DEFINE_COHORT_ACCESSOR(real,height)
DEFINE_COHORT_ACCESSOR(real,theta_av_phen)
DEFINE_COHORT_ACCESSOR(real,psist_av_phen)

DEFINE_COHORT_ACCESSOR(real,scav_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,scav_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,mine_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,mine_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,nfix_C_reservoir)
DEFINE_COHORT_ACCESSOR(real,nfix_N_reservoir)
DEFINE_COHORT_ACCESSOR(real,scav_C)
DEFINE_COHORT_ACCESSOR(real,scav_N)
DEFINE_COHORT_ACCESSOR(real,mine_C)
DEFINE_COHORT_ACCESSOR(real,mine_N)
DEFINE_COHORT_ACCESSOR(real,nfix_C)
DEFINE_COHORT_ACCESSOR(real,nfix_N)
DEFINE_COHORT_ACCESSOR(real,stored_N)
DEFINE_COHORT_ACCESSOR(real,wood_N)
DEFINE_COHORT_ACCESSOR(real,sapwood_N)
DEFINE_COHORT_ACCESSOR(real,leaf_N)
DEFINE_COHORT_ACCESSOR(real,seed_N)
DEFINE_COHORT_ACCESSOR(real,root_N)
DEFINE_COHORT_ACCESSOR(real,total_N)
DEFINE_COHORT_ACCESSOR(real,nitrogen_stress)
DEFINE_COHORT_ACCESSOR(real,scav_mgain_smoothed)
DEFINE_COHORT_ACCESSOR(real,mine_mgain_smoothed)
DEFINE_COHORT_ACCESSOR(real,nfix_mgain_smoothed)
DEFINE_COHORT_ACCESSOR(real,exud_mgain_smoothed)

DEFINE_COHORT_ACCESSOR(real,max_monthly_scav_alloc)
DEFINE_COHORT_ACCESSOR(real,max_monthly_mine_alloc)
DEFINE_COHORT_ACCESSOR(real,max_monthly_Nfix_alloc)
DEFINE_COHORT_ACCESSOR(real,Nfix_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,mine_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,scav_alloc_accum)
DEFINE_COHORT_ACCESSOR(real,max_Nfix_alloc)
DEFINE_COHORT_ACCESSOR(real,max_mine_alloc)
DEFINE_COHORT_ACCESSOR(real,max_scav_alloc)
DEFINE_COHORT_ACCESSOR(real,nitrogen_stress_smoothed)

! wolf
DEFINE_COHORT_ACCESSOR(real,psi_r)
DEFINE_COHORT_ACCESSOR(real,psi_x)
DEFINE_COHORT_ACCESSOR(real,psi_l)
DEFINE_COHORT_ACCESSOR(real,Kxa)
DEFINE_COHORT_ACCESSOR(real,Kla)

! ens for branch sapwood
DEFINE_COHORT_ACCESSOR(real,brsw)

end module vegn_accessors_mod