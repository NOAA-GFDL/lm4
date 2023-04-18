module snow_accessors_mod

use land_tile_mod, only: land_tile_type
use snowpack_mod, only: snow_layer_type

implicit none
public

contains

! ============================================================================
! tile existence detector: returns a logical value indicating wether component
! model tile exists or nt
! logical function snow_tile_exists(tile)
!    type(land_tile_type), pointer :: tile
!    snow_tile_exists = associated(tile%snow)
! end function snow_tile_exists

! ============================================================================
! cohort accessor functions: given a pointer to a snowlayer, return a pointer to a
! specific member of the cohort structure

#define DEFINE_SNOWPACK_ACCESSOR_0D(xtype,x) subroutine snowtile_ ## x ## _ptr(t,p);\
type(land_tile_type),pointer::t;xtype,pointer::p;p=>NULL();if(associated(t))then;if(associated(t%snow))p=>t%snow%sp%x;endif;end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR(xtype,x) subroutine snowlayer_ ## x ## _ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x;end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_BC(xtype,x) subroutine snowlayer_ ## x ## _bc_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(1);end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_MD(xtype,x) subroutine snowlayer_ ## x ## _md_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(2);end subroutine

#define DEFINE_SNOWLAYER_ACCESSOR_OM(xtype,x) subroutine snowlayer_ ## x ## _om_ptr(c,p);\
type(snow_layer_type),pointer::c;xtype,pointer::p;p=>NULL();if(associated(c))p=>c%x(3);end subroutine



DEFINE_SNOWPACK_ACCESSOR_0D(integer, nlayers)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topwater)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topwheat)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topsnowdeficit)
DEFINE_SNOWPACK_ACCESSOR_0D(real, topsnowheatdeficit)

DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_tot)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_em)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_bceq_im)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_dendr)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_optd)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_sph)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_rho)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_age)
DEFINE_SNOWPACK_ACCESSOR_0D(real, nearsurf_T)


DEFINE_SNOWLAYER_ACCESSOR(real,T)
DEFINE_SNOWLAYER_ACCESSOR(real,wl)
DEFINE_SNOWLAYER_ACCESSOR(real,ws)
DEFINE_SNOWLAYER_ACCESSOR(real,dz)
DEFINE_SNOWLAYER_ACCESSOR(real,age)
DEFINE_SNOWLAYER_ACCESSOR(real,dendr)
DEFINE_SNOWLAYER_ACCESSOR(real,optd)
DEFINE_SNOWLAYER_ACCESSOR(real,sph)



DEFINE_SNOWLAYER_ACCESSOR_BC(real,wc_im) 
DEFINE_SNOWLAYER_ACCESSOR_MD(real,wc_im) 
DEFINE_SNOWLAYER_ACCESSOR_OM(real,wc_im) 
DEFINE_SNOWLAYER_ACCESSOR_BC(real,wc_em) 
DEFINE_SNOWLAYER_ACCESSOR_MD(real,wc_em) 
DEFINE_SNOWLAYER_ACCESSOR_OM(real,wc_em) 





end module snow_accessors_mod