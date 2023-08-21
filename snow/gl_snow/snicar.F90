
module snicar_mod

#include <fms_platform.h>
#include "../../shared/debug.inc"
#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
use fms_mod, only : file_exist, check_nml_error, &
   close_file, stdlog, read_data, error_mesg, FATAL, WARNING, NOTE, field_size, write_data, mpp_pe, mpp_root_pe
use land_debug_mod, only:  is_watch_point, is_watch_cell,check_var_range, set_current_point, land_error_message 
use land_data_mod, only : lnd, log_version
use snow_constants_mod
use constants_mod, only : PI
use snowpack_mod
! use snow_evolution_mod, only: compute_snow_grain_shape

implicit none
private

public :: SNICAR_RT           
public :: SNICAR_AD_RT        
public :: read_snicar_optics_data 
public :: compute_snicar_albedo 
public :: read_snow_snicar_namelist

!! ENRICO EZSNOW NOTES
!! //TODO check optical properties loaded - old file is ok for new snicar version?
!! //TODO check units for input data (radius, conc. of imp.)
!! //TODO cleanup old snicar versions
!! //TODO check sno_fs and sno_AR values, see that shape is read correctly layer-by-layer
!! //TODO update flx_wgt_dif - flx_wgt_dir?  in RT_HE
!! //TODO is_dust_internal_mixing, is_BC_internal_mixing <=> snicar_snobc_intmix, snicar_snodst_intmix 
!! //TODO: ask to make sure order of layering vs variables loaded here passed from lm4p2 (BOTH START FROM TOP?)
!!  IN GLASS ordering of layers is from the TOP
!! //TODO: IN GLASS-SNICAR, add options for clean snow or only-dust snow
!!
!! He version obtained from
!! https://github.com/cenlinhe/CTSM/blob/snicar_allupdate/src/biogeophys/SnowSnicarMod.F90
!!

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'snicar_mod' 
#include "../../shared/version_variable.inc"


! ---- namelist
! integer, parameter :: snow_shape_defined = 1 ! IF SET TO ZERO, USE ACTUAL VALUE. ELSE, 1-2-3-4 IMPOSE SHAPE
! logical, parameter :: use_snicar_ad = .TRUE.
! logical, parameter :: is_dust_internal_mixing = .FALSE. ! FOR NOW FORCE IT
! logical, parameter :: is_BC_internal_mixing = .FALSE.   ! FOR NOW FORCE IT
! integer, parameter :: snicar_atm_type = 0 ! default

integer :: snow_shape_defined = 0 ! IF SET TO ZERO, USE ACTUAL VALUE. ELSE, 1-2-3-4 IMPOSE SHAPE
logical :: use_snicar_ad = .TRUE.
logical :: is_dust_internal_mixing = .TRUE. ! FOR NOW FORCE IT
logical :: is_BC_internal_mixing = .TRUE.   ! FOR NOW FORCE IT
integer :: snicar_atm_type = 0 ! default (midlatitude winter)
CHARACTER(LEN=22) :: ncid = "snicar_optics.nc"

namelist /snow_snicar_nml/ &
   snow_shape_defined , use_snicar_ad, is_BC_internal_mixing, is_dust_internal_mixing, snicar_atm_type, ncid
! ---- end of namelist


! // TODO: remove
integer, parameter :: num_nourbanc = 1 ! EZDEV

  ! EZSNOW - added for new HE version of snicar:: (when we are sure they are the same, update)
  logical,  public :: snicar_snobc_intmix =    .TRUE.   ! internal mixing of BC? ! EZSNOW ASK
  logical,  public :: snicar_snodst_intmix =   .TRUE.    ! internal mixing of DUST? ! EZSNOW ASK
  !
  integer,  public, parameter :: sno_nbr_aer =   8        ! number of aerosol species in snowpack
  logical,  public, parameter :: DO_SNO_OC =    .false.   ! parameter to include organic carbon (OC)
  logical,  public, parameter :: DO_SNO_AER =   .true.    ! parameter to include aerosols 
  integer,  parameter :: numrad_snw  =   5               ! number of spectral bands used in snow model [nbr]
  integer,  parameter :: nir_bnd_bgn =   2               ! first band index in near-IR spectrum [idx]
  integer,  parameter :: nir_bnd_end =   5               ! ending near-IR band index [idx]
  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]

  real, public, parameter :: snw_rds_min = 54.526 ! minimum allowed snow effective radius (also "fresh snow" value) [microns]
  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_int = nint(snw_rds_min) ! minimum allowed snow effective radius as integer [microns]
!   real, parameter :: snw_rds_max     = 1500.0      ! maximum allowed snow effective radius [microns]
  real, parameter :: min_snw = 1.0E-30            ! minimum snow mass required for SNICAR RT calculation [kg m-2]
!   real, parameter :: tim_cns_bc_rmv  = 2.2E-8     ! time constant for removal of BC in snow on sea-ice [s-1] (50% mass removal/year)
!   real, parameter :: tim_cns_oc_rmv  = 2.2E-8     ! time constant for removal of OC in snow on sea-ice [s-1] (50% mass removal/year)
!   real, parameter :: tim_cns_dst_rmv = 2.2E-8     ! time constant for removal of dust in snow on sea-ice [s-1] (50% mass removal/year)

  ! direct-beam weighted ice optical properties
  real :: ss_alb_snw_drc     (idx_Mie_snw_mx,numrad_snw);
  real :: asm_prm_snw_drc    (idx_Mie_snw_mx,numrad_snw);
  real :: ext_cff_mss_snw_drc(idx_Mie_snw_mx,numrad_snw);

  ! diffuse radiation weighted ice optical properties
  real :: ss_alb_snw_dfs     (idx_Mie_snw_mx,numrad_snw);
  real :: asm_prm_snw_dfs    (idx_Mie_snw_mx,numrad_snw);
  real :: ext_cff_mss_snw_dfs(idx_Mie_snw_mx,numrad_snw);

  !!! direct & diffuse
  real :: flx_wgt_dif    (6, numrad_snw); ! 6 atmospheric types
  real :: flx_wgt_dir    (6,90,numrad_snw); ! 6 atmospheric types, 0-89 SZA

  ! hydrophiliic BC
  real :: ss_alb_bc1     (numrad_snw);
  real :: asm_prm_bc1    (numrad_snw);
  real :: ext_cff_mss_bc1(numrad_snw);
  ! hydrophobic BC
  real :: ss_alb_bc2     (numrad_snw);
  real :: asm_prm_bc2    (numrad_snw);
  real :: ext_cff_mss_bc2(numrad_snw);
  ! hydrophobic OC
  real :: ss_alb_oc1     (numrad_snw);
  real :: asm_prm_oc1    (numrad_snw);
  real :: ext_cff_mss_oc1(numrad_snw);
  ! hydrophilic OC
  real :: ss_alb_oc2     (numrad_snw);
  real :: asm_prm_oc2    (numrad_snw);
  real :: ext_cff_mss_oc2(numrad_snw);
  ! dust species 1:
  real :: ss_alb_dst1     (numrad_snw);
  real :: asm_prm_dst1    (numrad_snw);
  real :: ext_cff_mss_dst1(numrad_snw);
  ! dust species 2:
  real :: ss_alb_dst2     (numrad_snw);
  real :: asm_prm_dst2    (numrad_snw);
  real :: ext_cff_mss_dst2(numrad_snw);
  ! dust species 3:
  real :: ss_alb_dst3     (numrad_snw);
  real :: asm_prm_dst3    (numrad_snw);
  real :: ext_cff_mss_dst3(numrad_snw);
  ! dust species 4:
  real :: ss_alb_dst4     (numrad_snw);
  real :: asm_prm_dst4    (numrad_snw);
  real :: ext_cff_mss_dst4(numrad_snw);


contains


subroutine read_snow_snicar_namelist()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines
  integer :: l            ! layer iterator

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=snow_snicar_nml, iostat=io)
  ierr = check_nml_error(io, 'snow_snicar_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=snow_snicar_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'snow_snicar_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=snow_snicar_nml)
  endif

  ! read optical properties now
  call read_snicar_optics_data()

end subroutine read_snow_snicar_namelist


   subroutine compute_snicar_albedo(s, cosz, subs_adif)
   class(snowpack_t), intent(inout) :: s !< state of snowpack
   real, intent(in) :: cosz ! cosine of solar zenith angle
   real, dimension(NBANDS), intent(in) :: subs_adif ! [diffuse] albedo of substrate

   ! LOCAL INPUT TO SNICAR
   integer flg_slr ! 1= direct light; 2=diffuse light
   integer flg_snw_ice ! 1=CLM, 2=CSIM ! //TODO remove
   real :: wsmat(1, s%nlayers)
   real :: wlmat(1, s%nlayers)
   integer :: remat(1, s%nlayers)
   integer :: shmat(1, s%nlayers) ! grain shape
   ! logical :: iemat(1, s%nlayers) ! internal vs external mixing
   real :: trmat(1, s%nlayers, 8) ! default snicar number of tracers => sno_nbr_aer = 8
   real :: subs_adif_mat(1, NBANDS)
   ! LOCAL OUTPUT FROM SNICAR
   real albsnd(1, NBANDS) ! direct albedo of snowpack
   real albsni(1, NBANDS) ! diffuse albedo of snowpack
   real flx_absd_snw(1, s%nlayers, NBANDS) ! fraction of SWdown absorbed in each layer
   real flx_absi_snw(1, s%nlayers, NBANDS) ! fraction of SWdown absorbed in each layer

   integer il, ib,it, NL
   ! real sum_absd_snw(NBANDS), sum_absi_snw(NBANDS)

   NL = s%nlayers
   flg_snw_ice = 1 ! set to CLM - remove
   do ib=1, NBANDS
      subs_adif_mat(1,ib)=subs_adif(ib)
   enddo
   do il=1, NL ! SNICAR STARTS FROM THE BOTTOM
      ! ill = il
      ! ill = NL-il+1 ! TO INVERT LAYER ORDERING - NOT NEEDED HERE
      wsmat(1,il) = s%snow(il)%ws
      wlmat(1,il) = s%snow(il)%wl
      ! remat(1,il) = int(s%snow(il)%optd/2.0*1E6) ! from optd [m] to GRAIN RADIUS IN [/mu m]
      remat(1,il) = min(snw_rds_max_tbl, max(snw_rds_min_tbl, int(s%snow(il)%optd/2.0*1E6))) ! from optd [m] to GRAIN RADIUS IN [/mu m]
      ! // TODO: for now sum internally mixed and externally mixed tracers
      ! // wc[mg]-> compute concentration in kg/kg

      ! init conc values - let's start with clean snow
         trmat(1, il, 1) = 0.0
         trmat(1, il, 2) = 0.0
         trmat(1, il, 3) = 0.0
         trmat(1, il, 4) = 0.0
         trmat(1, il, 5) = 0.0
         trmat(1, il, 6) = 0.0
         trmat(1, il, 7) = 0.0
         trmat(1, il, 8) = 0.0

      ! BLACK CARBON (1=PHI, 2=PHO)
      if (lap_albedo_include_bc) then
      trmat(1, il, 1) = 0.5 * 1E-6 * (s%snow(il)%wc_em(1) + s%snow(il)%wc_im(1) )/(s%snow(il)%ws+s%snow(il)%wl)
      trmat(1, il, 2) = 0.5 * 1E-6 * (s%snow(il)%wc_em(1) + s%snow(il)%wc_im(1) )/(s%snow(il)%ws+s%snow(il)%wl)
      endif
      ! ! ORGANIC CARBON (1=PHI, 2=PHO)
      if (lap_albedo_include_om) then
      trmat(1, il, 3) = 0.5 * 1E-6 * (s%snow(il)%wc_em(3) + s%snow(il)%wc_im(3) )/(s%snow(il)%ws+s%snow(il)%wl)
      trmat(1, il, 4) = 0.5 * 1E-6 * (s%snow(il)%wc_em(3) + s%snow(il)%wc_im(3) )/(s%snow(il)%ws+s%snow(il)%wl)
      endif
      ! ! MINERAL DUST for various size bins
      if (lap_albedo_include_md) then
      trmat(1, il, 5) = 0.25 * 1E-6 * (s%snow(il)%wc_em(2) + s%snow(il)%wc_im(2) )/(s%snow(il)%ws+s%snow(il)%wl)
      trmat(1, il, 6) = 0.25 * 1E-6 * (s%snow(il)%wc_em(2) + s%snow(il)%wc_im(2) )/(s%snow(il)%ws+s%snow(il)%wl)
      trmat(1, il, 7) = 0.25 * 1E-6 * (s%snow(il)%wc_em(2) + s%snow(il)%wc_im(2) )/(s%snow(il)%ws+s%snow(il)%wl)
      trmat(1, il, 8) = 0.25 * 1E-6 * (s%snow(il)%wc_em(2) + s%snow(il)%wc_im(2) )/(s%snow(il)%ws+s%snow(il)%wl)
      endif





      ! GET INTERNAL VS EXTERNAL MIXING RATIO FOR EACH TRACER
      ! SKIP FOR NOW, //TODO
      ! USE SET VALUE IN NAMELIST FOR NOW


      ! GET SNOW SHAPE FOR EACH LAYER
      if (snow_shape_defined==0) then
         ! IN THIS CASE USE ACTUAL COMPUTED GRAIN SHAPE
         call compute_snow_grain_shape(s%snow(il)%dendr, s%snow(il)%sph, shmat(1, il) )
      else ! case equal 1-2-3-4: FORCE GRAIN SHAPE
         shmat(1, il) = snow_shape_defined
      endif


   enddo

      if(is_watch_point()) then
         write(*,*) "##### compute_snicar_albedo - start checkpoint 1 #####"
         write(*,*) "snow_shape_defined = ", snow_shape_defined
         write(*,*) "use_snicar_ad = ", use_snicar_ad
         write(*,*) "is_dust_internal_mixing = ", is_dust_internal_mixing
         write(*,*) "is_BC_internal_mixing = ", is_BC_internal_mixing
         write(*,*) "snicar_atm_type = ", snicar_atm_type
         write(*,*) "cosz = ", cosz
         write(*,*) "subs_adif = ", subs_adif
         write(*,*) "trmat = ", trmat
         write(*,*) "wlmat = ", wlmat
         write(*,*) "wsmat = ", wsmat
         write(*,*) "shmat = ", shmat
         write(*,*) "remat = ", remat
         write(*,*) "print the state of snowpack:"
         call s%print()
         write(*,*) "##### compute_snicar_albedo - end checkpoint 1 #####"
     endif


     if(is_watch_point()) then
      write(*,*) "##### compute_snicar_albedo - start checkpoint 1B #####"


      write(*,*) "size of Mie parameters arrays:"
      write(*,*) 'size of ss_alb_ice_drc', shape(ss_alb_snw_drc)
      write(*,*) 'size of asm_prm_ice_drc', shape(asm_prm_snw_drc)
      write(*,*) 'size of ext_cff_mss_ice_drc', shape(ext_cff_mss_snw_drc)
      write(*,*) "direct-beam snow Mie parameters - min radius:"
      write(*,*) 'ss_alb_ice_drc', ss_alb_snw_drc(1,:)
      write(*,*) 'asm_prm_ice_drc', asm_prm_snw_drc(1,:)
      write(*,*) 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc(1,:)

      write(*,*) "diffuse snow Mie parameters - min radius:"
      write(*,*)'ss_alb_ice_dfs', ss_alb_snw_dfs(1,:)           
      write(*,*)'asm_prm_ice_dfs', asm_prm_snw_dfs(1,:)         
      write(*,*)'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs(1,:) 

      write(*,*) "direct-beam snow Mie parameters - max radius:"
      write(*,*) 'ss_alb_ice_drc', ss_alb_snw_drc(idx_Mie_snw_mx,:)
      write(*,*) 'asm_prm_ice_drc', asm_prm_snw_drc(idx_Mie_snw_mx,:)
      write(*,*) 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc(idx_Mie_snw_mx,:)

      write(*,*) "diffuse snow Mie parameters - max radius:"
      write(*,*)'ss_alb_ice_dfs', ss_alb_snw_dfs(idx_Mie_snw_mx,:)           
      write(*,*)'asm_prm_ice_dfs', asm_prm_snw_dfs(idx_Mie_snw_mx,:)         
      write(*,*)'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs(idx_Mie_snw_mx,:) 

   if (snicar_atm_type > 0)then
      write(*,*) "Solar spectrum weights:"
      write(*,*) "flx_wgt_dir", flx_wgt_dir
      write(*,*) "flx_wgt_dif", flx_wgt_dif
   endif

   write(*,*) "hydrophiliic BC:"
  write(*,*) ss_alb_bc1     
  write(*,*) asm_prm_bc1    
  write(*,*) ext_cff_mss_bc1
  write(*,*) "hydrophobic BC:"
  write(*,*) ss_alb_bc2     
  write(*,*) asm_prm_bc2    
  write(*,*) ext_cff_mss_bc2
  write(*,*) "hydrophobic OC:"
  write(*,*) ss_alb_oc1     
  write(*,*) asm_prm_oc1    
  write(*,*) ext_cff_mss_oc1
  write(*,*) "hydrophilic OC:"
  write(*,*) ss_alb_oc2     
  write(*,*) asm_prm_oc2    
  write(*,*) ext_cff_mss_oc2
  write(*,*) "dust species 1::"
  write(*,*) ss_alb_dst1     
  write(*,*) asm_prm_dst1    
  write(*,*) ext_cff_mss_dst1
  write(*,*) "dust species 2::"
  write(*,*) ss_alb_dst2     
  write(*,*) asm_prm_dst2    
  write(*,*) ext_cff_mss_dst2
  write(*,*) "dust species 3::"
  write(*,*) ss_alb_dst3     
  write(*,*) asm_prm_dst3    
  write(*,*) ext_cff_mss_dst3
  write(*,*) "dust species 4::"
  write(*,*) ss_alb_dst4     
  write(*,*) asm_prm_dst4    
  write(*,*) ext_cff_mss_dst4
      write(*,*) "##### compute_snicar_albedo - end checkpoint 1B #####"
  endif



      ! if ((s%depth()>1.0).and.((s%nlayers > 5).and.(cosz > 0.3))) then
      !    write(*,*) "snow depth = ", s%depth()
      !    write(*,*) "snow nlayers = ", s%nlayers
      !    call s%print()
      !    call land_error_message("EZTEMP: here snow larger than threshold!", severity=FATAL)
      ! endif


      !  flg_slr = 1; ! direct light
      !  if (use_snicar_ad) then
      !      call SNICAR_AD_RT(s%nlayers, flg_snw_ice, (/ cosz /), flg_slr, &   ! INPUT
      !                        wlmat, wsmat, remat, shmat, trmat, subs_adif_mat, &     ! INPUT
      !                        albsnd, flx_absd_snw)                            ! OUTPUT
      !  else
      !      call SNICAR_RT(s%nlayers, flg_snw_ice, (/ cosz /), flg_slr, &      ! INPUT
      !                        wlmat, wsmat, remat, trmat, subs_adif_mat, &     ! INPUT
      !                        albsnd, flx_absd_snw)                            ! OUTPUT
      !  endif 
      !  flg_slr = 2; ! diffuse light
      !  if (use_snicar_ad) then
      !      call SNICAR_AD_RT(s%nlayers, flg_snw_ice, (/ cosz /), flg_slr, &   ! INPUT
      !                        wlmat, wsmat, remat, shmat, trmat, subs_adif_mat, &     ! INPUT
      !                        albsni, flx_absi_snw)                            ! OUTPUT
      !  else
      !      call SNICAR_RT(s%nlayers, flg_snw_ice, (/ cosz /), flg_slr, &      ! INPUT
      !                        wlmat, wsmat, remat, trmat, subs_adif_mat, &     ! INPUT
      !                        albsni, flx_absi_snw)                            ! OUTPUT
      !  endif 


       flg_slr = 1; ! direct light
       call SNICAR_RT_HE(s%nlayers, flg_snw_ice, &
         (/ cosz /), flg_slr, wlmat, wsmat, remat, shmat,   &
         trmat, subs_adif_mat, albsnd, flx_absd_snw)
       flg_slr = 2; ! diffuse light
       call SNICAR_RT_HE(s%nlayers, flg_snw_ice, &
         (/ cosz /), flg_slr, wlmat, wsmat, remat, shmat,   &
         trmat, subs_adif_mat, albsni, flx_absi_snw)


   ! write(*,*) "COMPUTING SNICAR SNOW ALBEDO!"
   ! write(*,*) "snow nlayers, depth = ", s%nlayers, s%depth()
   ! write(*,*) "snicar direct albedo = ", albsnd
   ! write(*,*) "snicar diffuse albedo = ", albsni
   !    write(*,*) "snicar  direct absorbed flux fractions = ", flx_absd_snw(1, :, 1)
   ! write(*,*) "snicar  diffuse absorbed flux fractions = ", flx_absi_snw(1, :, 1)
   ! write(*,*) "snicar VIS direct absorbed flux fractions = ", albsnd(1,1) + sum(flx_absd_snw(1, :, 1))
   ! write(*,*) "snicar VIS diffuse absorbed flux fractions = ", albsni(1,1) + sum(flx_absi_snw(1, :, 1))
   !    write(*,*) "snicar NIR direct absorbed flux fractions = ",albsnd(1,2) + sum(flx_absd_snw(1, :, 2))
   ! write(*,*) "snicar NIR diffuse absorbed flux fractions = ",albsni(1,2) + sum(flx_absi_snw(1, :, 2))

         ! SAVE SNICAR RESULTS IN SNOWPACK OBJECT
         s%snow_refl_dir = albsnd(1,:) 
         s%snow_refl_dif = albsni(1,:) 
         if (allocated(s%sw_frac_dir)) then
            ! write(*,*) "original size frac_dir = ", size(s%sw_frac_dir), NBANDS*s%nlayers
            if (size(s%sw_frac_dir)<NBANDS*s%nlayers) then
               deallocate(s%sw_frac_dir)
            endif
         endif
         if (allocated(s%sw_frac_dif)) then
            ! write(*,*) "original size frac_dir = ", size(s%sw_frac_dir), NBANDS*s%nlayers
            if (size(s%sw_frac_dif)<NBANDS*s%nlayers) then
               deallocate(s%sw_frac_dif)
            endif
         endif
         if (.not. allocated(s%sw_frac_dir)) then
            allocate(s%sw_frac_dir(s%nlayers, NBANDS))
         endif
         if (.not. allocated(s%sw_frac_dif)) then
            allocate(s%sw_frac_dif(s%nlayers, NBANDS))
         endif
         ! write(*,*) "init sums"
         ! sum_absd_snw = 0.0
         ! sum_absi_snw = 0.0
         ! do il=1, NL
         !    do ib=1,NBANDS
         !       sum_absd_snw(ib) = sum_absd_snw(ib) + flx_absd_snw(1, il, ib) 
         !       sum_absi_snw(ib) = sum_absi_snw(ib) + flx_absi_snw(1, il, ib) 
         !    enddo
         ! enddo
         do il=1, NL
            do ib=1,NBANDS
               ! if (sum_absd_snw(ib) >(1E-7)) then
               !    s%sw_frac_dir(il,ib) = flx_absd_snw(1, il, ib) / sum_absd_snw(ib)
               ! else
               !    s%sw_frac_dir(il,ib) = flx_absd_snw(1, il, ib) 
               ! endif
               ! if (sum_absi_snw(ib) >(1E-7)) then
               !    s%sw_frac_dif(il,ib) = flx_absi_snw(1, il, ib) / sum_absi_snw(ib)
               ! else
               !    s%sw_frac_dif(il,ib) = flx_absi_snw(1, il, ib) 
               ! endif
               s%sw_frac_dir(il,ib) = flx_absd_snw(1, il, ib) 
               s%sw_frac_dif(il,ib) = flx_absi_snw(1, il, ib) 
            enddo
         enddo
         ! write(*,*) "end SNICAR albedo subroutine"

         if(is_watch_point()) then
            write(*,*) "##### compute_snicar_albedo - start checkpoint 2 #####"
            write(*,*) "Energy absorbed in the snowpack [for unit source]:"
            write(*,*) "sw_frac_dir VIS = ", flx_absd_snw(1, :, 1) 
            write(*,*) "sw_frac_dir NIR = ", flx_absd_snw(1, :, 2) 
            write(*,*) "sw_frac_dif VIS = ", flx_absi_snw(1, :, 1) 
            write(*,*) "sw_frac_dif NIR = ", flx_absi_snw(1, :, 2) 
            write(*,*) "Direct albedo  VIS - NIR",          albsnd(1,:) 
            write(*,*) "Diffuse albedo VIS - NIR",          albsni(1,:) 
            write(*,*) "##### compute_snicar_albedo - end checkpoint 2 #####"
        endif
   end subroutine compute_snicar_albedo


  subroutine SNICAR_RT (nlevsno, flg_snw_ice, &
                        coszen, flg_slr_in, h2osno_liq, h2osno_ice, snw_rds,   &
                        mss_cnc_aer_in, albsfc, albout, flx_abs)
    !
    ! !DESCRIPTION:
    ! Determine reflectance of, and vertically-resolved solar absorption in,
    ! snow with impurities.
    !
    ! Original references on physical models of snow reflectance include:
    ! Wiscombe and Warren [1980] and Warren and Wiscombe [1980],
    ! Journal of Atmospheric Sciences, 37,
    !
    ! The multi-layer solution for multiple-scattering used here is from:
    ! Toon et al. [1989], Rapid calculation of radiative heating rates
    ! and photodissociation rates in inhomogeneous multiple scattering atmospheres,
    ! J. Geophys. Res., 94, D13, 16287-16301
    !
    ! The implementation of the SNICAR model in CLM/CSIM is described in:
    ! Flanner, M., C. Zender, J. Randerson, and P. Rasch [2007],
    ! Present-day climate forcing and response from black carbon in snow,
    ! J. Geophys. Res., 112, D11202, doi: 10.1029/2006JD008003
    !

    integer, INTENT(IN) :: nlevsno
    integer           , intent(in)  :: flg_snw_ice                                        ! flag: =1 when called from CLM, =2 when called from CSIM
   !  type (bounds_type), intent(in)  :: bounds
   !  integer           , intent(in)  :: num_nourbanc                                       ! number of columns in non-urban filter
   !  integer           , intent(in)  :: filter_nourbanc(:)                                 ! column filter for non-urban points
    real          , intent(in)  :: coszen         ( 1: )                    ! cosine of solar zenith angle for next time step (col) [unitless]
    integer           , intent(in)  :: flg_slr_in                                         ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux
    real          , intent(in)  :: h2osno_liq     ( 1: , -nlevsno+1: )      ! liquid water content (col,lyr) [kg/m2]
    real          , intent(in)  :: h2osno_ice     ( 1: , -nlevsno+1: )      ! ice content (col,lyr) [kg/m2]
    integer           , intent(in)  :: snw_rds        ( 1: , -nlevsno+1: )      ! snow effective radius (col,lyr) [microns, m^-6]
    real          , intent(in)  :: mss_cnc_aer_in ( 1: , -nlevsno+1: , 1: ) ! mass concentration of all aerosol species (col,lyr,aer) [kg/kg]
    real          , intent(in)  :: albsfc         ( 1: , 1: )               ! albedo of surface underlying snow (col,bnd) [frc]
    real          , intent(out) :: albout         ( 1: , 1: )               ! snow albedo, averaged into 2 bands (=0 if no sun or no snow) (col,bnd) [frc]
    real          , intent(out) :: flx_abs        ( 1: , -nlevsno+1: , 1: ) ! absorbed flux in each layer per unit flux incident (col, lyr, bnd)
    !
    ! !LOCAL VARIABLES:
    !
    ! variables for snow radiative transfer calculations

    ! Local variables representing single-column values of arrays:
    integer :: snl_lcl                            ! negative number of snow layers [nbr]
    integer :: snw_rds_lcl(-nlevsno+1:0)          ! snow effective radius [m^-6]
    real:: flx_slrd_lcl(1:numrad_snw)         ! direct beam incident irradiance [W/m2] (set to 1)
    real:: flx_slri_lcl(1:numrad_snw)         ! diffuse incident irradiance [W/m2] (set to 1)
    real:: mss_cnc_aer_lcl(-nlevsno+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
    real:: h2osno_lcl                         ! total column snow mass [kg/m2]
    real:: h2osno_liq_lcl(-nlevsno+1:0)       ! liquid water mass [kg/m2]
    real:: h2osno_ice_lcl(-nlevsno+1:0)       ! ice mass [kg/m2]
    real:: albsfc_lcl(1:numrad_snw)           ! albedo of underlying surface [frc]
    real:: ss_alb_snw_lcl(-nlevsno+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
    real:: asm_prm_snw_lcl(-nlevsno+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
    real:: ext_cff_mss_snw_lcl(-nlevsno+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
    real:: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc]
    real:: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
    real:: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]

    ! Other local variables
    integer :: APRX_TYP                           ! two-stream approximation type
                                                  ! (1=Eddington, 2=Quadrature, 3=Hemispheric Mean) [nbr]
    integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
                                                  ! (1= use, 0= don't use)
    real:: flx_wgt(1:numrad_snw)              ! weights applied to spectral bands,
                                                  ! specific to direct and diffuse cases (bnd) [frc]

    integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
                                                  ! =0 if at least 1 snow layer [flg]
    integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
    integer :: flg_dover                          ! defines conditions for RT redo (explained below)

    real:: albedo                             ! temporary snow albedo [frc]
    real:: flx_sum                            ! temporary summation variable for NIR weighting
    real:: albout_lcl(numrad_snw)             ! snow albedo by band [frc]
    real:: flx_abs_lcl(-nlevsno+1:1,numrad_snw)! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]

    real:: L_snw(-nlevsno+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
    real:: tau_snw(-nlevsno+1:0)              ! snow optical depth (lyr) [unitless]
    real:: L_aer(-nlevsno+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2]
    real:: tau_aer(-nlevsno+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
    real:: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
    real:: tau_elm(-nlevsno+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless]
    real:: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
    real:: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]

    real:: tau(-nlevsno+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
    real:: omega(-nlevsno+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
    real:: g(-nlevsno+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
    real:: tau_star(-nlevsno+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                                                  ! (lyr) [unitless]
    real:: omega_star(-nlevsno+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
    real:: g_star(-nlevsno+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                                                  ! (lyr) [frc]

    integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
    integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= numrad_snw) [idx]
    integer :: rds_idx                            ! snow effective radius index for retrieving
                                                  ! Mie parameters from lookup table [idx]
    integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
    integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
    integer :: fc                                 ! column filter index
    integer :: i                                  ! layer index [idx]
    integer :: j                                  ! aerosol number index [idx]
    integer :: n                                  ! tridiagonal matrix index [idx]
    integer :: m                                  ! secondary layer index [idx]
    integer :: nint_snw_rds_min                   ! nearest integer value of snw_rds_min
    
    real:: F_direct(-nlevsno+1:0)             ! direct-beam radiation at bottom of layer interface (lyr) [W/m^2]
    real:: F_net(-nlevsno+1:0)                ! net radiative flux at bottom of layer interface (lyr) [W/m^2]
    real:: F_abs(-nlevsno+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
    real:: F_abs_sum                          ! total absorbed energy in column [W/m^2]
    real:: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
    real:: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]
    real:: F_sfc_net                          ! net flux at top of snowpack [W/m^2]
    real:: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
    real:: F_direct_btm                       ! direct-beam radiation at bottom of snowpack [W/m^2]
    real:: mu_not                             ! cosine of solar zenith angle (used locally) [frc]

    integer :: err_idx                            ! counter for number of times through error loop [nbr]
    real:: lat_coord                          ! gridcell latitude (debugging only)
    real:: lon_coord                          ! gridcell longitude (debugging only)
    integer :: sfctype                            ! underlying surface type (debugging only)
   !  real:: pi                                 ! 3.1415...

    integer :: nstep

    ! intermediate variables for radiative transfer approximation:
    real:: gamma1(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: gamma2(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: gamma3(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: gamma4(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: lambda(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: GAMMA(-nlevsno+1:0)                ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: mu_one                             ! two-stream coefficient from Toon et al. (lyr) [unitless]
    real:: e1(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr)
    real:: e2(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr)
    real:: e3(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr)
    real:: e4(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr)
    real:: C_pls_btm(-nlevsno+1:0)            ! intermediate variable: upward flux at bottom interface (lyr) [W/m2]
    real:: C_mns_btm(-nlevsno+1:0)            ! intermediate variable: downward flux at bottom interface (lyr) [W/m2]
    real:: C_pls_top(-nlevsno+1:0)            ! intermediate variable: upward flux at top interface (lyr) [W/m2]
    real:: C_mns_top(-nlevsno+1:0)            ! intermediate variable: downward flux at top interface (lyr) [W/m2]
    real:: A(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: B(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: D(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: E(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: AS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: DS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: X(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    real:: Y(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
    !-----------------------------------------------------------------------

   real :: h2osno(num_nourbanc)
   real :: snl(num_nourbanc)

    ! Enforce expected array sizes

   !  associate(&
   !       snl         =>   col_pp%snl     , & ! Input:  [integer (:)]  negative number of snow layers (col) [nbr]

   !       h2osno      =>   col_ws%h2osno        , & ! Input:  [real(r8) (:)]  snow liquid water equivalent (col) [kg/m2]
   !       frac_sno    =>   col_ws%frac_sno_eff    & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   !       )

   h2osno(1) = sum(h2osno_ice + h2osno_liq) ! EZSNOW
   snl(1) = -nlevsno

      ! Define constants
      ! pi = SHR_CONST_PI
      nint_snw_rds_min = nint(snw_rds_min)

      ! always use Delta approximation for snow
      DELTA = 1

      ! Get current timestep
      ! nstep = nstep_mod

      ! Loop over all non-urban columns
      ! (when called from CSIM, there is only one column)
      do fc = 1,num_nourbanc
         ! c_idx = filter_nourbanc(fc)
         c_idx = fc


         ! Zero absorbed radiative fluxes:
         do i=-nlevsno+1,1,1
            flx_abs_lcl(:,:)   = 0.0
            flx_abs(c_idx,i,:) = 0.0
         enddo

         ! set snow/ice mass to be used for RT:
         if (flg_snw_ice == 1) then
            h2osno_lcl = h2osno(c_idx)
         else
            h2osno_lcl = h2osno_ice(c_idx,0)
         endif


         ! Qualifier for computing snow RT:
         !  1) sunlight from atmosphere model
         !  2) minimum amount of snow on ground.
         !     Otherwise, set snow albedo to zero
         if ((coszen(c_idx) > 0.0) .and. (h2osno_lcl > min_snw)) then

            ! Set variables specific to CLM
            if (flg_snw_ice == 1) then
               ! If there is snow, but zero snow layers, we must create a layer locally.
               ! This layer is presumed to have the fresh snow effective radius.
               if (snl(c_idx) > -1) then
                  flg_nosnl         =  1
                  snl_lcl           =  -1
                  h2osno_ice_lcl(0) =  h2osno_lcl
                  h2osno_liq_lcl(0) =  0.0
                  snw_rds_lcl(0)    =  nint_snw_rds_min
               else
                  flg_nosnl         =  0
                  snl_lcl           =  snl(c_idx)
                  h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
                  h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
                  snw_rds_lcl(:)    =  snw_rds(c_idx,:)
               endif

               snl_btm   = 0
               snl_top   = snl_lcl+1

               ! for debugging only
               ! l_idx     = col_pp%landunit(c_idx)
               ! g_idx     = col_pp%gridcell(c_idx)
               ! sfctype   = lun_pp%itype(l_idx)
               ! lat_coord = grc_pp%latdeg(g_idx)
               ! lon_coord = grc_pp%londeg(g_idx)
               l_idx     = 1.0
               g_idx     = 1.0
               sfctype   = 1.0
               lat_coord = 1.0
               lon_coord = 1.0


               ! Set variables specific to CSIM
            else
               flg_nosnl         = 0
               snl_lcl           = -1
               h2osno_liq_lcl(:) = h2osno_liq(c_idx,:)
               h2osno_ice_lcl(:) = h2osno_ice(c_idx,:)
               snw_rds_lcl(:)    = snw_rds(c_idx,:)
               snl_btm           = 0
               snl_top           = 0
               sfctype           = -1
               lat_coord         = -90
               lon_coord         = 0
            endif



            ! Set local aerosol array
            do j=1,sno_nbr_aer
               mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
            enddo


            ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
            albsfc_lcl(1)                       = albsfc(c_idx,1)
            albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,2)


            ! Error check for snow grain size:
! #ifndef _OPENACC
            do i=snl_top,snl_btm,1
               if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                  write (*,*)  "SNICAR ERROR: snow grain radius of out of bounds."
                  write (*,*) "NSTEP= ", nstep
                  write (*,*) "flg_snw_ice= ", flg_snw_ice
                  write (*,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
                  write (*,*) "lat= ", lat_coord, " lon= ", lon_coord
                  write (*,*) "h2osno(c)= ", h2osno_lcl
                  ! call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
                  call land_error_message("SNICAR_RT in snicar_mod: snow grain radius out of bounds!", severity=FATAL)
               endif
            enddo
! #endif _OPENACC

            ! Incident flux weighting parameters
            !  - sum of all VIS bands must equal 1
            !  - sum of all NIR bands must equal 1
            !
            ! Spectral bands (5-band case)
            !  Band 1: 0.3-0.7um (VIS)
            !  Band 2: 0.7-1.0um (NIR)
            !  Band 3: 1.0-1.2um (NIR)
            !  Band 4: 1.2-1.5um (NIR)
            !  Band 5: 1.5-5.0um (NIR)
            !
            ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
            !
            ! 3-band weights
            if (numrad_snw==3) then
               ! Direct:
               if (flg_slr_in == 1) then
                  flx_wgt(1) = 1.
                  flx_wgt(2) = 0.66628670195247
                  flx_wgt(3) = 0.33371329804753
                  ! Diffuse:
               elseif (flg_slr_in == 2) then
                  flx_wgt(1) = 1.
                  flx_wgt(2) = 0.77887652162877
                  flx_wgt(3) = 0.22112347837123
               endif

               ! 5-band weights
            elseif(numrad_snw==5) then
               ! Direct:
               if (flg_slr_in == 1) then
                  flx_wgt(1) = 1.
                  flx_wgt(2) = 0.49352158521175
                  flx_wgt(3) = 0.18099494230665
                  flx_wgt(4) = 0.12094898498813
                  flx_wgt(5) = 0.20453448749347
                  ! Diffuse:
               elseif (flg_slr_in == 2) then
                  flx_wgt(1) = 1.
                  flx_wgt(2) = 0.58581507618433
                  flx_wgt(3) = 0.20156903770812
                  flx_wgt(4) = 0.10917889346386
                  flx_wgt(5) = 0.10343699264369
               endif
            endif

            ! Loop over snow spectral bands
            do bnd_idx = 1,numrad_snw

               mu_not    = coszen(c_idx)  ! must set here, because of error handling
               flg_dover = 1              ! default is to redo
               err_idx   = 0              ! number of times through loop

               do while (flg_dover > 0)

                  ! DEFAULT APPROXIMATIONS:
                  !  VIS:       Delta-Eddington
                  !  NIR (all): Delta-Hemispheric Mean
                  !  WARNING:   DO NOT USE DELTA-EDDINGTON FOR NIR DIFFUSE - this sometimes results in negative albedo
                  !
                  ! ERROR CONDITIONS:
                  !  Conditions which cause "trip", resulting in redo of RT approximation:
                  !   1. negative absorbed flux
                  !   2. total absorbed flux greater than incident flux
                  !   3. negative albedo
                  !   NOTE: These errors have only been encountered in spectral bands 4 and 5
                  !
                  ! ERROR HANDLING
                  !  1st error (flg_dover=2): switch approximation (Edd->HM or HM->Edd)
                  !  2nd error (flg_dover=3): change zenith angle by 0.02 (this happens about 1 in 10^6 cases)
                  !  3rd error (flg_dover=4): switch approximation with new zenith
                  !  Subsequent errors: repeatedly change zenith and approximations...

                  if (bnd_idx == 1) then
                     if (flg_dover == 2) then
                        APRX_TYP = 3
                     elseif (flg_dover == 3) then
                        APRX_TYP = 1
                        if (coszen(c_idx) > 0.5) then
                           mu_not = mu_not - 0.02
                        else
                           mu_not = mu_not + 0.02
                        endif
                     elseif (flg_dover == 4) then
                        APRX_TYP = 3
                     else
                        APRX_TYP = 1
                     endif

                  else
                     if (flg_dover == 2) then
                        APRX_TYP = 1
                     elseif (flg_dover == 3) then
                        APRX_TYP = 3
                        if (coszen(c_idx) > 0.5) then
                           mu_not = mu_not - 0.02
                        else
                           mu_not = mu_not + 0.02
                        endif
                     elseif (flg_dover == 4) then
                        APRX_TYP = 1
                     else
                        APRX_TYP = 3
                     endif

                  endif

                  ! Set direct or diffuse incident irradiance to 1
                  ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
                  if (flg_slr_in == 1) then
                     flx_slrd_lcl(bnd_idx) = 1./(mu_not*pi) ! this corresponds to incident irradiance of 1.0
                     flx_slri_lcl(bnd_idx) = 0.
                  else
                     flx_slrd_lcl(bnd_idx) = 0.
                     flx_slri_lcl(bnd_idx) = 1.
                  endif

                  ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                  ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                  if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
                     mss_cnc_aer_lcl(:,:) = 0.
                  endif

                  if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
                     mss_cnc_aer_lcl(:,:) = 0.
                  endif

                  ! Define local Mie parameters based on snow grain size and aerosol species,
                  !  retrieved from a lookup table.
                  if (flg_slr_in == 1) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (direct radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
                        asm_prm_snw_lcl(i)     = asm_prm_snw_drc(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
                     enddo
                  elseif (flg_slr_in == 2) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (diffuse radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
                        asm_prm_snw_lcl(i)     = asm_prm_snw_dfs(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
                     enddo
                  endif


                  ! aerosol species 1 optical properties
                 ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                 asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                 ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                  ! aerosol species 2 optical properties
                 ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                 asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                 ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)

                  ! aerosol species 3 optical properties
                  ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)
                  asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
                  ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

                  ! aerosol species 4 optical properties
                  ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)
                  asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
                  ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

                  ! aerosol species 5 optical properties
                  ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)
                  asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
                  ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)

                  ! aerosol species 6 optical properties
                  ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)
                  asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
                  ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)

                  ! aerosol species 7 optical properties
                  ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)
                  asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
                  ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)

                  ! aerosol species 8 optical properties
                  ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)
                  asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
                  ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)


                  ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                  ! 2. optical Depths (tau_snw, tau_aer)
                  ! 3. weighted Mie properties (tau, omega, g)

                  ! Weighted Mie parameters of each layer
                  do i=snl_top,snl_btm,1
                     L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                     tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                     do j=1,sno_nbr_aer
                        L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                        tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                     enddo

                     tau_sum   = 0.
                     omega_sum = 0.
                     g_sum     = 0.

                     do j=1,sno_nbr_aer
                        tau_sum    = tau_sum + tau_aer(i,j)
                        omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                        g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                     enddo

                     tau(i)    = tau_sum + tau_snw(i)
                     omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                     g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                  enddo

                  ! DELTA transformations, if requested
                  if (DELTA == 1) then
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)/(1+g(i))
                        omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                        tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
                     enddo
                  else
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)
                        omega_star(i) = omega(i)
                        tau_star(i)   = tau(i)
                     enddo
                  endif

                  ! Total column optical depth:
                  ! tau_elm(i) = total optical depth above the bottom of layer i
                  tau_elm(snl_top) = 0.
                  do i=snl_top+1,snl_btm,1
                     tau_elm(i) = tau_elm(i-1)+tau_star(i-1)
                  enddo

                  ! Direct radiation at bottom of snowpack:
                  F_direct_btm = albsfc_lcl(bnd_idx)*mu_not * &
                       exp(-(tau_elm(snl_btm)+tau_star(snl_btm))/mu_not)*pi*flx_slrd_lcl(bnd_idx)

                  ! Intermediates
                  ! Gamma values are approximation-specific.

                  ! Eddington
                  if (APRX_TYP==1) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = (7-(omega_star(i)*(4+(3*g_star(i)))))/4
                        gamma2(i) = -(1-(omega_star(i)*(4-(3*g_star(i)))))/4
                        gamma3(i) = (2-(3*g_star(i)*mu_not))/4
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 0.5
                     enddo

                     ! Quadrature
                  elseif (APRX_TYP==2) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = (3**0.5)*(2-(omega_star(i)*(1+g_star(i))))/2
                        gamma2(i) = omega_star(i)*(3**0.5)*(1-g_star(i))/2
                        gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 1/(3**0.5)
                     enddo

                     ! Hemispheric Mean
                  elseif (APRX_TYP==3) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = 2 - (omega_star(i)*(1+g_star(i)))
                        gamma2(i) = omega_star(i)*(1-g_star(i))
                        gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 0.5
                     enddo
                  endif

                  ! Intermediates for tri-diagonal solution
                  do i=snl_top,snl_btm,1
                     lambda(i) = sqrt(abs((gamma1(i)**2) - (gamma2(i)**2)))
                     GAMMA(i)  = gamma2(i)/(gamma1(i)+lambda(i))

                     e1(i)     = 1+(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                     e2(i)     = 1-(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                     e3(i)     = GAMMA(i) + exp(-lambda(i)*tau_star(i))
                     e4(i)     = GAMMA(i) - exp(-lambda(i)*tau_star(i))
                  enddo !enddo over snow layers


                  ! Intermediates for tri-diagonal solution
                  do i=snl_top,snl_btm,1
                     if (flg_slr_in == 1) then

                        C_pls_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-(tau_elm(i)+tau_star(i))/mu_not)*   &
                             (((gamma1(i)-(1/mu_not))*gamma3(i))+     &
                             (gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_mns_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-(tau_elm(i)+tau_star(i))/mu_not)*   &
                             (((gamma1(i)+(1/mu_not))*gamma4(i))+     &
                             (gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_pls_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-tau_elm(i)/mu_not)*(((gamma1(i)-(1/mu_not))* &
                             gamma3(i))+(gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_mns_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-tau_elm(i)/mu_not)*(((gamma1(i)+(1/mu_not))* &
                             gamma4(i))+(gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                     else
                        C_pls_btm(i) = 0.
                        C_mns_btm(i) = 0.
                        C_pls_top(i) = 0.
                        C_mns_top(i) = 0.
                     endif
                  enddo

                  ! Coefficients for tridiaganol matrix solution
                  do i=2*snl_lcl+1,0,1

                     !Boundary values for i=1 and i=2*snl_lcl, specifics for i=odd and i=even
                     if (i==(2*snl_lcl+1)) then
                        A(i) = 0
                        B(i) = e1(snl_top)
                        D(i) = -e2(snl_top)
                        E(i) = flx_slri_lcl(bnd_idx)-C_mns_top(snl_top)

                     elseif(i==0) then
                        A(i) = e1(snl_btm)-(albsfc_lcl(bnd_idx)*e3(snl_btm))
                        B(i) = e2(snl_btm)-(albsfc_lcl(bnd_idx)*e4(snl_btm))
                        D(i) = 0
                        E(i) = F_direct_btm-C_pls_btm(snl_btm)+(albsfc_lcl(bnd_idx)*C_mns_btm(snl_btm))

                     elseif(mod(i,2)==-1) then   ! If odd and i>=3 (n=1 for i=3)
                        n=floor(i/2.0)
                        A(i) = (e2(n)*e3(n))-(e4(n)*e1(n))
                        B(i) = (e1(n)*e1(n+1))-(e3(n)*e3(n+1))
                        D(i) = (e3(n)*e4(n+1))-(e1(n)*e2(n+1))
                        E(i) = (e3(n)*(C_pls_top(n+1)-C_pls_btm(n)))+(e1(n)*(C_mns_btm(n)-C_mns_top(n+1)))

                     elseif(mod(i,2)==0) then    ! If even and i<=2*snl_lcl
                        n=(i/2)
                        A(i) = (e2(n+1)*e1(n))-(e3(n)*e4(n+1))
                        B(i) = (e2(n)*e2(n+1))-(e4(n)*e4(n+1))
                        D(i) = (e1(n+1)*e4(n+1))-(e2(n+1)*e3(n+1))
                        E(i) = (e2(n+1)*(C_pls_top(n+1)-C_pls_btm(n)))+(e4(n+1)*(C_mns_top(n+1)-C_mns_btm(n)))
                     endif
                  enddo

                  AS(0) = A(0)/B(0)
                  DS(0) = E(0)/B(0)

                  do i=-1,(2*snl_lcl+1),-1
                     X(i)  = 1/(B(i)-(D(i)*AS(i+1)))
                     AS(i) = A(i)*X(i)
                     DS(i) = (E(i)-(D(i)*DS(i+1)))*X(i)
                  enddo

                  Y(2*snl_lcl+1) = DS(2*snl_lcl+1)
                  do i=(2*snl_lcl+2),0,1
                     Y(i) = DS(i)-(AS(i)*Y(i-1))
                  enddo

                  ! Downward direct-beam and net flux (F_net) at the base of each layer:
                  do i=snl_top,snl_btm,1
                     F_direct(i) = mu_not*pi*flx_slrd_lcl(bnd_idx)*exp(-(tau_elm(i)+tau_star(i))/mu_not)
                     F_net(i)    = (Y(2*i-1)*(e1(i)-e3(i))) + (Y(2*i)*(e2(i)-e4(i))) + &
                          C_pls_btm(i) - C_mns_btm(i) - F_direct(i)
                  enddo

                  ! Upward flux at snowpack top:
                  F_sfc_pls = (Y(2*snl_lcl+1)*(exp(-lambda(snl_top)*tau_star(snl_top))+ &
                       GAMMA(snl_top))) + (Y(2*snl_lcl+2)*(exp(-lambda(snl_top)* &
                       tau_star(snl_top))-GAMMA(snl_top))) + C_pls_top(snl_top)

                  ! Net flux at bottom = absorbed radiation by underlying surface:
                  F_btm_net = -F_net(snl_btm)


                  ! Bulk column albedo and surface net flux
                  albedo    = F_sfc_pls/((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))
                  F_sfc_net = F_sfc_pls - ((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))

                  trip = 0
                  ! Absorbed flux in each layer
                  do i=snl_top,snl_btm,1
                     if(i==snl_top) then
                        F_abs(i) = F_net(i)-F_sfc_net
                     else
                        F_abs(i) = F_net(i)-F_net(i-1)
                     endif
                     flx_abs_lcl(i,bnd_idx) = F_abs(i)


                     ! ERROR check: negative absorption
                     if (flx_abs_lcl(i,bnd_idx) < -0.00001) then
                        trip = 1
                     endif
                  enddo

                  flx_abs_lcl(1,bnd_idx) = F_btm_net

                  if (flg_nosnl == 1) then
                     ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                     !flx_abs_lcl(:,bnd_idx) = 0.
                     !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

                     ! changed on 20070408:
                     ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                     ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                     ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                     flx_abs_lcl(0,bnd_idx) = F_abs(0)
                     flx_abs_lcl(1,bnd_idx) = F_btm_net

                  endif

                  !Underflow check (we've already tripped the error condition above)
                  do i=snl_top,1,1
                     if (flx_abs_lcl(i,bnd_idx) < 0.) then
                        flx_abs_lcl(i,bnd_idx) = 0.
                     endif
                  enddo

                  F_abs_sum = 0.
                  do i=snl_top,snl_btm,1
                     F_abs_sum = F_abs_sum + F_abs(i)
                  enddo


                  !ERROR check: absorption greater than incident flux
                  ! (should make condition more generic than "1._r8")
                  if (F_abs_sum > 1.) then
                     trip = 1
                  endif

                  !ERROR check:
                  if ((albedo < 0.).and.(trip==0)) then
                     trip = 1
                  endif

                  ! Set conditions for redoing RT calculation
                  if ((trip == 1).and.(flg_dover == 1)) then
                     flg_dover = 2
                  elseif ((trip == 1).and.(flg_dover == 2)) then
                     flg_dover = 3
                  elseif ((trip == 1).and.(flg_dover == 3)) then
                     flg_dover = 4
                  elseif((trip == 1).and.(flg_dover == 4).and.(err_idx < 20)) then
                     flg_dover = 3
                     err_idx = err_idx + 1
                  elseif((trip == 1).and.(flg_dover == 4).and.(err_idx >= 20)) then
                     flg_dover = 0

                     write(*,*) "SNICAR ERROR: FOUND A WORMHOLE. STUCK IN INFINITE LOOP! Called from: ", flg_snw_ice
                     write(*,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)
                     write(*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)
                     write(*,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                     write(*,*) "SNICAR STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
                     write(*,*) "SNICAR STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
                     write(*,*) "SNICAR STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
                     write(*,*) "SNICAR STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
                     write(*,*) "SNICAR STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
                     write(*,*) "SNICAR STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
                     ! l_idx     = col_pp%landunit(c_idx)
                     ! write(*,*) "column index: ", c_idx
                     ! write(*,*) "landunit type", lun_pp%itype(l_idx)
                     ! write(*,*) "frac_sno: ", frac_sno(c_idx)
                  call land_error_message("SNICAR_RT in snicar_mod: FOUND A WORMHOLE. STUCK IN INFINITE LOOP!", severity=FATAL)

                  else
                     flg_dover = 0
                  endif

               enddo !enddo while (flg_dover > 0)

               ! Energy conservation check:
               ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
               energy_sum = (mu_not*pi*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
               if (abs(energy_sum) > 0.00001) then
                    write(*,*) "SNICAR ERROR: Energy conservation error of : ", energy_sum
                  call land_error_message("SNICAR_RT in snicar_mod: Energy conservation error!", severity=FATAL)
               endif

               albout_lcl(bnd_idx) = albedo

               ! Check that albedo is less than 1
               if (albout_lcl(bnd_idx) > 1.0) then
                  write(*,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx
                  write(*,*) "SNICAR STATS: bnd_idx= ",bnd_idx
                  write (*,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), &
                       " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
                  write (*,*) "SNICAR STATS: landtype= ", sfctype
                  write (*,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                  write (*,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in

                  write (*,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
                  write (*,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
                  write (*,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
                  write (*,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
                  write (*,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)

                  write (*,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
                  write (*,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
                  write (*,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
                  write (*,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
                  write (*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

                  write (*,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
                  write (*,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
                  write (*,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
                  write (*,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
                  write (*,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)
                  call land_error_message("SNICAR_RT in snicar_mod: Snow albedo out of bounds!", severity=FATAL)
               endif

            enddo   ! loop over wvl bands


            ! Weight output NIR albedo appropriately
            albout(c_idx,1) = albout_lcl(1)
            flx_sum         = 0.
            do bnd_idx= nir_bnd_bgn,nir_bnd_end
               flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
            end do
            albout(c_idx,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

            ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
            flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
            do i=snl_top,1,1
               flx_sum = 0.
               do bnd_idx= nir_bnd_bgn,nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
               enddo
               flx_abs(c_idx,i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
            end do

            ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
         elseif ( (coszen(c_idx) > 0.) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0.) ) then
            albout(c_idx,1) = albsfc(c_idx,1)
            albout(c_idx,2) = albsfc(c_idx,2)

            ! There is either zero snow, or no sun
         else
            albout(c_idx,1) = 0.
            albout(c_idx,2) = 0.
         endif    ! if column has snow and coszen > 0

      enddo    ! loop over all columns

   !  end associate

  end subroutine SNICAR_RT


  !-----------------------------------------------------------------------
     subroutine read_snicar_optics_data()

     integer :: ier      
     logical :: readvar      
      ! write(*,*) "Reading SNICAR optics data"
      ! LM4p2 READ:
      ! direct-beam snow Mie parameters:
      call read_data( ncid, 'ss_alb_ice_drc', ss_alb_snw_drc, no_domain=.true.)
      call read_data( ncid, 'asm_prm_ice_drc', asm_prm_snw_drc, no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc, no_domain=.true.)
      ! diffuse snow Mie parameters:
      call read_data( ncid, 'ss_alb_ice_dfs', ss_alb_snw_dfs,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_ice_dfs', asm_prm_snw_dfs,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs, no_domain=.true.)

      if (snicar_atm_type > 0)then
         call read_data( ncid, 'flx_wgt_dir', flx_wgt_dir, no_domain=.true.) ! direct-beam incident spectral flux: 
         call read_data( ncid, 'flx_wgt_dif', flx_wgt_dif, no_domain=.true.) ! diffuse incident spectral flux:
      endif

      ! BC species 1 Mie parameters
      call read_data( ncid, 'ss_alb_bcphil', ss_alb_bc1,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_bcphil', asm_prm_bc1,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_bcphil', ext_cff_mss_bc1, no_domain=.true.)
      ! ! BC species 2 Mie parameters
      call read_data( ncid, 'ss_alb_bcphob', ss_alb_bc2,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_bcphob', asm_prm_bc2,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_bcphob', ext_cff_mss_bc2, no_domain=.true.)

      ! OC species 1 Mie parameters
      call read_data( ncid, 'ss_alb_ocphil',      ss_alb_oc1,      no_domain=.true.)
      call read_data( ncid, 'asm_prm_ocphil',     asm_prm_oc1,     no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_ocphil', ext_cff_mss_oc1, no_domain=.true.)
      !
      ! OC species 2 Mie parameters
      call read_data( ncid, 'ss_alb_ocphob', ss_alb_oc2,            no_domain=.true.)
      call read_data( ncid, 'asm_prm_ocphob', asm_prm_oc2,          no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_ocphob', ext_cff_mss_oc2,  no_domain=.true.)
      !
      ! dust species 1 Mie parameters
      call read_data( ncid, 'ss_alb_dust01', ss_alb_dst1,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_dust01', asm_prm_dst1,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_dust01', ext_cff_mss_dst1, no_domain=.true.)
      !
      ! dust species 2 Mie parameters
      call read_data( ncid, 'ss_alb_dust02', ss_alb_dst2,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_dust02', asm_prm_dst2,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_dust02', ext_cff_mss_dst2, no_domain=.true.)
      !
      ! dust species 3 Mie parameters
      call read_data( ncid, 'ss_alb_dust03', ss_alb_dst3,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_dust03', asm_prm_dst3,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_dust03', ext_cff_mss_dst3, no_domain=.true.)
      !
      ! dust species 4 Mie parameters
      call read_data( ncid, 'ss_alb_dust04', ss_alb_dst4,           no_domain=.true.)
      call read_data( ncid, 'asm_prm_dust04', asm_prm_dst4,         no_domain=.true.)
      call read_data( ncid, 'ext_cff_mss_dust04', ext_cff_mss_dst4, no_domain=.true.)



      if(is_watch_point()) then
         write(*,*) "##### read_snicar_optics_data - start checkpoint 1 #####"

         write(*,*) "direct-beam snow Mie parameters:"
         write(*,*) 'ss_alb_ice_drc', ss_alb_snw_drc
         write(*,*) 'asm_prm_ice_drc', asm_prm_snw_drc
         write(*,*) 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc
         write(*,*) "diffuse snow Mie parameters:"
         write(*,*)'ss_alb_ice_dfs', ss_alb_snw_dfs           
         write(*,*)'asm_prm_ice_dfs', asm_prm_snw_dfs         
         write(*,*)'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs 

      if (snicar_atm_type > 0)then
         write(*,*) "Solar spectrum weights:"
         write(*,*) "flx_wgt_dir", flx_wgt_dir
         write(*,*) "flx_wgt_dif", flx_wgt_dif
      endif
         write(*,*) "##### read_snicar_optics_data - end checkpoint 1 #####"
     endif
      !
      !    write(*,*) 'Successfully read snow optical properties'
      !    ! print some diagnostics:
      !    write (*,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
      !         ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
      !         ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
      !    write (*,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
      !         ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
      !         ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
      !    if (DO_SNO_OC) then
      !       write (*,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
      !    else
      !       write (*,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
      !    endif

      !    write (*,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
      !         ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
      !    write (*,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
      !         ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
      !   !
      !    if (DO_SNO_OC) then
      !       write (*,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
      !            ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
      !       write (*,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
      !            ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
      !    endif
      !    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
      !         ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
      !    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
      !         ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
      !    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
      !         ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
      !    write (*,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
      !         ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
      !    write(*,*)
      !
    end subroutine read_snicar_optics_data



   !-----------------------------------------------------------------------
   subroutine SNICAR_AD_RT (nlevsno, flg_snw_ice,  &
                         coszen, flg_slr_in, h2osno_liq, h2osno_ice, snw_rds, snw_shp,   &
                         mss_cnc_aer_in, albsfc, albout, flx_abs)
     !
     ! !DESCRIPTION:
     ! Determine reflectance of, and vertically-resolved solar absorption in,
     ! snow with impurities, with updated shortwave scheme
     !
     ! The multi-layer solution for multiple-scattering used here is from:
     ! Briegleb, P. and Light, B.: A Delta-Eddington mutiple scattering
     ! parameterization for solar radiation in the sea ice component of the
     ! community climate system model, 2007.
     !
     ! The implementation of the SNICAR-AD model in ELM is described in:
     ! Dang et al., Inter-comparison and improvement of 2-stream shortwave
     ! radiative transfer models for unified treatment of cryospheric surfaces
     ! in ESMs, in review, 2019
 
     integer       , intent(in)  :: nlevsno    ! number of layer in snowpack
     integer       , intent(in)  :: flg_snw_ice                                        ! flag: =1 when called from CLM, =2 when called from CSIM
     real          , intent(in)  :: coszen         ( 1: )                    ! cosine of solar zenith angle for next time step (col) [unitless]
     integer       , intent(in)  :: flg_slr_in                                         ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux
     real          , intent(in)  :: h2osno_liq     ( 1: , -nlevsno+1: )      ! liquid water content (col,lyr) [kg/m2]
     real          , intent(in)  :: h2osno_ice     ( 1: , -nlevsno+1: )      ! ice content (col,lyr) [kg/m2]
     integer       , intent(in)  :: snw_rds        ( 1: , -nlevsno+1: )      ! snow effective radius (col,lyr) [microns, m^-6]
     integer       , intent(in)  :: snw_shp        ( 1: , -nlevsno+1: )      ! snow shape layer by layer (col,lyr) [nbr]
     real          , intent(in)  :: mss_cnc_aer_in ( 1: , -nlevsno+1: , 1: ) ! mass concentration of all aerosol species (col,lyr,aer) [kg/kg]
     real          , intent(in)  :: albsfc         ( 1: , 1: )               ! albedo of surface underlying snow (col,bnd) [frc]
     real          , intent(out) :: albout         ( 1: , 1: )               ! snow albedo, averaged into 2 bands (=0 if no sun or no snow) (col,bnd) [frc]
     real          , intent(out) :: flx_abs        ( 1: , -nlevsno+1: , 1: ) ! absorbed flux in each layer per unit flux incident (col, lyr, bnd)
     !
     ! !LOCAL VARIABLES:

     ! Local variables representing single-column values of arrays:
     integer :: snl_lcl                            ! negative number of snow layers [nbr]
     integer :: snw_rds_lcl(-nlevsno+1:0)          ! snow effective radius [m^-6]
     real:: flx_slrd_lcl(1:numrad_snw)         ! direct beam incident irradiance [W/m2] (set to 1)
     real:: flx_slri_lcl(1:numrad_snw)         ! diffuse incident irradiance [W/m2] (set to 1)
     real:: mss_cnc_aer_lcl(-nlevsno+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
     real:: h2osno_lcl                         ! total column snow mass [kg/m2]
     real:: h2osno_liq_lcl(-nlevsno+1:0)       ! liquid water mass [kg/m2]
     real:: h2osno_ice_lcl(-nlevsno+1:0)       ! ice mass [kg/m2]
     real:: albsfc_lcl(1:numrad_snw)           ! albedo of underlying surface [frc]
     real:: ss_alb_snw_lcl(-nlevsno+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
     real:: asm_prm_snw_lcl(-nlevsno+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
     real:: ext_cff_mss_snw_lcl(-nlevsno+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
     real:: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc]
     real:: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
     real:: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]

     ! Other local variables
     integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
                                                   ! (1= use, 0= don't use)
     real:: flx_wgt(1:numrad_snw)              ! weights applied to spectral bands,
                                                   ! specific to direct and diffuse cases (bnd) [frc]
     integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
                                                   ! =0 if at least 1 snow layer [flg]
     !integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
     !integer :: flg_dover                          ! defines conditions for RT redo (explained below)

     real:: albedo                             ! temporary snow albedo [frc]
     real:: flx_sum                            ! temporary summation variable for NIR weighting
     real:: albout_lcl(numrad_snw)             ! snow albedo by band [frc]
     real:: flx_abs_lcl(-nlevsno+1:1,numrad_snw)! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]

     real:: L_snw(-nlevsno+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
     real:: tau_snw(-nlevsno+1:0)              ! snow optical depth (lyr) [unitless]
     real:: L_aer(-nlevsno+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2]
     real:: tau_aer(-nlevsno+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
     real:: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
     real:: tau_elm(-nlevsno+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless]
     real:: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
     real:: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]

     real:: tau(-nlevsno+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
     real:: omega(-nlevsno+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
     real:: g(-nlevsno+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
     real:: tau_star(-nlevsno+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                                                   ! (lyr) [unitless]
     real:: omega_star(-nlevsno+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
     real:: g_star(-nlevsno+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                                                   ! (lyr) [frc]

     integer :: nstep                              ! current timestep [nbr] (debugging only)
     integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
     integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= numrad_snw) [idx]
     integer :: rds_idx                            ! snow effective radius index for retrieving
                                                   ! Mie parameters from lookup table [idx]
     integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
     integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
     integer :: fc                                 ! column filter index
     integer :: i                                  ! layer index [idx]
     integer :: j                                  ! aerosol number index [idx]
     integer :: m                                  ! secondary layer index [idx]
     integer :: nint_snw_rds_min                   ! nearest integer value of snw_rds_min

     real:: F_abs(-nlevsno+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
     real:: F_abs_sum                          ! total absorbed energy in column [W/m^2]
     real:: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
     real:: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]
     real:: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
     real:: mu_not                             ! cosine of solar zenith angle (used locally) [frc]

     integer :: err_idx                            ! counter for number of times through error loop [nbr]
     real:: lat_coord                          ! gridcell latitude (debugging only)
     real:: lon_coord                          ! gridcell longitude (debugging only)
     integer :: sfctype                            ! underlying surface type (debugging only)
   !   real(r8):: pi                                 ! 3.1415...

   !!!!!!!!!!!!!!!!!!!!!!!!
     ! New variales for non-spherical snow shape   ! He et al., 2017
     integer :: snw_shp_lcl(-nlevsno+1:0)          ! Snow grain shape option:
                                                   ! 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake
     integer :: snw_fs_lcl(-nlevsno+1:0)           ! Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                                   ! 0=use recommended default value (He et al. 2017);
                                                   ! others(0<fs<1)= use user-specified value
                                                   ! only activated when sno_shp > 1 (i.e. nonspherical)
     integer :: snw_ar_lcl(-nlevsno+1:0)           ! % Aspect ratio: ratio of grain width to length
                                                   ! 0=use recommended default value (He et al. 2017);
                                                   ! others(0.1<fs<20)= use user-specified value
                                                   ! only activated when sno_shp > 1 (i.e. nonspherical)
     real:: &
         diam_ice  , & !
         fs_sphd  , & !
         fs_hex0  , & ! 
         fs_hex  , & ! 
         fs_koch  , & ! 
         AR_tmp  , & ! 
         g_ice_Cg_tmp(7)  , & !
         gg_ice_F07_tmp(7)  , & !
         g_ice_F07  , & !
         g_ice  , & !
         gg_F07_intp  , & !
         g_Cg_intp, & !
         R_1_omega_tmp , & ! !!! BC internal mixing
         C_BC_total , & ! !!! BC concentrtion
         C_dust_total !! dust concentration
     integer :: slr_zen
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
     ! SNICAR_AD new variables, follow sea-ice shortwave conventions
     real:: &
        trndir(-nlevsno+1:1)  , & ! solar beam down transmission from top
        trntdr(-nlevsno+1:1)  , & ! total transmission to direct beam for layers above
        trndif(-nlevsno+1:1)  , & ! diffuse transmission to diffuse beam for layers above
        rupdir(-nlevsno+1:1)  , & ! reflectivity to direct radiation for layers below
        rupdif(-nlevsno+1:1)  , & ! reflectivity to diffuse radiation for layers below
        rdndif(-nlevsno+1:1)  , & ! reflectivity to diffuse radiation for layers above
        dfdir(-nlevsno+1:1)   , & ! down-up flux at interface due to direct beam at top surface
        dfdif(-nlevsno+1:1)   , & ! down-up flux at interface due to diffuse beam at top surface
        dftmp(-nlevsno+1:1)       ! temporary variable for down-up flux at interface

     real:: &
        rdir(-nlevsno+1:0)       , & ! layer reflectivity to direct radiation
        rdif_a(-nlevsno+1:0)     , & ! layer reflectivity to diffuse radiation from above
        rdif_b(-nlevsno+1:0)     , & ! layer reflectivity to diffuse radiation from below
        tdir(-nlevsno+1:0)       , & ! layer transmission to direct radiation (solar beam + diffuse)
        tdif_a(-nlevsno+1:0)     , & ! layer transmission to diffuse radiation from above
        tdif_b(-nlevsno+1:0)     , & ! layer transmission to diffuse radiation from below
        trnlay(-nlevsno+1:0)         ! solar beam transm for layer (direct beam only)

     real:: &
         ts       , & ! layer delta-scaled extinction optical depth
         ws       , & ! layer delta-scaled single scattering albedo
         gs       , & ! layer delta-scaled asymmetry parameter
         extins   , & ! extinction
         alp      , & ! temporary for alpha
         gam      , & ! temporary for agamm
         amg      , & ! alp - gam
         apg      , & ! alp + gam
         ue       , & ! temporary for u
         refk     , & ! interface multiple scattering
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref
         tdndif       ! total down diffuse = tot tran - direct tran

     real :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         agamm    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,agamm,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         mu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n

     ! perpendicular and parallel relative to plane of incidence and scattering
     real :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below

     real :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt      , & ! accumulator for tdif gaussian integration
         exp_min      ! minimum exponential value

     integer :: &
         ng             , & ! gaussian integration index
         snl_btm_itf    , & ! index of bottom snow layer interfaces (1) [idx]
         ngmax = 8          ! gaussian integration index

     ! Gaussian integration angle and coefficients
     real :: &
         difgauspt(1:8)  , &
         difgauswt(1:8)
     ! real(r8),  dimension (1:8) :: &
     !     dif_gauspt     & ! gaussian angles (radians)
     !       = (/ 0.9894009_r8,  0.9445750_r8, &
     !            0.8656312_r8,  0.7554044_r8, &
     !            0.6178762_r8,  0.4580168_r8, &
     !            0.2816036_r8,  0.0950125_r8/) , &
     !     dif_gauswt     & ! gaussian weights
     !       = (/ 0.0271525_r8,  0.0622535_r8, &
     !            0.0951585_r8,  0.1246290_r8, &
     !            0.1495960_r8,  0.1691565_r8, &
     !            0.1826034_r8,  0.1894506_r8/)

     ! constants used in algorithm
     real :: &
         c0      = 0.0     , &
         c1      = 1.0     , &
         c3      = 3.0     , &
         c4      = 4.0     , &
         c6      = 6.0     , &
         cp01    = 0.01    , &
         cp5     = 0.5     , &
         cp75    = 0.75    , &
         c1p5    = 1.5     , &
         trmin   = 0.001   , &
         argmax  = 10.0       ! maximum argument of exponential

     ! cconstant coefficients used for SZA parameterization
     real :: &
         sza_a0 =  0.085730 , &
         sza_a1 = -0.630883 , &
         sza_a2 =  1.303723 , &
         sza_b0 =  1.467291 , &
         sza_b1 = -3.338043 , &
         sza_b2 =  6.807489 , &
         puny   =  1.0e-11  , &
         mu_75  =  0.2588       ! cosine of 75 degree

     ! coefficients used for SZA parameterization
     real :: &
         sza_c1          , & ! coefficient, SZA parameteirzation
         sza_c0          , & ! coefficient, SZA parameterization
         sza_factor      , & ! factor used to adjust NIR direct albedo
         flx_sza_adjust  , & ! direct NIR flux adjustment from sza_factor
         mu0                 ! incident solar zenith angle

     ! Delta-Eddington solution expressions
     ! alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
     ! agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
     ! n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
     ! u(w,gg,e)        = c1p5*(c1 - w*gg)/e
     ! el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))

     !-----------------------------------------------------------------------


 ! Constants for aspherical ice particles %%%
      ! g_snw asymmetry factor parameterization coefficients (6 bands) from
      !  Table 3 & Eqs. 6-7 in He et al. (2017)
      ! assume same values for 4-5 um band, which leads to very small biases (<3%)
      
      real :: g_b2(7)
      real :: g_b1(7)
      real :: g_b0(7)
      real :: g_F07_c2(7)
      real :: g_F07_c1(7)
      real :: g_F07_c0(7)
      real :: g_F07_p2(7)
      real :: g_F07_p1(7)
      real :: g_F07_p0(7)
      real :: BC_d0(3)
      real :: BC_d1(3)
      real :: BC_d2(3)
      real :: dust_clear_d0(3)
      real :: dust_clear_d1(3)
      real :: dust_clear_d2(3)
      real :: dust_cloudy_d0(3)
      real :: dust_cloudy_d1(3)
      real :: dust_cloudy_d2(3)

   real :: h2osno(num_nourbanc)
   real :: snl(num_nourbanc)
      
     ! Enforce expected array sizes

   !   associate(&
   !        snl         =>   col_pp%snl           , & ! Input:  [integer (:)]  negative number of snow layers (col) [nbr]
   !        h2osno      =>   col_ws%h2osno        , & ! Input:  [real(r8) (:)]  snow liquid water equivalent (col) [kg/m2]
   !        frac_sno    =>   col_ws%frac_sno_eff    & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   !        )


   h2osno(1) = sum(h2osno_ice + h2osno_liq) ! EZSNOW
   snl(1) = -nlevsno

       ! Define constants
      !  pi = SHR_CONST_PI
       nint_snw_rds_min = nint(snw_rds_min)

       ! always use Delta approximation for snow
       DELTA = 1

       ! Get current timestep
      !  nstep = get_nstep()
       nstep = 1 ! NOT USED in LM4p2

       !Gaussian integration angle and coefficients for diffuse radiation
       difgauspt(1:8)     & ! gaussian angles (radians)
         = (/ 0.9894009,  0.9445750, &
              0.8656312,  0.7554044, &
              0.6178762,  0.4580168, &
              0.2816036,  0.0950125/)
       difgauswt(1:8)     & ! gaussian weights
         = (/ 0.0271525,  0.0622535, &
              0.0951585,  0.1246290, &
              0.1495960,  0.1691565, &
              0.1826034,  0.1894506/)
 !!!!!!!!!!!! snow albedo improvement by Dalei Hao 2021
      !!!!!!!!! snow shape
      ! define snow shape
      ! snw_shp_lcl(:) = snow_shape_defined
      ! EZDEV - read level-by-level snow shape:
      snw_shp_lcl(:) = snw_shp(1,:)
      snw_fs_lcl(:)  = 0. 
      snw_ar_lcl(:)  = 0. 
      !data g_wvl(:) /0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00/ ! wavelength (um) division point
      !g_wvl_center = g_wvl(2:8)/2 + g_wvl(1:7)/2 ; ! center point for wavelength band
      data g_b0(:) /9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01/
      data g_b1(:) /5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02/
      data g_b2(:) /-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04/
      ! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007
      data g_F07_c2(:) /1.349959E-1,1.115697E-1,9.853958E-2,5.557793E-2,-1.233493E-1,0.0,0.0/
      data g_F07_c1(:) /-3.987320E-1,-3.723287E-1,-3.924784E-1,-3.259404E-1,4.429054E-2,-1.726586E-1,-1.726586E-1/
      data g_F07_c0(:) /7.938904E-1,8.030084E-1,8.513932E-1,8.692241E-1,7.085850E-1,6.412701E-1,6.412701E-1/
      data g_F07_p2(:) /3.165543E-3,2.014810E-3,1.780838E-3,6.987734E-4,-1.882932E-2,-2.277872E-2,-2.277872E-2/
      data g_F07_p1(:) /1.140557E-1,1.143152E-1,1.143814E-1,1.071238E-1,1.353873E-1,1.914431E-1,1.914431E-1/
      data g_F07_p0(:) /5.292852E-1,5.425909E-1,5.601598E-1,6.023407E-1,6.473899E-1,4.634944E-1,4.634944E-1/
      
      !!! BC internal mixing
      data BC_d0(:) /3.50098E-02,6.51688E-03,7.96544E-01/
      data BC_d1(:) /9.91050E-01,7.36315E-01,4.36649E-02/
      data BC_d2(:) /3.00370E+01,9.52134E+02,2.57288E+02/
      
      !!! dust internal mixing
      data dust_clear_d0(:) /1.0413E+00,1.0168E+00,1.0189E+00/
      data dust_clear_d1(:) /1.0016E+00,1.0070E+00,1.0840E+00/
      data dust_clear_d2(:) /2.4208E-01,1.5300E-03,1.1230E-04/
      
      data dust_cloudy_d0(:) /1.0388E+00,1.0167E+00,1.0189E+00/
      data dust_cloudy_d1(:) /1.0015E+00,1.0061E+00,1.0823E+00/
      data dust_cloudy_d2(:) /2.5973E-01,1.6200E-03,1.1721E-04/
!!!!!!!!!!!!!


      
      ! Loop over all non-urban columns
      ! (when called from CSIM, there is only one column)
       do fc = 1,num_nourbanc
         !  c_idx = filter_nourbanc(fc)
         c_idx = fc

          ! Zero absorbed radiative fluxes:
          do i=-nlevsno+1,1,1
             flx_abs_lcl(:,:)   = 0.
             flx_abs(c_idx,i,:) = 0.
          enddo

          ! set snow/ice mass to be used for RT:
          if (flg_snw_ice == 1) then
             h2osno_lcl = h2osno(c_idx)
          else
             h2osno_lcl = h2osno_ice(c_idx,0)
          endif


          ! Qualifier for computing snow RT:
          !  1) sunlight from atmosphere model
          !  2) minimum amount of snow on ground.
          !     Otherwise, set snow albedo to zero
          if ((coszen(c_idx) > 0.) .and. (h2osno_lcl > min_snw) ) then

             ! Set variables specific to ELM
             if (flg_snw_ice == 1) then
                ! If there is snow, but zero snow layers, we must create a layer locally.
                ! This layer is presumed to have the fresh snow effective radius.
                if (snl(c_idx) > -1) then
                   flg_nosnl         =  1
                   snl_lcl           =  -1
                   h2osno_ice_lcl(0) =  h2osno_lcl
                   h2osno_liq_lcl(0) =  0.
                   snw_rds_lcl(0)    =  nint_snw_rds_min
                else
                   flg_nosnl         =  0
                   snl_lcl           =  snl(c_idx)
                   h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
                   h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
                   snw_rds_lcl(:)    =  snw_rds(c_idx,:)
                endif

                snl_btm   = 0
                snl_top   = snl_lcl+1

                ! for debugging only
               !  l_idx     = col_pp%landunit(c_idx)
               !  g_idx     = col_pp%gridcell(c_idx)
               !  sfctype   = lun_pp%itype(l_idx)
               !  lat_coord = grc_pp%latdeg(g_idx)
               !  lon_coord = grc_pp%londeg(g_idx)
                l_idx     = 1.0
                g_idx     = 1.0
                sfctype   = 1.0
                lat_coord = 1.0
                lon_coord = 1.0


                ! Set variables specific to CSIM
             else
                flg_nosnl         = 0
                snl_lcl           = -1
                h2osno_liq_lcl(:) = h2osno_liq(c_idx,:)
                h2osno_ice_lcl(:) = h2osno_ice(c_idx,:)
                snw_rds_lcl(:)    = snw_rds(c_idx,:)
                snl_btm           = 0
                snl_top           = 0
                sfctype           = -1
                lat_coord         = -90
                lon_coord         = 0
             endif ! end if flg_snw_ice == 1

             ! Set local aerosol array
             do j=1,sno_nbr_aer
                mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
             enddo


             ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
             albsfc_lcl(1)                       = albsfc(c_idx,1)
             albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,2)


             ! Error check for snow grain size:
             do i=snl_top,snl_btm,1
                if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                   write (*,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
                  !  write (*,*) "NSTEP= ", nstep
                   write (*,*) "flg_snw_ice= ", flg_snw_ice
                   write (*,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
                   write (*,*) "lat= ", lat_coord, " lon= ", lon_coord
                   write (*,*) "h2osno(c)= ", h2osno_lcl
                  call land_error_message("SNICAR_AD_RT in snicar_mod: Snow grain radius out of bounds!", severity=FATAL)
                endif
             enddo

             ! Incident flux weighting parameters
             !  - sum of all VIS bands must equal 1
             !  - sum of all NIR bands must equal 1
             !
             ! Spectral bands (5-band case)
             !  Band 1: 0.3-0.7um (VIS)
             !  Band 2: 0.7-1.0um (NIR)
             !  Band 3: 1.0-1.2um (NIR)
             !  Band 4: 1.2-1.5um (NIR)
             !  Band 5: 1.5-5.0um (NIR)
             !
             ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
             !
             ! 3-band weights
             if (numrad_snw==3) then
                ! Direct:
                if (flg_slr_in == 1) then
                   flx_wgt(1) = 1.
                   flx_wgt(2) = 0.66628670195247
                   flx_wgt(3) = 0.33371329804753
                   ! Diffuse:
                elseif (flg_slr_in == 2) then
                   flx_wgt(1) = 1.
                   flx_wgt(2) = 0.77887652162877
                   flx_wgt(3) = 0.22112347837123
                endif

                ! 5-band weights
             elseif(numrad_snw==5) then
                ! Direct:
                if (flg_slr_in == 1) then
                   if (snicar_atm_type == 0) then
                     flx_wgt(1) = 1.
                     flx_wgt(2) = 0.49352158521175
                     flx_wgt(3) = 0.18099494230665
                     flx_wgt(4) = 0.12094898498813
                     flx_wgt(5) = 0.20453448749347
                  else                 
                     slr_zen = nint(acosd(coszen(c_idx)))
                     if (slr_zen>89) then
                        slr_zen = 89
                     endif
                     flx_wgt(1) = 1.
                     flx_wgt(2) = flx_wgt_dir(snicar_atm_type, slr_zen+1, 2)
                     flx_wgt(3) = flx_wgt_dir(snicar_atm_type, slr_zen+1, 3)
                     flx_wgt(4) = flx_wgt_dir(snicar_atm_type, slr_zen+1, 4)
                     flx_wgt(5) = flx_wgt_dir(snicar_atm_type, slr_zen+1, 5)  
                     
                    ! write(iulog,*) "SNICAR_AD STATS: coszen(c_idx) (0)= ", coszen(c_idx) ! add by Dalei check
                   !  write(iulog,*) "SNICAR_AD STATS: slr_zen (0)= ", slr_zen ! add by Dalei check
                    ! write(iulog,*) "SNICAR_AD STATS: flx_wgt (2)= ", flx_wgt(2) ! add by Dalei check
                    ! write(iulog,*) "SNICAR_AD STATS: flx_wgt (4)= ", flx_wgt(4) ! add by Dalei check
                  endif
                   ! Diffuse:
                elseif (flg_slr_in == 2) then
                     if  (snicar_atm_type == 0) then
                     flx_wgt(1) = 1.
                     flx_wgt(2) = 0.58581507618433
                     flx_wgt(3) = 0.20156903770812
                     flx_wgt(4) = 0.10917889346386
                     flx_wgt(5) = 0.10343699264369
                  else
                     flx_wgt(1) = 1.
                     flx_wgt(2) = flx_wgt_dif(snicar_atm_type, 2)
                     flx_wgt(3) = flx_wgt_dif(snicar_atm_type, 3)
                     flx_wgt(4) = flx_wgt_dif(snicar_atm_type, 4)
                     flx_wgt(5) = flx_wgt_dif(snicar_atm_type, 5)
                     
                     !write(iulog,*) "SNICAR_AD STATS: flx_wgt (0)= ", flx_wgt(2) ! add by Dalei check
                     !write(iulog,*) "SNICAR_AD STATS: flx_wgt (0)= ", flx_wgt(4) ! add by Dalei check
                  endif
                endif
             endif ! end if numrad_snw

             ! Loop over snow spectral bands

             exp_min = exp(-argmax)
             do bnd_idx = 1,numrad_snw

               ! note that we can remove flg_dover since this algorithm is
               ! stable for mu_not > 0.01

               ! mu_not is cosine solar zenith angle above the fresnel level; make
               ! sure mu_not is large enough for stable and meaningful radiation
               ! solution: .01 is like sun just touching horizon with its lower edge
               ! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
                mu_not = max(coszen(c_idx), cp01)


                   ! Set direct or diffuse incident irradiance to 1
                   ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
                   if (flg_slr_in == 1) then
                      flx_slrd_lcl(bnd_idx) = 1./(mu_not*pi) ! this corresponds to incident irradiance of 1.0
                      flx_slri_lcl(bnd_idx) = 0.
                   else
                      flx_slrd_lcl(bnd_idx) = 0.
                      flx_slri_lcl(bnd_idx) = 1.
                   endif

                   ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                   ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                   if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
                      mss_cnc_aer_lcl(:,:) = 0.
                   endif

                   if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
                      mss_cnc_aer_lcl(:,:) = 0.
                   endif

                   ! Define local Mie parameters based on snow grain size and aerosol species,
                   !  retrieved from a lookup table.
                   if (flg_slr_in == 1) then
                      do i=snl_top,snl_btm,1
                         rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                         ! snow optical properties (direct radiation)
                         ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
                         asm_prm_snw_lcl(i)     = asm_prm_snw_drc(rds_idx,bnd_idx)
                         ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
                      enddo
                   elseif (flg_slr_in == 2) then
                      do i=snl_top,snl_btm,1
                         rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                         ! snow optical properties (diffuse radiation)
                         ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
                         asm_prm_snw_lcl(i)     = asm_prm_snw_dfs(rds_idx,bnd_idx)
                         ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
                      enddo
                   endif

  !!! Dalei Hao 
                  ! shape-dependent asymetry factors (He et al., 2017)
                  do i=snl_top,snl_btm,1
                     if(snw_shp_lcl(i) == 2) then ! spheroid
                     
                       diam_ice = 2.*snw_rds_lcl(i)
                        if(snw_fs_lcl(i) == 0) then
                           fs_sphd = 0.929
                        else
                           fs_sphd = snw_fs_lcl(i)               
                        endif
                        fs_hex = 0.788 
                        if(snw_ar_lcl(i) == 0) then
                           AR_tmp = 0.5
                        else
                           AR_tmp = snw_ar_lcl(i)              
                        endif
                        g_ice_Cg_tmp = g_b0 * ((fs_sphd/fs_hex)**g_b1) * (diam_ice**g_b2) ! Eq.7, He et al. (2017)
                        gg_ice_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * (AR_tmp**2) ! Eqn. 3.1 in Fu (2007)                           
            
                     elseif(snw_shp_lcl(i) == 3) then ! hexagonal plate
                          diam_ice = 2.*snw_rds_lcl(i)
                        if(snw_fs_lcl(i) == 0) then
                           fs_hex0 = 0.788
                        else
                           fs_hex0 = snw_fs_lcl(i)               
                        endif
                        fs_hex = 0.788 
                        if(snw_ar_lcl(i) == 0) then
                           AR_tmp = 2.5
                        else
                           AR_tmp = snw_ar_lcl(i)              
                        endif
                        g_ice_Cg_tmp = g_b0 * ((fs_hex0/fs_hex)**g_b1) * (diam_ice**g_b2) ! Eq.7, He et al. (2017)
                        gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2) ! Eqn. 3.3 in Fu (2007)
            
                     elseif(snw_shp_lcl(i) == 4) then ! koch snowflake
                     diam_ice = 2. * snw_rds_lcl(i) /0.544
                        if(snw_fs_lcl(i) == 0) then
                           fs_koch = 0.712
                        else
                           fs_koch = snw_fs_lcl(i)               
                        endif
                        fs_hex = 0.788 
                        if(snw_ar_lcl(i) == 0) then
                           AR_tmp = 2.5
                        else
                           AR_tmp = snw_ar_lcl(i)              
                        endif
                        
                        g_ice_Cg_tmp = g_b0 * ((fs_koch/fs_hex)**g_b1) * (diam_ice**g_b2) ! Eq.7, He et al. (2017)
                        gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2) ! Eqn. 3.3 in Fu (2007)
        
                     endif
                     
                     ! 6 wavelength bands for g_ice to be interpolated into 480-bands of SNICAR
                     ! shape-preserving piecewise interpolation into 480-bands
                     if(snw_shp_lcl(i) > 1) then
                        !g_Cg_intp = pchip(g_wvl_center,g_ice_Cg_tmp,wvl) ;
                        !gg_F07_intp = pchip(g_wvl_center,gg_ice_F07_tmp,wvl) ;
                        !data g_wvl(:) /0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00/ ! wavelength (um) division point
                        !g_wvl_center = g_wvl(2:8)/2 + g_wvl(1:7)/2 ; ! center point for wavelength band
                        ! elm wavelength/ /0.3,0.7,1.0,1.2,1.5,5/
                        !wvl_5bd = [0.5 0.85 1.1 1.35 3.25];
                        ! He /0.475 1.055 1.655 2.2 3 3.75 4.5
                     ! linear interpolation to get the Cg and G_f80 for band_idx.
                     if(bnd_idx == 1) then
                        g_Cg_intp = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055-0.475)*(0.5-0.475)+g_ice_Cg_tmp(1);
                        gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055-0.475)*(0.5-0.475)+gg_ice_F07_tmp(1);
                     elseif(bnd_idx == 2) then 
                        g_Cg_intp = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055-0.475)*(0.85-0.475)+g_ice_Cg_tmp(1);
                        gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055-0.475)*(0.85-0.475)+gg_ice_F07_tmp(1);
 
                     elseif(bnd_idx == 3) then 
                        g_Cg_intp = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655-1.055)*(1.1-1.055)+g_ice_Cg_tmp(2);
                        gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655-1.055)*(1.1-1.055)+gg_ice_F07_tmp(2);
                     elseif(bnd_idx == 4) then 
                        g_Cg_intp = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655-1.055)*(1.35-1.055)+g_ice_Cg_tmp(2);
                        gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655-1.055)*(1.35-1.055)+gg_ice_F07_tmp(2);
                     elseif(bnd_idx == 5) then
                        g_Cg_intp = (g_ice_Cg_tmp(6)-g_ice_Cg_tmp(5))/(3.75-3.0)*(3.25-3.0)+g_ice_Cg_tmp(5);
                        gg_F07_intp = (gg_ice_F07_tmp(6)-gg_ice_F07_tmp(5))/(3.75-3.0)*(3.25-3.0)+gg_ice_F07_tmp(5);
                     endif
      
                        g_ice_F07 = gg_F07_intp + (1. - gg_F07_intp) / ss_alb_snw_lcl(i) / 2. ! Eq.2.2 in Fu (2007)
                        g_ice = g_ice_F07 * g_Cg_intp ! Eq.6, He et al. (2017)
                        asm_prm_snw_lcl(i) = g_ice;
                     endif
                     
                     if(asm_prm_snw_lcl(i) > 0.99) then 
                      asm_prm_snw_lcl(i) = 0.99
                     endif                        
                      
                  enddo
                  
                  !  aerosol species 1 optical properties
                  ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                  asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                  ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                  !  aerosol species 2 optical properties
                  ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                  asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                  ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)

                   ! aerosol species 3 optical properties
                   ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)
                   asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
                   ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

                   ! aerosol species 4 optical properties
                   ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)
                   asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
                   ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

                   ! aerosol species 5 optical properties
                   ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)
                   asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
                   ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)

                   ! aerosol species 6 optical properties
                   ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)
                   asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
                   ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)

                   ! aerosol species 7 optical properties
                   ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)
                   asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
                   ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)

                   ! aerosol species 8 optical properties
                   ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)
                   asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
                   ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)


                   ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                   ! 2. optical Depths (tau_snw, tau_aer)
                   ! 3. weighted Mie properties (tau, omega, g)

                   ! Weighted Mie parameters of each layer
                   do i=snl_top,snl_btm,1

                    
                    if (is_BC_internal_mixing) then
                    
                     if(bnd_idx < 4) then
                     ! R_1_omega: BC-induced enhancement in snow single-scattering coalbedo
                     ! R_1_omega from Eq.8b in He et al.(2017,JC) is based on BC Re=0.1um &
                     ! MAC=6.81 m2/g (@550 nm) & BC density=1.7g/cm3.
                     ! To be consistent with SNICAR default (BC MAC=7.5 m2/g @550nm), we
                     ! made adjustments on BC size & density as follows to get MAC=7.5m2/g.
                     ! (1) We use BC Re=0.045um [geometric mean diameter=0.06um (Dentener et al.2006, 
                     ! Yu and Luo,2009) & geometric std=1.5 (Flanner et al.2007;Aoki et al., 2011)]
                     ! (2) We tune BC density from 1.7 to 1.49 g/cm3 (Aoki et al., 2011) to match BC MAC=7.5 m2/g @550 nm. 100
                        C_BC_total = mss_cnc_aer_lcl(i,1) * 1.7/1.49 * 1.0E+09; ! kg/kg to ng/g
                        
                        if (C_BC_total > 0) then
                           R_1_omega_tmp = BC_d0(bnd_idx) * ((C_BC_total + BC_d2(bnd_idx))**BC_d1(bnd_idx)) ! Eq. 8b in He et al.2017,JC
                     ! Adust R_1_omega_tmp due to BC Re from 0.1 to 0.045um based on 
                     ! Eq. 1 & Table S1 in He et al.2018 (GRL)
                           if(bnd_idx == 1) then
                              R_1_omega_tmp = (R_1_omega_tmp / ((0.1/0.05)**(-0.1866)))**((0.1/0.05)**0.1918)  ! visible
                              R_1_omega_tmp = ((0.045/0.05)**(-0.1866)) * (R_1_omega_tmp ** ((0.045/0.05)**(-0.1918))) ! visible
                           else
                              R_1_omega_tmp = (R_1_omega_tmp / ((0.1/0.05)**(-0.0046)))** ((0.1/0.05)**0.5177)  ! NIR
                              R_1_omega_tmp = ((0.045/0.05)**(-0.0046)) * (R_1_omega_tmp ** ((0.045/0.05)**(-0.5177))) ! NIR
                           endif
                     ! new omega for entire BC-snow internal mixture
                           ss_alb_snw_lcl(i) = 1.0 - (1.0 - ss_alb_snw_lcl(i))*R_1_omega_tmp
                                              
                        endif
                     endif
                     
                     ss_alb_aer_lcl(1)     = 0.
                     asm_prm_aer_lcl(1)       = 0.
                     ext_cff_mss_aer_lcl(1)   = 0.
                     !mss_cnc_aer_lcl(i,1) = 0._r8
                    endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! dust internal mixing
                     if (is_dust_internal_mixing) then
                        if (bnd_idx < 4) then
                           C_dust_total = mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) + mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8)
                           C_dust_total = C_dust_total * 1.0E+06 ! kg/kg to ug/g
                           if(C_dust_total > 0) then
                    ! Direct:
                              if (flg_slr_in == 1) then
                                 R_1_omega_tmp = dust_clear_d0(bnd_idx) + dust_clear_d2(bnd_idx)*(C_dust_total**dust_clear_d1(bnd_idx)) ! Eq. 1 in He et al.2019,JAMES                     
                              else
                                 R_1_omega_tmp = dust_cloudy_d0(bnd_idx) + dust_cloudy_d2(bnd_idx)*(C_dust_total**dust_cloudy_d1(bnd_idx)) ! Eq. 1 in He et al.2019,JAMES   
                              endif
                           
                   
                           ! new omega for entire BC-snow internal mixture
                              ss_alb_snw_lcl(i) = 1.0 - (1.0 - ss_alb_snw_lcl(i)) *R_1_omega_tmp
             
                           endif
                        endif
                        do j = 5,8,1
                           ss_alb_aer_lcl(j)     = 0.
                           asm_prm_aer_lcl(j)       = 0.
                           ext_cff_mss_aer_lcl(j)   = 0.
                           !mss_cnc_aer_lcl(i,j) = 0._r8
                        enddo
                    endif
                    

                      L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                      tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                      do j=1,sno_nbr_aer
                         if (is_dust_internal_mixing .and. (j >= 5)) then
                           L_aer(i,j)  = 0.
                         else
                           L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                         endif
                         tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                      enddo

                      tau_sum   = 0.
                      omega_sum = 0.
                      g_sum     = 0.

                      do j=1,sno_nbr_aer
                         tau_sum    = tau_sum + tau_aer(i,j)
                         omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                         g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                      enddo

                      tau(i)    = tau_sum + tau_snw(i)
                      omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                      g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                   enddo ! endWeighted Mie parameters of each layer

                   ! DELTA transformations, if requested
                   if (DELTA == 1) then
                      do i=snl_top,snl_btm,1
                         g_star(i)     = g(i)/(1+g(i))
                         omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                         tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
                      enddo
                   else
                      do i=snl_top,snl_btm,1
                         g_star(i)     = g(i)
                         omega_star(i) = omega(i)
                         tau_star(i)   = tau(i)
                      enddo
                   endif

                   ! Begin radiative transfer solver
                   ! Given input vertical profiles of optical properties, evaluate the
                   ! monochromatic Delta-Eddington adding-doubling solution

                   ! note that trndir, trntdr, trndif, rupdir, rupdif, rdndif
                   ! are variables at the layer interface,
                   ! for snow with layers rangeing from snl_top to snl_btm
                   ! there are snl_top to snl_btm+1 layer interface
                   snl_btm_itf = snl_btm + 1

                   do i = snl_top,snl_btm_itf,1
                      trndir(i) = c0
                      trntdr(i) = c0
                      trndif(i) = c0
                      rupdir(i) = c0
                      rupdif(i) = c0
                      rdndif(i) = c0
                   enddo

                   ! initialize top interface of top layer
                   trndir(snl_top) = c1
                   trntdr(snl_top) = c1
                   trndif(snl_top) = c1
                   rdndif(snl_top) = c0

                  ! begin main level loop
                  ! for layer interfaces except for the very bottom
                  do i = snl_top,snl_btm,1

                     ! initialize all layer apparent optical properties to 0
                     rdir  (i) = c0
                     rdif_a(i) = c0
                     rdif_b(i) = c0
                     tdir  (i) = c0
                     tdif_a(i) = c0
                     tdif_b(i) = c0
                     trnlay(i) = c0

                     ! compute next layer Delta-eddington solution only if total transmission
                     ! of radiation to the interface just above the layer exceeds trmin.

                     if (trntdr(i) > trmin ) then

                        ! calculation over layers with penetrating radiation

                        ! delta-transformed single-scattering properties
                        ! of this layer
                        ts = tau_star(i)
                        ws = omega_star(i)
                        gs = g_star(i)

                       ! Delta-Eddington solution expressions
                        ! n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
                        ! u(w,gg,e)        = c1p5*(c1 - w*gg)/e
                        ! el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
                        lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))  !lm = el(ws,gs)
                        ue = c1p5*(c1 - ws*gs)/lm           !ue = u(ws,gs,lm)
                        extins = max(exp_min, exp(-lm*ts))
                        ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins) !ne = n(ue,extins)

                        ! first calculation of rdif, tdif using Delta-Eddington formulas
                        ! rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
                        rdif_a(i) = (ue**2-c1)*(c1/extins - extins)/ne
                        tdif_a(i) = c4*ue/ne

                        ! evaluate rdir,tdir for direct beam
                        trnlay(i) = max(exp_min, exp(-ts/mu_not))

                        ! Delta-Eddington solution expressions
                        ! alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
                        ! agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
                        ! alp = alpha(ws,mu_not,gs,lm)
                        ! gam = agamm(ws,mu_not,gs,lm)
                        alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
                        gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
                        apg = alp + gam
                        amg = alp - gam

                        rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
                        tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)

                        ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
                        ! since Delta-Eddington rdif formula is not well-behaved (it is usually
                        ! biased low and can even be negative); use ngmax angles and gaussian
                        ! integration for most accuracy:
                        R1 = rdif_a(i) ! use R1 as temporary
                        T1 = tdif_a(i) ! use T1 as temporary
                        swt = c0
                        smr = c0
                        smt = c0
                        do ng=1,ngmax
                           mu  = difgauspt(ng)
                           gwt = difgauswt(ng)
                           swt = swt + mu*gwt
                           trn = max(exp_min, exp(-ts/mu))
                           ! alp = alpha(ws,mu,gs,lm)
                           ! gam = agamm(ws,mu,gs,lm)
                           alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
                           gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
                           apg = alp + gam
                           amg = alp - gam
                           rdr = apg*R1 + amg*T1*trn - amg
                           tdr = apg*T1 + amg*R1*trn - apg*trn + trn
                           smr = smr + mu*rdr*gwt
                           smt = smt + mu*tdr*gwt
                        enddo      ! ng
                        rdif_a(i) = smr/swt
                        tdif_a(i) = smt/swt

                        ! homogeneous layer
                        rdif_b(i) = rdif_a(i)
                        tdif_b(i) = tdif_a(i)

                      endif ! trntdr(k) > trmin

                      ! Calculate the solar beam transmission, total transmission, and
                      ! reflectivity for diffuse radiation from below at interface i,
                      ! the top of the current layer k:
                      !
                      !              layers       interface
                      !
                      !       ---------------------  i-1
                      !                i-1
                      !       ---------------------  i
                      !                 i
                      !       ---------------------

                      trndir(i+1) = trndir(i)*trnlay(i)
                      refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))
                      tdrrdir     = trndir(i)*rdir(i)
                      tdndif      = trntdr(i) - trndir(i)
                      trntdr(i+1) = trndir(i)*tdir(i) + &
                           (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
                      rdndif(i+1) = rdif_b(i) + &
                           (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
                      trndif(i+1) = trndif(i)*refkm1*tdif_a(i)

                  enddo       ! i    end main level loop


                  ! compute reflectivity to direct and diffuse radiation for layers
                  ! below by adding succesive layers starting from the underlying
                  ! ground and working upwards:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------  i+1
                  !                i+1
                  !       ---------------------

                  ! set the underlying ground albedo == albedo of near-IR
                  ! unless bnd_idx == 1, for visible
                  rupdir(snl_btm_itf) = albsfc(c_idx,2)
                  rupdif(snl_btm_itf) = albsfc(c_idx,2)
                  if (bnd_idx == 1) then
                      rupdir(snl_btm_itf) = albsfc(c_idx,1)
                      rupdif(snl_btm_itf) = albsfc(c_idx,1)
                  endif

                  do i=snl_btm,snl_top,-1
                     ! interface scattering
                     refkp1        = c1/( c1 - rdif_b(i)*rupdif(i+1))
                     ! dir from top layer plus exp tran ref from lower layer, interface
                     ! scattered and tran thru top layer from below, plus diff tran ref
                     ! from lower layer with interface scattering tran thru top from below
                     rupdir(i) = rdir(i) &
                          + (        trnlay(i)  *rupdir(i+1) &
                          +  (tdir(i)-trnlay(i))*rupdif(i+1))*refkp1*tdif_b(i)
                     ! dif from top layer from above, plus dif tran upwards reflected and
                     ! interface scattered which tran top from below
                     rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
                  enddo       ! i

                  ! net flux (down-up) at each layer interface from the
                  ! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
                  ! the interface reflectivities and transmissivities required
                  ! to evaluate interface fluxes are returned from solution_dEdd;
                  ! now compute up and down fluxes for each interface, using the
                  ! combined layer properties at each interface:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------

                  do i = snl_top, snl_btm_itf
                     ! interface scattering
                     refk          = c1/(c1 - rdndif(i)*rupdif(i))
                     ! dir tran ref from below times interface scattering, plus diff
                     ! tran and ref from below times interface scattering
                     ! fdirup(i) = (trndir(i)*rupdir(i) + &
                     !                 (trntdr(i)-trndir(i))  &
                     !                 *rupdif(i))*refk
                     ! dir tran plus total diff trans times interface scattering plus
                     ! dir tran with up dir ref and down dif ref times interface scattering
                     ! fdirdn(i) = trndir(i) + (trntdr(i) &
                     !               - trndir(i) + trndir(i)  &
                     !               *rupdir(i)*rdndif(i))*refk
                     ! diffuse tran ref from below times interface scattering
                     ! fdifup(i) = trndif(i)*rupdif(i)*refk
                     ! diffuse tran times interface scattering
                     ! fdifdn(i) = trndif(i)*refk

                     ! netflux, down - up
                     ! dfdir = fdirdn - fdirup
                     dfdir(i) = trndir(i) &
                                 + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
                                 -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
                     if (dfdir(i) < puny) dfdir(i) = c0
                     ! dfdif = fdifdn - fdifup
                     dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
                     if (dfdif(i) < puny) dfdif(i) = c0
                  enddo       ! k

                  ! SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
                  ! direct incident
                  if (flg_slr_in == 1) then
                    albedo = rupdir(snl_top)
                    dftmp  = dfdir
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = (trndir(snl_top)*rupdir(snl_top) + &
                                (trntdr(snl_top)-trndir(snl_top))  &
                                 *rupdif(snl_top))*refk
                  !diffuse incident
                  else
                    albedo = rupdif(snl_top)
                    dftmp  = dfdif
                    refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
                    F_sfc_pls = trndif(snl_top)*rupdif(snl_top)*refk
                  endif

                  ! Absorbed flux in each layer
                  do i=snl_top,snl_btm,1
                    F_abs(i) = dftmp(i)-dftmp(i+1)
                    flx_abs_lcl(i,bnd_idx) = F_abs(i)

                    ! ERROR check: negative absorption
                    if (flx_abs_lcl(i,bnd_idx) < -0.00001) then
                      write (*,"(a,e13.6,a,i6,a,i6)") "SNICAR ERROR: negative absoption : ", flx_abs_lcl(i,bnd_idx), &
                           " at timestep: ", nstep, " at column: ", c_idx
                      write(*,*) "SNICAR_AD STATS: snw_rds(0)= ", snw_rds(c_idx,0)
                      write(*,*) "SNICAR_AD STATS: L_snw(0)= ", L_snw(0)
                      write(*,*) "SNICAR_AD STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                      write(*,*) "SNICAR_AD STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
                      write(*,*) "SNICAR_AD STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
                      write(*,*) "SNICAR_AD STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
                      write(*,*) "SNICAR_AD STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
                      write(*,*) "SNICAR_AD STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
                      write(*,*) "SNICAR_AD STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
                      write(*,*) "SNICAR_AD STATS: ss_alb_snw_lcl (0)= ", ss_alb_snw_lcl(i) ! add by Dalei check
                      write(*,*) "SNICAR_AD STATS: asm_prm_snw_lcl (0)= ", asm_prm_snw_lcl(i) ! add by Dalei check
                      write(*,*) "SNICAR_AD STATS: ext_cff_mss_snw_lcl (0)= ", ext_cff_mss_snw_lcl(i) ! add by Dalei check
                      write(*,*) "SNICAR_AD STATS: g_star (0)= ", g_star(i) ! add by Dalei check
                      write(*,*) "SNICAR_AD STATS: omega_star (0)= ", omega_star(i) ! add by Dalei check
                      write(*,*) "SNICAR_AD STATS: tau_star (0)= ", tau_star(i) ! add by Dalei check
                      write(*,*) "g_Cg_intp ", g_Cg_intp ! add by Dalei check
                      write(*,*) "gg_F07_intp ", gg_F07_intp ! add by Dalei check
                      write(*,*) "g_ice_F07 ", g_ice_F07 ! add by Dalei check
                      write(*,*) "bnd_idx ", bnd_idx ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(1) ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(2) ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(3) ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(4) ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(5) ! add by Dalei check
                      write(*,*) "gg_ice_F07_tmp ", gg_ice_F07_tmp(6) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(1) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(2) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(3) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(4) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(5) ! add by Dalei check
                      write(*,*) "g_ice_Cg_tmp ", g_ice_Cg_tmp(6) ! add by Dalei check
                      
                     call land_error_message("SNICAR_AD_RT in snicar_mod: Negative absorption!", severity=FATAL)
                    endif
                  enddo

                  ! absobed flux by the underlying ground
                  F_btm_net = dftmp(snl_btm_itf)

                  ! note here, snl_btm_itf = 1 by snow column set up in CLM
                  flx_abs_lcl(1,bnd_idx) = F_btm_net

                 if (flg_nosnl == 1) then
                    ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                    !flx_abs_lcl(:,bnd_idx) = 0._r8
                    !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

                    ! changed on 20070408:
                    ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                    ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                    ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                    flx_abs_lcl(0,bnd_idx) = F_abs(0)
                    flx_abs_lcl(1,bnd_idx) = F_btm_net
                 endif

                 !Underflow check (we've already tripped the error condition above)
                 do i=snl_top,1,1
                    if (flx_abs_lcl(i,bnd_idx) < 0.) then
                       flx_abs_lcl(i,bnd_idx) = 0.
                    endif
                 enddo

                 F_abs_sum = 0.
                 do i=snl_top,snl_btm,1
                    F_abs_sum = F_abs_sum + F_abs(i)
                 enddo

                !enddo !enddo while (flg_dover > 0)

                ! Energy conservation check:
                ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
                energy_sum = (mu_not*pi*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
                if (abs(energy_sum) > 0.00001) then
                   write (*,"(a,e13.6,a,i6,a,i6)") "SNICAR ERROR: Energy conservation error of : ", energy_sum, &
                        " at timestep: ", nstep, " at column: ", c_idx
                   write(*,*) "F_abs_sum: ",F_abs_sum
                   write(*,*) "F_btm_net: ",F_btm_net
                   write(*,*) "F_sfc_pls: ",F_sfc_pls
                   write(*,*) "mu_not*pi*flx_slrd_lcl(bnd_idx): ", mu_not*pi*flx_slrd_lcl(bnd_idx)
                   write(*,*) "flx_slri_lcl(bnd_idx)", flx_slri_lcl(bnd_idx)
                   write(*,*) "bnd_idx", bnd_idx
                   write(*,*) "F_abs", F_abs
                   write(*,*) "albedo", albedo
                  call land_error_message("SNICAR_AD_RT in snicar_mod: Energy conservation error!", severity=FATAL)
                endif

                albout_lcl(bnd_idx) = albedo
                ! Check that albedo is less than 1
                if (albout_lcl(bnd_idx) > 1.0) then
                   write (*,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx, " NSTEP= ",nstep
                   write (*,*) "SNICAR STATS: bnd_idx= ",bnd_idx
                   write (*,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), &
                        " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
                   write (*,*) "SNICAR STATS: landtype= ", sfctype
                   write (*,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                   write (*,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in

                   write (*,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
                   write (*,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
                   write (*,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
                   write (*,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
                   write (*,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)

                   write (*,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
                   write (*,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
                   write (*,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
                   write (*,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
                   write (*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

                   write (*,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
                   write (*,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
                   write (*,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
                   write (*,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
                   write (*,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)

                  call land_error_message("SNICAR_AD_RT in snicar_mod: Albedo out of bounds!", severity=FATAL)
                endif

             enddo   ! loop over wvl bands


             ! Weight output NIR albedo appropriately
             albout(c_idx,1) = albout_lcl(1)
             flx_sum         = 0.
             do bnd_idx= nir_bnd_bgn,nir_bnd_end
                flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
             enddo
             albout(c_idx,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

             ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
             flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
             do i=snl_top,1,1
                flx_sum = 0.
                do bnd_idx= nir_bnd_bgn,nir_bnd_end
                   flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
                enddo
                flx_abs(c_idx,i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
             enddo

             ! near-IR direct albedo/absorption adjustment for high solar zenith angles
             ! solar zenith angle parameterization
             ! calculate the scaling factor for NIR direct albedo if SZA>75 degree
             if ((mu_not < mu_75) .and. (flg_slr_in == 1)) then
                sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * mu_not**2
                sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * mu_not**2
                sza_factor = sza_c1 * (log10(snw_rds_lcl(snl_top) * c1) - c6) + sza_c0
                flx_sza_adjust  = albout(c_idx,2) * (sza_factor-c1) * sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
                albout(c_idx,2) = albout(c_idx,2) * sza_factor
                flx_abs(c_idx,snl_top,2) = flx_abs(c_idx,snl_top,2) - flx_sza_adjust
             endif

             ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
          elseif ( (coszen(c_idx) > 0.) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0.) ) then
             albout(c_idx,1) = albsfc(c_idx,1)
             albout(c_idx,2) = albsfc(c_idx,2)

             ! There is either zero snow, or no sun
          else
             albout(c_idx,1) = 0.
             albout(c_idx,2) = 0.
          endif    ! if column has snow and coszen > 0

       enddo    ! loop over all columns

   !   end associate

   end subroutine SNICAR_AD_RT



     !-----------------------------------------------------------------------
  subroutine SNICAR_RT_HE(nlevsno, flg_snw_ice, &
   coszen, flg_slr_in, h2osno_liq, h2osno_ice, snw_rds, snw_shp_input, &
   mss_cnc_aer_in, albsfc, albout, flx_abs)


! !DESCRIPTION:
! Determine reflectance of, and vertically-resolved solar absorption in, 
! snow with impurities.
!
! Original references on physical models of snow reflectance include: 
! Wiscombe and Warren [1980] and Warren and Wiscombe [1980],
! Journal of Atmospheric Sciences, 37,
!
! The multi-layer solution for multiple-scattering used here is from:
! Toon et al. [1989], Rapid calculation of radiative heating rates 
! and photodissociation rates in inhomogeneous multiple scattering atmospheres, 
! J. Geophys. Res., 94, D13, 16287-16301
!
! The implementation of the SNICAR model in CLM/CSIM is described in:
! Flanner, M., C. Zender, J. Randerson, and P. Rasch [2007], 
! Present-day climate forcing and response from black carbon in snow,
! J. Geophys. Res., 112, D11202, doi: 10.1029/2006JD008003
!
! Updated radiative transfer solver:
!
! The multi-layer solution for multiple-scattering used here is from:
! Briegleb, P. and Light, B.: A Delta-Eddington mutiple scattering
! parameterization for solar radiation in the sea ice component of the
! community climate system model, 2007.
!
! The implementation of the SNICAR-AD model in CLM is described in:
! Dang et al.2019, Inter-comparison and improvement of 2-stream shortwave
! radiative transfer models for unified treatment of cryospheric surfaces
! in ESMs; and Flanner et al. 2021, SNICAR-ADv3: a community tool for modeling 
! spectral snow albedo
!
! !USES:
! use clm_varpar       , only : nlevsno, numrad
! use clm_time_manager , only : get_nstep
! use shr_const_mod    , only : SHR_CONST_PI
!
! ENRICO ZORZETTO 2023:
! modifications for implementation in lm4p2:
! added as input number of snow layers (pos. integer nlevsno)
! snicar is called for a single column (bounds, num_nourbanc, filter_nourbanc -> not used )
!
! !ARGUMENTS:
integer, INTENT(IN) :: nlevsno ! number of snow layers
integer           , intent(in)  :: flg_snw_ice                                        ! flag: =1 when called from CLM, =2 when called from CSIM
! type (bounds_type), intent(in)  :: bounds                                    
! integer           , intent(in)  :: num_nourbanc                                       ! number of columns in non-urban filter
! integer           , intent(in)  :: filter_nourbanc(:)                                 ! column filter for non-urban points
real          , intent(in)  :: coszen         ( 1: )                    ! cosine of solar zenith angle for next time step (col) [unitless]
integer           , intent(in)  :: flg_slr_in                                         ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux
real          , intent(in)  :: h2osno_liq     ( 1: , -nlevsno+1: )      ! liquid water content (col,lyr) [kg/m2]
real          , intent(in)  :: h2osno_ice     ( 1: , -nlevsno+1: )      ! ice content (col,lyr) [kg/m2]
integer           , intent(in)  :: snw_rds        ( 1: , -nlevsno+1: )      ! snow effective radius (col,lyr) [microns, m^-6]
integer           , intent(in)  :: snw_shp_input        ( 1: , -nlevsno+1: )      ! snow shape (col, lyr) integer code to be converted to string ! EZSNOW
real          , intent(in)  :: mss_cnc_aer_in ( 1: , -nlevsno+1: , 1: ) ! mass concentration of all aerosol species (col,lyr,aer) [kg/kg]
real          , intent(in)  :: albsfc         ( 1: , 1: )               ! albedo of surface underlying snow (col,bnd) [frc]
real          , intent(out) :: albout         ( 1: , 1: )               ! snow albedo, averaged into 2 bands (=0 if no sun or no snow) (col,bnd) [frc]
real          , intent(out) :: flx_abs        ( 1: , -nlevsno+1: , 1: ) ! absorbed flux in each layer per unit flux incident (col, lyr, bnd)
! type(waterdiagnosticbulk_type) , intent(in)  :: waterdiagnosticbulk_inst
!
! !LOCAL VARIABLES:
!
! variables for snow radiative transfer calculations
! integer :: nir_bnd_bgn  ! first band index in near-IR spectrum [idx]
! integer :: nir_bnd_end  ! ending near-IR band index [idx]

! Local variables representing single-column values of arrays:
! real :: h2osno_total   ( 1: )                    ! total snow content (col) [kg/m2] ! EZSNOW moved here, compute from sum(liq + ice)
integer :: snl_lcl                            ! negative number of snow layers [nbr]
integer :: snw_rds_lcl(-nlevsno+1:0)          ! snow effective radius [m^-6]
real :: flx_slrd_lcl(1:numrad_snw)         ! direct beam incident irradiance [W/m2] (set to 1)
real :: flx_slri_lcl(1:numrad_snw)         ! diffuse incident irradiance [W/m2] (set to 1)
real :: mss_cnc_aer_lcl(-nlevsno+1:0,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
real :: h2osno_lcl                         ! total column snow mass [kg/m2]
real :: h2osno_liq_lcl(-nlevsno+1:0)       ! liquid water mass [kg/m2]
real :: h2osno_ice_lcl(-nlevsno+1:0)       ! ice mass [kg/m2]
real :: albsfc_lcl(1:numrad_snw)           ! albedo of underlying surface [frc]
real :: ss_alb_snw_lcl(-nlevsno+1:0)       ! single-scatter albedo of ice grains (lyr) [frc]
real :: asm_prm_snw_lcl(-nlevsno+1:0)      ! asymmetry parameter of ice grains (lyr) [frc]
real :: ext_cff_mss_snw_lcl(-nlevsno+1:0)  ! mass extinction coefficient of ice grains (lyr) [m2/kg]
real :: ss_alb_aer_lcl(sno_nbr_aer)        ! single-scatter albedo of aerosol species (aer_nbr) [frc] 
real :: asm_prm_aer_lcl(sno_nbr_aer)       ! asymmetry parameter of aerosol species (aer_nbr) [frc]
real :: ext_cff_mss_aer_lcl(sno_nbr_aer)   ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]

! Other local variables
integer :: APRX_TYP                           ! two-stream approximation type
                             ! (1=Eddington, 2=Quadrature, 3=Hemispheric Mean) [nbr]
integer :: DELTA                              ! flag to use Delta approximation (Joseph, 1976)
                             ! (1= use, 0= don't use)
real :: flx_wgt(1:numrad_snw)              ! weights applied to spectral bands,
                             ! specific to direct and diffuse cases (bnd) [frc] 
integer :: flg_nosnl                          ! flag: =1 if there is snow, but zero snow layers,
                             ! =0 if at least 1 snow layer [flg]   
integer :: trip                               ! flag: =1 to redo RT calculation if result is unrealistic
integer :: flg_dover                          ! defines conditions for RT redo (explained below)
real :: albedo                             ! temporary snow albedo [frc]
real :: flx_sum                            ! temporary summation variable for NIR weighting
real :: albout_lcl(numrad_snw)             ! snow albedo by band [frc]
real :: flx_abs_lcl(-nlevsno+1:1,numrad_snw)! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
real :: L_snw(-nlevsno+1:0)                ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
real :: tau_snw(-nlevsno+1:0)              ! snow optical depth (lyr) [unitless]
real :: L_aer(-nlevsno+1:0,sno_nbr_aer)    ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2] 
real :: tau_aer(-nlevsno+1:0,sno_nbr_aer)  ! aerosol optical depth (lyr,nbr_aer) [unitless]
real :: tau_sum                            ! cumulative (snow+aerosol) optical depth [unitless]
real :: tau_clm(-nlevsno+1:0)              ! column optical depth from layer bottom to snowpack top (lyr) [unitless] 
real :: omega_sum                          ! temporary summation of single-scatter albedo of all aerosols [frc]
real :: g_sum                              ! temporary summation of asymmetry parameter of all aerosols [frc]
real :: tau(-nlevsno+1:0)                  ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
real :: omega(-nlevsno+1:0)                ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
real :: g(-nlevsno+1:0)                    ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
real :: tau_star(-nlevsno+1:0)             ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                             ! (lyr) [unitless]
real :: omega_star(-nlevsno+1:0)           ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
real :: g_star(-nlevsno+1:0)               ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                             ! (lyr) [frc]
integer :: nstep                              ! current timestep [nbr] (debugging only)
integer :: g_idx, c_idx, l_idx                ! gridcell, column, and landunit indices [idx]
integer :: bnd_idx                            ! spectral band index (1 <= bnd_idx <= snicar_numrad_snw) [idx]
integer :: rds_idx                            ! snow effective radius index for retrieving
                             ! Mie parameters from lookup table [idx]
integer :: snl_btm                            ! index of bottom snow layer (0) [idx]
integer :: snl_top                            ! index of top snow layer (-4 to 0) [idx]
integer :: fc                                 ! column filter index
integer :: i                                  ! layer index [idx]
integer :: j                                  ! aerosol number index [idx]
integer :: n                                  ! tridiagonal matrix index [idx]
integer :: m                                  ! secondary layer index [idx]   
real :: F_direct(-nlevsno+1:0)             ! direct-beam radiation at bottom of layer interface (lyr) [W/m^2]
real :: F_net(-nlevsno+1:0)                ! net radiative flux at bottom of layer interface (lyr) [W/m^2]
real :: F_abs(-nlevsno+1:0)                ! net absorbed radiative energy (lyr) [W/m^2]
real :: F_abs_sum                          ! total absorbed energy in column [W/m^2]
real :: F_sfc_pls                          ! upward radiative flux at snowpack top [W/m^2]
real :: F_btm_net                          ! net flux at bottom of snowpack [W/m^2]                    
real :: F_sfc_net                          ! net flux at top of snowpack [W/m^2]
real :: energy_sum                         ! sum of all energy terms; should be 0.0 [W/m^2]
real :: F_direct_btm                       ! direct-beam radiation at bottom of snowpack [W/m^2]
real :: mu_not                             ! cosine of solar zenith angle (used locally) [frc]
integer :: err_idx                            ! counter for number of times through error loop [nbr]
real :: lat_coord                          ! gridcell latitude (debugging only)
real :: lon_coord                          ! gridcell longitude (debugging only)
integer :: sfctype                            ! underlying surface type (debugging only)
real :: pi                                 ! 3.1415...

!-----------------------------------------------------------------------
! variables used for Toon et al. 1989 2-stream solver (Flanner et al. 2007):
! intermediate variables for radiative transfer approximation:
real :: gamma1(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: gamma2(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: gamma3(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: gamma4(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: lambda(-nlevsno+1:0)               ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: GAMMA(-nlevsno+1:0)                ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: mu_one                             ! two-stream coefficient from Toon et al. (lyr) [unitless]
real :: e1(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
real :: e2(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
real :: e3(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
real :: e4(-nlevsno+1:0)                   ! tri-diag intermediate variable from Toon et al. (lyr) 
real :: C_pls_btm(-nlevsno+1:0)            ! intermediate variable: upward flux at bottom interface (lyr) [W/m2]
real :: C_mns_btm(-nlevsno+1:0)            ! intermediate variable: downward flux at bottom interface (lyr) [W/m2]
real :: C_pls_top(-nlevsno+1:0)            ! intermediate variable: upward flux at top interface (lyr) [W/m2]
real :: C_mns_top(-nlevsno+1:0)            ! intermediate variable: downward flux at top interface (lyr) [W/m2]
real :: A(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: B(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: D(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: E(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: AS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: DS(-2*nlevsno+1:0)                 ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: X(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)
real :: Y(-2*nlevsno+1:0)                  ! tri-diag intermediate variable from Toon et al. (2*lyr)

!-----------------------------------------------------------------------
! variables used for Adding-doubling 2-stream solver based on SNICAR-ADv3 version 
! (Dang et al. 2019; Flanner et al. 2021)
real :: trndir(-nlevsno+1:1)               ! solar beam down transmission from top
real :: trntdr(-nlevsno+1:1)               ! total transmission to direct beam for layers above
real :: trndif(-nlevsno+1:1)               ! diffuse transmission to diffuse beam for layers above
real :: rupdir(-nlevsno+1:1)               ! reflectivity to direct radiation for layers below
real :: rupdif(-nlevsno+1:1)               ! reflectivity to diffuse radiation for layers below
real :: rdndif(-nlevsno+1:1)               ! reflectivity to diffuse radiation for layers above
real :: dfdir(-nlevsno+1:1)                ! down-up flux at interface due to direct beam at top surface
real :: dfdif(-nlevsno+1:1)                ! down-up flux at interface due to diffuse beam at top surface
real :: dftmp(-nlevsno+1:1)                ! temporary variable for down-up flux at interface
real :: rdir(-nlevsno+1:0)                 ! layer reflectivity to direct radiation
real :: rdif_a(-nlevsno+1:0)               ! layer reflectivity to diffuse radiation from above
real :: rdif_b(-nlevsno+1:0)               ! layer reflectivity to diffuse radiation from below
real :: tdir(-nlevsno+1:0)                 ! layer transmission to direct radiation (solar beam + diffuse)
real :: tdif_a(-nlevsno+1:0)               ! layer transmission to diffuse radiation from above
real :: tdif_b(-nlevsno+1:0)               ! layer transmission to diffuse radiation from below
real :: trnlay(-nlevsno+1:0)               ! solar beam transm for layer (direct beam only)
real :: ts                                 ! layer delta-scaled extinction optical depth
real :: ws                                 ! layer delta-scaled single scattering albedo
real :: gs                                 ! layer delta-scaled asymmetry parameter
real :: extins                             ! extinction
real :: alp                                ! temporary for alpha
real :: gam                                ! temporary for agamm
real :: amg                                ! alp - gam
real :: apg                                ! alp + gam
real :: ue                                 ! temporary for u
real :: refk                               ! interface multiple scattering
real :: refkp1                             ! interface multiple scattering for k+1
real :: refkm1                             ! interface multiple scattering for k-1
real :: tdrrdir                            ! direct tran times layer direct ref
real :: tdndif                             ! total down diffuse = tot tran - direct tran
real :: taus                               ! scaled extinction optical depth
real :: omgs                               ! scaled single particle scattering albedo
real :: asys                               ! scaled asymmetry parameter
real :: lm                                 ! temporary for el
real :: mu                                 ! cosine solar zenith for either snow or water
real :: ne                                 ! temporary for n
real :: R1                                 ! perpendicular polarization reflection amplitude
real :: R2                                 ! parallel polarization reflection amplitude
real :: T1                                 ! perpendicular polarization transmission amplitude
real :: T2                                 ! parallel polarization transmission amplitude
real :: Rf_dir_a                           ! fresnel reflection to direct radiation
real :: Tf_dir_a                           ! fresnel transmission to direct radiation
real :: Rf_dif_a                           ! fresnel reflection to diff radiation from above
real :: Rf_dif_b                           ! fresnel reflection to diff radiation from below
real :: Tf_dif_a                           ! fresnel transmission to diff radiation from above
real :: Tf_dif_b                           ! fresnel transmission to diff radiation from below
real :: gwt                                ! gaussian weight
real :: swt                                ! sum of weights
real :: trn                                ! layer transmission
real :: rdr                                ! rdir for gaussian integration
real :: tdr                                ! tdir for gaussian integration
real :: smr                                ! accumulator for rdif gaussian integration
real :: smt                                ! accumulator for tdif gaussian integration
real :: exp_min                            ! minimum exponential value
real, allocatable :: difgauspt(:)         ! Gaussian integration angle
real, allocatable :: difgauswt(:)         ! Gaussian integration coefficients/weights
integer :: ng                                 ! gaussian integration index
integer :: ngmax = 8                          ! max gaussian integration index
integer :: snl_btm_itf                        ! index of bottom snow layer interfaces (1) [idx]
! constants used in algorithm
real :: c0     = 0.0
real :: c1     = 1.0
real :: c3     = 3.0
real :: c4     = 4.0
real :: c6     = 6.0
real :: cp01   = 0.01
real :: cp5    = 0.5
real :: cp75   = 0.75
real :: c1p5   = 1.5
real :: trmin  = 0.001
real :: argmax = 10.0                   ! maximum argument of exponential
! cconstant and coefficients used for SZA parameterization
real :: sza_a0 =  0.085730
real :: sza_a1 = -0.630883
real :: sza_a2 =  1.303723
real :: sza_b0 =  1.467291
real :: sza_b1 = -3.338043
real :: sza_b2 =  6.807489
real :: puny   =  1.0e-11
real :: mu_75  =  0.2588                ! cosine of 75 degree
real :: sza_c1                             ! coefficient, SZA parameteirzation
real :: sza_c0                             ! coefficient, SZA parameterization
real :: sza_factor                         ! factor used to adjust NIR direct albedo
real :: flx_sza_adjust                     ! direct NIR flux adjustment from sza_factor
real :: mu0                                ! incident solar zenith angle

!-----------------------------------------------------------------------
! variables used for nonspherical snow grain treatment (He et al. 2017 J of Climate):
character(len=15) :: sno_shp(-nlevsno+1:0)    ! Snow shape type: sphere, spheroid, hexagonal plate, koch snowflake
                             ! currently only assuming same shapes for all snow layers
real :: sno_fs(-nlevsno+1:0)              ! Snow shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                             ! only activated when snicar_snw_shape is nonspherical
                             ! 0=use recommended default value (He et al. 2017);
                             ! others(0<sno_fs<1)= user-specified value
real :: sno_AR(-nlevsno+1:0)              ! Snow grain aspect ratio: ratio of grain width to length
                             ! only activated when snicar_snw_shape is nonspherical
                             ! 0=use recommended default value (He et al. 2017);
                             ! others(0.1<fs<20)= use user-specified value
! Constants and parameters for aspherical ice particles    
! asymmetry factor parameterization coefficients (6 bands) from Table 3 & Eqs. 6-7 in He et al. (2017)
real :: g_wvl(1:8)                        ! wavelength (um) division point
real :: g_wvl_ct(1:7)                     ! center point for wavelength band (um)
real :: g_b0(1:7)
real :: g_b1(1:7)
real :: g_b2(1:7)
! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 JAS
real :: g_F07_c2(1:7)
real :: g_F07_c1(1:7)
real :: g_F07_c0(1:7)
real :: g_F07_p2(1:7)
real :: g_F07_p1(1:7)
real :: g_F07_p0(1:7)
! other temporary variables
real, allocatable :: wvl_ct(:)            ! band center wavelength (um) for 5 or 480-band case
real :: diam_ice                          ! effective snow grain diameter (SSA-equivalent) unit: microns
real :: fs_sphd                           ! shape factor for spheroid snow
real :: fs_hex                            ! shape factor for reference hexagonal snow
real :: fs_hex0                           ! shape factor for hexagonal plate
real :: fs_koch                           ! shape factor for Koch snowflake
real :: AR_tmp                            ! aspect ratio temporary
real :: g_ice_Cg_tmp(1:7)                 ! temporary asymmetry factor correction coeff
real :: gg_ice_F07_tmp(1:7)               ! temporary asymmetry factor related to geometric reflection & refraction
real :: g_Cg_intp                         ! interpolated asymmetry factor correction coeff to target bands 
real :: gg_F07_intp                       ! interpolated asymmetry factor related to geometric reflection & refraction
real :: g_ice_F07                         ! asymmetry factor for Fu 2007 parameterization value
integer  :: igb                               ! loop index

!-----------------------------------------------------------------------
! variables used for BC-snow internal mixing (He et al. 2017 J of Climate):
real :: enh_omg_bcint                     ! BC-induced enhancement in snow single-scattering co-albedo (1-omega)
real :: enh_omg_bcint_tmp(1:16)           ! temporary BC-induced enhancement in snow 1-omega
real :: enh_omg_bcint_tmp2(1:16)          ! temporary BC-induced enhancement in snow 1-omega
real :: bcint_wvl(1:17)                   ! Parameterization band (0.2-1.2um) for BC-induced enhancement in snow 1-omega
real :: bcint_wvl_ct(1:16)                ! Parameterization band center wavelength (um)
real :: bcint_d0(1:16)                    ! Parameterization coefficients at each band center wavelength
real :: bcint_d1(1:16)                    ! Parameterization coefficients at each band center wavelength
real :: bcint_d2(1:16)                    ! Parameterization coefficients at each band center wavelength
real :: den_bc = 1.49                  ! target BC particle density (g/cm3) used in BC MAC adjustment
real :: Re_bc = 0.045                     ! target BC effective radius (um) used in BC MAC adjustment
real :: bcint_m(1:3)                      ! Parameterization coefficients for BC size adjustment in BC-snow int mix
real :: bcint_n(1:3)                      ! Parameterization coefficients for BC size adjustment in BC-snow int mix
real :: bcint_m_tmp                       ! temporary of bcint_m
real :: bcint_n_tmp                       ! temporary of bcint_n
real :: bcint_dd                          ! intermediate parameter
real :: bcint_dd2                         ! intermediate parameter
real :: bcint_f                           ! intermediate parameter
real :: enh_omg_bcint_intp                ! BC-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
real :: enh_omg_bcint_intp2               ! BC-induced enhancement in snow 1-omega interpolated to CLM wavelength
real :: wvl_doint                         ! wavelength doing BC-snow int mixing (<=1.2um)
integer  :: ibb                               ! loop index

!-----------------------------------------------------------------------
! variables used for dust-snow internal mixing (He et al. 2019 JAMES):
real :: enh_omg_dstint                    ! dust-induced enhancement in snow single-scattering co-albedo (1-omega)
real :: enh_omg_dstint_tmp(1:6)           ! temporary dust-induced enhancement in snow 1-omega
real :: enh_omg_dstint_tmp2(1:6)          ! temporary dust-induced enhancement in snow 1-omega
real :: dstint_wvl(1:7)                   ! Parameterization band (0.2-1.2um) for dust-induced enhancement in snow 1-omega
real :: dstint_wvl_ct(1:6)                ! Parameterization band center wavelength (um)
real :: dstint_a1(1:6)                    ! Parameterization coefficients at each band center wavelength
real :: dstint_a2(1:6)                    ! Parameterization coefficients at each band center wavelength
real :: dstint_a3(1:6)                    ! Parameterization coefficients at each band center wavelength
real :: enh_omg_dstint_intp               ! dust-induced enhancement in snow 1-omega (logscale) interpolated to CLM wavelength
real :: enh_omg_dstint_intp2              ! dust-induced enhancement in snow 1-omega interpolated to CLM wavelength
real :: tot_dst_snw_conc                  ! total dust content in snow across all size bins (ppm=ug/g)
integer  :: idb                               ! loop index

real :: h2osno_total(num_nourbanc)
real :: snl(num_nourbanc)

!-----------------------------------------------------------------------

! Enforce expected array sizes
! SHR_ASSERT_ALL_FL((ubound(coszen)         == (/bounds%endc/)),                 sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(h2osno_liq)     == (/bounds%endc, 0/)),              sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(h2osno_ice)     == (/bounds%endc, 0/)),              sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(h2osno_total)   == (/bounds%endc/)),                 sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(snw_rds)        == (/bounds%endc, 0/)),              sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(mss_cnc_aer_in) == (/bounds%endc, 0, sno_nbr_aer/)), sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(albsfc)         == (/bounds%endc, numrad/)),         sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(albout)         == (/bounds%endc, numrad/)),         sourcefile, __LINE__)
! SHR_ASSERT_ALL_FL((ubound(flx_abs)        == (/bounds%endc, 1, numrad/)),      sourcefile, __LINE__)

! associate(& 
! snl         =>   col%snl                           , & ! Input:  [integer (:)]  negative number of snow layers (col) [nbr]
! frac_sno    =>   waterdiagnosticbulk_inst%frac_sno_eff_col    & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
! )

! snl = - nlevsno ! EZSNOW
! frac_sno = 1.0 ! EZSNOW

h2osno_total(1) = sum(h2osno_ice + h2osno_liq) ! EZSNOW
snl(1) = -nlevsno ! EZSNOW

! initialize parameter
! select case (numrad_snw)
! case (5)
! nir_bnd_bgn = 2
! case (480)
! nir_bnd_bgn = 51
! end select
! nir_bnd_end    = numrad_snw 

! initialize for adding-doubling solver
allocate(difgauspt(ngmax))
allocate(difgauswt(ngmax))
difgauspt(:) = &  ! gaussian angles (radians)
 (/ 0.9894009,  0.9445750, &
    0.8656312,  0.7554044, &
    0.6178762,  0.4580168, &
    0.2816036,  0.0950125/)
difgauswt(:) = &  ! gaussian weights
 (/ 0.0271525,  0.0622535, &
    0.0951585,  0.1246290, &
    0.1495960,  0.1691565, &
    0.1826034,  0.1894506/)

! initialize for nonspherical snow grains
! sno_shp(:) = snicar_snw_shape ! currently only assuming same shapes for all snow layers
! sno_fs(:)  = 0.0
! sno_AR(:)  = 0.0
!!! EZSNOW: read prognostic snow grain size, and make variable layer-to-layer
        ! idxshp = 1 => SPHERE
    ! idxshp = 1 => SPHEROID
    ! idxshp = 3 => HEXAGONAL
    ! idxshp = 4 => KOCH

    ! note: grain shap array is for each layer, no for each column, as in lm4p2 we run single col here anyway.
do i=-nlevsno+1,0,1
! do i=snl_top,snl_btm,1
   if ((snw_shp_input(1,i) < 1).or.(snw_shp_input(1,i) > 4)) then
      write(*,*) "detected snow shape out of bounds :: = ", snw_shp_input(1,i)
      call land_error_message("SNICAR_RT_HE in snicar_mod: Snow shape error, value out of bounds!", severity=FATAL)
   endif
   ! if (snw_shp_input(i) == 0) sno_shp(i) = snow_shape_defined
   if (snw_shp_input(1,i) == 1) sno_shp(i) = 'sphere' 
   if (snw_shp_input(1,i) == 2) sno_shp(i) = 'spheroid' 
   if (snw_shp_input(1,i) == 3) sno_shp(i) = 'hexagonal_plate' 
   if (snw_shp_input(1,i) == 4) sno_shp(i) = 'koch_snowflake' 
enddo
! sno_shp(:) = 'sphere'
sno_fs(:)  = 0.0 ! ASK
sno_AR(:)  = 0.0 ! ASK
!!! end EZSNOW addition

! Table 3 of He et al 2017 JC
g_wvl(1:8)    = (/ 0.25, 0.70, 1.41, 1.90, &
    2.50, 3.50, 4.00, 5.00 /)
g_wvl_ct(1:7) = g_wvl(2:8) * 0.5 + g_wvl(1:7) * 0.5
g_b0(1:7)     = (/  9.76029E-1,  9.67798E-1,  1.00111, 1.00224, &
     9.64295E-1,  9.97475E-1,  9.97475E-1 /)
g_b1(1:7)     = (/  5.21042E-1,  4.96181E-1,  1.83711E-1,  1.37082E-1, &
     5.50598E-2,  8.48743E-2,  8.48743E-2 /)
g_b2(1:7)     = (/ -2.66792E-4,  1.14088E-3,  2.37011E-4, -2.35905E-4, &
     8.40449E-4, -4.71484E-4, -4.71484E-4 /)
! Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 JAS
g_F07_c2(1:7) = (/  1.349959E-1,  1.115697E-1,  9.853958E-2,  5.557793E-2, &
    -1.233493E-1,  0.0        ,  0.0         /)
g_F07_c1(1:7) = (/ -3.987320E-1, -3.723287E-1, -3.924784E-1, -3.259404E-1, &
     4.429054E-2, -1.726586E-1, -1.726586E-1 /)
g_F07_c0(1:7) = (/  7.938904E-1,  8.030084E-1,  8.513932E-1,  8.692241E-1, &
     7.085850E-1,  6.412701E-1,  6.412701E-1 /)
g_F07_p2(1:7) = (/  3.165543E-3,  2.014810E-3,  1.780838E-3,  6.987734E-4, &
    -1.882932E-2, -2.277872E-2, -2.277872E-2 /)
g_F07_p1(1:7) = (/  1.140557E-1,  1.143152E-1,  1.143814E-1,  1.071238E-1, &
     1.353873E-1,  1.914431E-1,  1.914431E-1 /)
g_F07_p0(1:7) = (/  5.292852E-1,  5.425909E-1,  5.601598E-1,  6.023407E-1, &
     6.473899E-1,  4.634944E-1,  4.634944E-1 /)

! initialize for BC-snow internal mixing
! Eq. 8b & Table 4 in He et al., 2017 J. Climate (wavelength>1.2um, no BC-snow int mixing effect)
bcint_wvl(1:17) = (/ 0.20, 0.25, 0.30, 0.33, 0.36, 0.40, 0.44, 0.48, &
      0.52, 0.57, 0.64, 0.69, 0.75, 0.78, 0.87, 1.0, 1.2 /)
bcint_wvl_ct(1:16) = bcint_wvl(2:17) * 0.5 + bcint_wvl(1:16) * 0.5
bcint_d0(1:16)  = (/ 2.48045   , 4.70305   , 4.68619   , 4.67369   , 4.65040   , &
      2.40364   , 7.95408E-1, 2.92745E-1, 8.63396E-2, 2.76299E-2, &
      1.40864E-2, 8.65705E-3, 6.12971E-3, 4.45697E-3, 3.06648E-2, &
      7.96544E-1 /)
bcint_d1(1:16)  = (/ 9.77209E-1, 9.73317E-1, 9.79650E-1, 9.84579E-1, 9.93537E-1, &
      9.95955E-1, 9.95218E-1, 9.74284E-1, 9.81193E-1, 9.81239E-1, &
      9.55515E-1, 9.10491E-1, 8.74196E-1, 8.27238E-1, 4.82870E-1, &
      4.36649E-2 /)
bcint_d2(1:16)  = (/ 3.95960E-1, 2.04820E-1, 2.07410E-1, 2.09390E-1, 2.13030E-1, &
      4.18570E-1, 1.29682   , 3.75514   , 1.27372E+1, 3.93293E+1, &
      8.78918E+1, 1.86969E+2, 3.45600E+2, 7.08637E+2, 1.41067E+3, &
      2.57288E+2 /)
! Eq. 1a,1b and Table S1 in He et al. 2018 GRL
bcint_m(1:3)    = (/ -0.8724, -0.1866, -0.0046 /)
bcint_n(1:3)    = (/ -0.0072, -0.1918, -0.5177 /)

! initialize for dust-snow internal mixing
! Eq. 1 and Table 1 in He et al. 2019 JAMES (wavelength>1.2um, no dust-snow int mixing effect)
dstint_wvl(1:7) = (/ 0.2, 0.2632, 0.3448, 0.4415, 0.625, 0.7782, 1.2422/)
dstint_wvl_ct(1:6) = dstint_wvl(2:7) * 0.5 + dstint_wvl(1:6) * 0.5
dstint_a1(1:6) = (/ -2.1307E+1, -1.5815E+1, -9.2880   , 1.1115   , 1.0307   , 1.0185    /)
dstint_a2(1:6) = (/  1.1746E+2,  9.3241E+1,  4.0605E+1, 3.7389E-1, 1.4800E-2, 2.8921E-4 /)
dstint_a3(1:6) = (/  9.9701E-1,  9.9781E-1,  9.9848E-1, 1.0035   , 1.0024   , 1.0356    /)

! SNICAR/CLM snow band center wavelength (um)
allocate(wvl_ct(numrad_snw))
! select case (numrad_snw)
! case (5)
wvl_ct(:)  = (/ 0.5, 0.85, 1.1, 1.35, 3.25 /)  ! 5-band
! case (480)
! do igb = 1, snicar_numrad_snw
! wvl_ct(igb) = 0.205 + 0.01 * (igb - 1.0)  ! 480-band
! enddo
! end select

! Define constants
! pi = SHR_CONST_PI ! EZSNOW commented

! always use Delta approximation for snow
DELTA = 1

! Get current timestep
! nstep = get_nstep()
nstep = 1 ! EZSNOW :: NOT USED in LM4p2

! Loop over all non-urban columns
! (when called from CSIM, there is only one column)
do fc = 1,num_nourbanc
! c_idx = filter_nourbanc(fc)
c_idx = fc ! EZSNOW -> single column

! Zero absorbed radiative fluxes:
do i=-nlevsno+1,1,1
flx_abs_lcl(i,:)   = 0.0
flx_abs(c_idx,i,:) = 0.0
enddo

! set snow/ice mass to be used for RT:
if (flg_snw_ice == 1) then
h2osno_lcl = h2osno_total(c_idx)
else
h2osno_lcl = h2osno_ice(c_idx,0)
endif

! Qualifier for computing snow RT: 
!  1) sunlight from atmosphere model 
!  2) minimum amount of snow on ground. 
!     Otherwise, set snow albedo to zero
if ((coszen(c_idx) > 0.0) .and. (h2osno_lcl > min_snw)) then     

! Set variables specific to CLM
if (flg_snw_ice == 1) then
! If there is snow, but zero snow layers, we must create a layer locally.
! This layer is presumed to have the fresh snow effective radius.
if (snl(c_idx) > -1) then 
flg_nosnl         =  1
snl_lcl           =  -1
h2osno_ice_lcl(0) =  h2osno_lcl
h2osno_liq_lcl(0) =  0.0
snw_rds_lcl(0)    =  snw_rds_min_int
else
flg_nosnl         =  0
snl_lcl           =  snl(c_idx)
h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
snw_rds_lcl(:)    =  snw_rds(c_idx,:)
endif

snl_btm   = 0
snl_top   = snl_lcl+1

! for debugging only ! EZSNOW commented
! l_idx     = col%landunit(c_idx)
! g_idx     = col%gridcell(c_idx)
! sfctype   = lun%itype(l_idx)
! lat_coord = grc%latdeg(g_idx)
! lon_coord = grc%londeg(g_idx)

! Set variables specific to CSIM
else
flg_nosnl         = 0
snl_lcl           = -1
h2osno_liq_lcl(:) = h2osno_liq(c_idx,:)
h2osno_ice_lcl(:) = h2osno_ice(c_idx,:)
snw_rds_lcl(:)    = snw_rds(c_idx,:)
snl_btm           = 0
snl_top           = 0
sfctype           = -1
lat_coord         = -90
lon_coord         = 0
endif ! end if flg_snw_ice == 1


! Set local aerosol array
do j=1,sno_nbr_aer
mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
enddo


! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
albsfc_lcl(1:(nir_bnd_bgn-1))       = albsfc(c_idx,1)
albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,2)


! Error check for snow grain size:
do i=snl_top,snl_btm,1
if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
write (*,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
write (*,*) "NSTEP= ", nstep
write (*,*) "flg_snw_ice= ", flg_snw_ice
write (*,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
! write (*,*) "lat= ", lat_coord, " lon= ", lon_coord
write (*,*) "h2osno_total(c)= ", h2osno_lcl
! call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
call land_error_message("SNICAR_RT_HE in snicar_mod: Snow grain radius out of bounds!", severity=FATAL)
endif
enddo


! Incident flux weighting parameters
!  - sum of all VIS bands must equal 1
!  - sum of all NIR bands must equal 1
!
! Spectral bands (5-band case)
!  Band 1: 0.3-0.7um (VIS)
!  Band 2: 0.7-1.0um (NIR)
!  Band 3: 1.0-1.2um (NIR)
!  Band 4: 1.2-1.5um (NIR)
!  Band 5: 1.5-5.0um (NIR)
!
! Hyperspectral (10-nm) bands (480-band case)
! Bands 1~50  : 0.2-0.7um (VIS)
! Bands 51~480: 0.7~5.0um (NIR)
!
! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
!
!           ! works for both 5-band & 480-band, flux weights directly read from input data
! Direct: 
! if (flg_slr_in == 1) then
! ! flx_wgt(1:numrad_snw) = flx_wgt_dir(1:numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
! flx_wgt(1:numrad_snw) = flx_wgt_dir(1,1:numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
! ! Diffuse:
! elseif (flg_slr_in == 2) then
! ! flx_wgt(1:numrad_snw) = flx_wgt_dif(1:numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
! flx_wgt(1:numrad_snw) = flx_wgt_dif(1,1,1:numrad_snw)  ! VIS or NIR band sum is already normalized to 1.0 in input data
! endif

! EZSNOW : until we read the correct data, and decide type of atm to use, use the old defult values::
! These are mid-latitude winter, from SNICAR_RT
if (flg_slr_in == 1) then ! direct 
     flx_wgt(1) = 1.0
     flx_wgt(2) = 0.49352158521175
     flx_wgt(3) = 0.18099494230665
     flx_wgt(4) = 0.12094898498813
     flx_wgt(5) = 0.20453448749347
elseif (flg_slr_in == 2) then ! diffuse
   flx_wgt(1) = 1.0
   flx_wgt(2) = 0.58581507618433
   flx_wgt(3) = 0.20156903770812
   flx_wgt(4) = 0.10917889346386
   flx_wgt(5) = 0.10343699264369
endif

! flx_wgt(1) = 1.
! flx_wgt(2) = 0.58581507618433
! flx_wgt(3) = 0.20156903770812
! flx_wgt(4) = 0.10917889346386
! flx_wgt(5) = 0.10343699264369

exp_min = exp(-argmax)

! Loop over snow spectral bands
do bnd_idx = 1,numrad_snw
! flg_dover is not used since this algorithm is stable for mu_not > 0.01
! mu_not is cosine solar zenith angle above the fresnel level; make
! sure mu_not is large enough for stable and meaningful radiation
! solution: .01 is like sun just touching horizon with its lower edge
! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
mu_not = max(coszen(c_idx), cp01)

flg_dover = 1    ! default is to redo
err_idx   = 0    ! number of times through loop

do while (flg_dover > 0)

! Set direct or diffuse incident irradiance to 1
! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
if (flg_slr_in == 1) then
flx_slrd_lcl(bnd_idx) = 1.0/(mu_not*PI) ! this corresponds to incident irradiance of 1.0
flx_slri_lcl(bnd_idx) = 0.0
else
flx_slrd_lcl(bnd_idx) = 0.0
flx_slri_lcl(bnd_idx) = 1.0
endif

! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
mss_cnc_aer_lcl(:,:) = 0.0
endif

if ( (numrad_snw == 480).and.(bnd_idx > 100) ) then ! >1.2um
mss_cnc_aer_lcl(:,:) = 0.0
endif


!--------------------------- Start snow & aerosol optics --------------------------------
! Define local Mie parameters based on snow grain size and aerosol species retrieved from a lookup table.

! Spherical snow: single-scatter albedo, mass extinction coefficient, asymmetry factor
if (flg_slr_in == 1) then
do i=snl_top,snl_btm,1
   rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
   ! snow optical properties (direct radiation)
   ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
   ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
   if (sno_shp(i) == 'sphere') asm_prm_snw_lcl(i) = asm_prm_snw_drc(rds_idx,bnd_idx)
enddo
elseif (flg_slr_in == 2) then
do i=snl_top,snl_btm,1
   rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
   ! snow optical properties (diffuse radiation)
   ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
   ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
   if (sno_shp(i) == 'sphere') asm_prm_snw_lcl(i) = asm_prm_snw_dfs(rds_idx,bnd_idx)
enddo
endif

! Nonspherical snow: shape-dependent asymmetry factors
do i=snl_top,snl_btm,1

select case (sno_shp(i))
case ('spheroid')
   diam_ice = 2.0 * snw_rds_lcl(i)   ! unit: microns
   if (sno_fs(i) == 0.0) then
      fs_sphd = 0.929  ! default; He et al. (2017), Table 1
   else
      fs_sphd = sno_fs(i) ! user specified value
   endif
   fs_hex = 0.788      ! reference shape factor
   if (sno_AR(i) == 0.0) then
      AR_tmp = 0.5     ! default; He et al. (2017), Table 1
   else
      AR_tmp = sno_AR(i)  ! user specified value
   endif
   do igb = 1,7
      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_sphd/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
      gg_ice_F07_tmp(igb) = g_F07_c0(igb) + g_F07_c1(igb) * AR_tmp + g_F07_c2(igb) * (AR_tmp * AR_tmp)  ! Eqn. 3.1 in Fu (2007)
   enddo

case ('hexagonal_plate')
   diam_ice = 2.0 * snw_rds_lcl(i)   ! unit: microns
   if (sno_fs(i) == 0.0) then
      fs_hex0 = 0.788  ! default; He et al. (2017), Table 1
   else
      fs_hex0 = sno_fs(i) ! user specified value
   endif
   fs_hex = 0.788      ! reference shape factor
   if (sno_AR(i) == 0.0) then
      AR_tmp = 2.5     ! default; He et al. (2017), Table 1
   else
      AR_tmp = sno_AR(i)  ! user specified value
   endif
   do igb = 1,7
      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_hex0/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
      gg_ice_F07_tmp(igb) = g_F07_p0(igb) + g_F07_p1(igb) * log(AR_tmp) + g_F07_p2(igb) * (log(AR_tmp) * log(AR_tmp)) ! Eqn. 3.3 in Fu (2007)
   enddo

case ('koch_snowflake')
   diam_ice = 2.0 * snw_rds_lcl(i) / 0.544  ! unit: microns
   if (sno_fs(i) == 0.0) then
      fs_koch = 0.712  ! default; He et al. (2017), Table 1
   else
      fs_koch = sno_fs(i) ! user specified value
   endif
   fs_hex = 0.788      ! reference shape factor
   if (sno_AR(i) == 0.0) then
      AR_tmp = 2.5     ! default; He et al. (2017), Table 1
   else
      AR_tmp = sno_AR(i)  ! user specified value
   endif
   do igb = 1,7
      g_ice_Cg_tmp(igb) = g_b0(igb) * ((fs_koch/fs_hex)**g_b1(igb)) * (diam_ice**g_b2(igb))   ! Eq.7, He et al. (2017)
      gg_ice_F07_tmp(igb) = g_F07_p0(igb) + g_F07_p1(igb) * log(AR_tmp) + g_F07_p2(igb) * (log(AR_tmp) * log(AR_tmp)) ! Eqn. 3.3 in Fu (2007)
   enddo

end select

! compute nonspherical snow asymmetry factor
if (sno_shp(i) /= 'sphere') then
   ! 7 wavelength bands for g_ice to be interpolated into targeted SNICAR bands here
   ! use the piecewise linear interpolation subroutine created at the end of this module
   ! tests showed the piecewise linear interpolation has similar results as pchip interpolation
   call piecewise_linear_interp1d(7, g_wvl_ct, g_ice_Cg_tmp, wvl_ct(bnd_idx), g_Cg_intp)
   call piecewise_linear_interp1d(7, g_wvl_ct, gg_ice_F07_tmp, wvl_ct(bnd_idx), gg_F07_intp)
   g_ice_F07 = gg_F07_intp + 0.5 * (1.0 - gg_F07_intp) / ss_alb_snw_lcl(i)  ! Eq.2.2 in Fu (2007)
   asm_prm_snw_lcl(i) = g_ice_F07 * g_Cg_intp     ! Eq.6, He et al. (2017)
endif

if (asm_prm_snw_lcl(i) > 0.99) asm_prm_snw_lcl(i) = 0.99 !avoid unreasonable values (rarely occur in large-size spheroid cases)

enddo ! snow layer loop

! aerosol species 2 optical properties, hydrophobic BC
ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)      
asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)

! aerosol species 3 optical properties, hydrophilic OC
ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)      
asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

! aerosol species 4 optical properties, hydrophobic OC
ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)      
asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
! 2. optical Depths (tau_snw, tau_aer)
! 3. weighted Mie properties (tau, omega, g)

! Weighted Mie parameters of each layer
do i=snl_top,snl_btm,1

! Optics for BC/dust-snow external mixing:
! aerosol species 1 optical properties, hydrophilic BC
ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)
! aerosol species 5 optical properties, dust size1
ss_alb_aer_lcl(5)      = ss_alb_dst1(bnd_idx)
asm_prm_aer_lcl(5)     = asm_prm_dst1(bnd_idx)
ext_cff_mss_aer_lcl(5) = ext_cff_mss_dst1(bnd_idx)
! aerosol species 6 optical properties, dust size2
ss_alb_aer_lcl(6)      = ss_alb_dst2(bnd_idx)
asm_prm_aer_lcl(6)     = asm_prm_dst2(bnd_idx)
ext_cff_mss_aer_lcl(6) = ext_cff_mss_dst2(bnd_idx)
! aerosol species 7 optical properties, dust size3
ss_alb_aer_lcl(7)      = ss_alb_dst3(bnd_idx)
asm_prm_aer_lcl(7)     = asm_prm_dst3(bnd_idx)
ext_cff_mss_aer_lcl(7) = ext_cff_mss_dst3(bnd_idx)
! aerosol species 8 optical properties, dust size4
ss_alb_aer_lcl(8)      = ss_alb_dst4(bnd_idx)
asm_prm_aer_lcl(8)     = asm_prm_dst4(bnd_idx)
ext_cff_mss_aer_lcl(8) = ext_cff_mss_dst4(bnd_idx)

! Start BC/dust-snow internal mixing for wavelength<=1.2um
wvl_doint = wvl_ct(bnd_idx)

if (wvl_doint <= 1.2) then

   ! BC-snow internal mixing applied to hydrophilic BC if activated
   ! BC-snow internal mixing primarily affect snow single-scattering albedo
   if ( snicar_snobc_intmix .and. (mss_cnc_aer_lcl(i,1) > 0.0) ) then
      ! result from Eq.8b in He et al.(2017) is based on BC Re=0.1um &
      ! MAC=6.81 m2/g (@550 nm) & BC density=1.7g/cm3.
      ! To be consistent with Bond et al. 2006 recommeded value (BC MAC=7.5 m2/g @550nm)
      ! we made adjustments on BC size & density as follows to get MAC=7.5m2/g:
      ! (1) We use BC Re=0.045um [geometric mean diameter=0.06um (Dentener et al.2006, 
      ! Yu and Luo,2009) & geometric std=1.5 (Flanner et al.2007;Aoki et al., 2011)].
      ! (2) We tune BC density from 1.7 to 1.49 g/cm3 (Aoki et al., 2011).
      ! These adjustments also lead to consistent results with Flanner et al. 2012 (ACP) lookup table
      ! for BC-snow internal mixing enhancement in albedo reduction (He et al. 2018 ACP)
      do ibb=1,16
         enh_omg_bcint_tmp(ibb) = bcint_d0(ibb) * &
            ( (mss_cnc_aer_lcl(i,1)*1.0E9*1.7/den_bc + bcint_d2(ibb)) **bcint_d1(ibb) )
         ! adjust enhancment factor for BC effective size from 0.1um to Re_bc (He et al. 2018 GRL Eqs.1a,1b)
         if (ibb < 3) then ! near-UV
            bcint_m_tmp = bcint_m(1)
            bcint_n_tmp = bcint_n(1)
         else if (ibb >= 3 .and. ibb <= 11) then ! visible
            bcint_m_tmp = bcint_m(2)
            bcint_n_tmp = bcint_n(2)
         else  ! ibb > 11, NIR
            bcint_m_tmp = bcint_m(3)
            bcint_n_tmp = bcint_n(3)
         endif
         bcint_dd  = (Re_bc * 20.0)**bcint_m_tmp
         bcint_dd2 = (0.1 * 20.0)**bcint_m_tmp
         bcint_f  = (Re_bc * 10.0)**bcint_n_tmp

         enh_omg_bcint_tmp2(ibb)=LOG10(max(1.0,bcint_dd*((enh_omg_bcint_tmp(ibb)/bcint_dd2)**bcint_f)))
      enddo
      ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
      call piecewise_linear_interp1d(16,bcint_wvl_ct,enh_omg_bcint_tmp2,wvl_doint,enh_omg_bcint_intp)
      ! update snow single-scattering albedo
      enh_omg_bcint_intp2 = 10.0 ** enh_omg_bcint_intp                           
      enh_omg_bcint_intp2 = min(1.0E5, max(enh_omg_bcint_intp2,1.0)) ! constrain enhancement to a reasonable range
      ss_alb_snw_lcl(i)   = 1.0 - (1.0 - ss_alb_snw_lcl(i)) * enh_omg_bcint_intp2
      ss_alb_snw_lcl(i)   = max(0.5, min(ss_alb_snw_lcl(i),1.0))
      ! reset hydrophilic BC property to 0 since it is accounted by updated snow ss_alb above
      ss_alb_aer_lcl(1)       = 0.0
      asm_prm_aer_lcl(1)      = 0.0
      ext_cff_mss_aer_lcl(1)  = 0.0
   endif ! end if BC-snow mixing type

   ! Dust-snow internal mixing applied to all size bins if activated
   ! Dust-snow internal mixing primarily affect snow single-scattering albedo
   ! default optics of externally mixed dust at 4 size bins based on effective
   ! radius of 1.38um and sigma=2.0 with truncation to each size bin (Flanner et al. 2021 GMD)
   ! parameterized dust-snow int mix results based on effective radius of 1.1um and sigma=2.0
   ! from (He et al. 2019 JAMES). Thus, the parameterization can be approximately applied to
   ! all dust size bins here.
   tot_dst_snw_conc = (mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) + &
                       mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8)) * 1.0E6 !kg/kg->ppm
   if ( snicar_snodst_intmix .and. (tot_dst_snw_conc > 0.0) ) then
      do idb=1,6
         enh_omg_dstint_tmp(idb) = dstint_a1(idb)+dstint_a2(idb)*(tot_dst_snw_conc**dstint_a3(idb))
         enh_omg_dstint_tmp2(idb) = LOG10(max(enh_omg_dstint_tmp(idb),1.0))
      enddo
      ! piecewise linear interpolate into targeted SNICAR bands in a logscale space
      call piecewise_linear_interp1d(6,dstint_wvl_ct,enh_omg_dstint_tmp2,wvl_doint,enh_omg_dstint_intp)
      ! update snow single-scattering albedo
      enh_omg_dstint_intp2 = 10.0 ** enh_omg_dstint_intp
      enh_omg_dstint_intp2 = min(1.0E5, max(enh_omg_dstint_intp2,1.0)) ! constrain enhancement to a reasonable range
      ss_alb_snw_lcl(i) = 1.0 - (1.0 - ss_alb_snw_lcl(i)) * enh_omg_dstint_intp2
      ss_alb_snw_lcl(i) = max(0.5, min(ss_alb_snw_lcl(i),1.0))
      ! reset all dust optics to zero  since it is accounted by updated snow ss_alb above
      ss_alb_aer_lcl(5:8)      = 0.0
      asm_prm_aer_lcl(5:8)     = 0.0
      ext_cff_mss_aer_lcl(5:8) = 0.0
   endif ! end if dust-snow internal mixing

endif ! end if BC/dust-snow internal mixing (bands<1.2um)

L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

do j=1,sno_nbr_aer
   L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
   tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
enddo

tau_sum   = 0.0
omega_sum = 0.0
g_sum     = 0.0

do j=1,sno_nbr_aer
   tau_sum    = tau_sum + tau_aer(i,j) 
   omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
   g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
enddo

tau(i)    = tau_sum + tau_snw(i)
omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))

enddo  ! end do snow layers

! DELTA transformations, if requested
if (DELTA == 1) then
do i=snl_top,snl_btm,1
   g_star(i)     = g(i)/(1+g(i))
   omega_star(i) = (1.0 - g(i) * g(i)) * omega(i) / (1.0 - omega(i) * (g(i) * g(i)))
   tau_star(i)   = (1.0 - omega(i) * (g(i) * g(i))) * tau(i)
enddo
else
do i=snl_top,snl_btm,1
   g_star(i)     = g(i)
   omega_star(i) = omega(i)
   tau_star(i)   = tau(i)
enddo
endif
!--------------------------- End of snow & aerosol optics --------------------------------

!--------------------------- Start Adding-doubling RT solver  --------------------------------

! Given input vertical profiles of optical properties, evaluate the
! monochromatic Delta-Eddington adding-doubling solution

! trndir, trntdr, trndif, rupdir, rupdif, rdndif are variables at the layer interface,
! for snow with layers from snl_top to snl_btm there are snl_top to snl_btm+1 layer interface
snl_btm_itf = snl_btm + 1

! initialization for layer interface
do i = snl_top,snl_btm_itf,1
trndir(i) = c0
trntdr(i) = c0
trndif(i) = c0
rupdir(i) = c0
rupdif(i) = c0
rdndif(i) = c0
enddo
! initialize top interface of top layer
trndir(snl_top) = c1
trntdr(snl_top) = c1
trndif(snl_top) = c1
rdndif(snl_top) = c0

! begin main level loop for snow layer interfaces except for the very bottom
do i = snl_top,snl_btm,1

! initialize all layer apparent optical properties to 0
rdir  (i) = c0
rdif_a(i) = c0
rdif_b(i) = c0
tdir  (i) = c0
tdif_a(i) = c0
tdif_b(i) = c0
trnlay(i) = c0

! compute next layer Delta-eddington solution only if total transmission
! of radiation to the interface just above the layer exceeds trmin.
if (trntdr(i) > trmin ) then

 ! delta-transformed single-scattering properties of this layer
 ts = tau_star(i)
 ws = omega_star(i)
 gs = g_star(i)

 ! Delta-Eddington solution expressions, Eq. 50: Briegleb and Light 2007
 lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))
 ue = c1p5*(c1 - ws*gs)/lm
 extins = max(exp_min, exp(-lm*ts))
 ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins)

 ! first calculation of rdif, tdif using Delta-Eddington formulas
 ! Eq.: Briegleb 1992; alpha and gamma for direct radiation
 rdif_a(i) = (ue * ue - c1) * (c1 / extins - extins) / ne
 tdif_a(i) = c4*ue/ne

 ! evaluate rdir,tdir for direct beam
 trnlay(i) = max(exp_min, exp(-ts/mu_not))

 ! Delta-Eddington solution expressions
 ! Eq. 50: Briegleb and Light 2007; alpha and gamma for direct radiation
 alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
 gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
 apg = alp + gam
 amg = alp - gam
 rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
 tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)

 ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
 ! since Delta-Eddington rdif formula is not well-behaved (it is usually
 ! biased low and can even be negative); use ngmax angles and gaussian
 ! integration for most accuracy:
 R1 = rdif_a(i) ! use R1 as temporary
 T1 = tdif_a(i) ! use T1 as temporary
 swt = c0
 smr = c0
 smt = c0
 ! gaussian angles for the AD integral
 do ng=1,ngmax
    mu  = difgauspt(ng)
    gwt = difgauswt(ng)
    swt = swt + mu*gwt
    trn = max(exp_min, exp(-ts/mu))
    alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
    gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
    apg = alp + gam
    amg = alp - gam
    rdr = apg*R1 + amg*T1*trn - amg
    tdr = apg*T1 + amg*R1*trn - apg*trn + trn
    smr = smr + mu*rdr*gwt
    smt = smt + mu*tdr*gwt
 enddo      ! ng
 rdif_a(i) = smr/swt
 tdif_a(i) = smt/swt

 ! homogeneous layer
 rdif_b(i) = rdif_a(i)
 tdif_b(i) = tdif_a(i)

endif ! trntdr(k) > trmin

! Calculate the solar beam transmission, total transmission, and
! reflectivity for diffuse radiation from below at interface i,
! the top of the current layer k:
!
!              layers       interface
!
!       ---------------------  i-1
!                i-1
!       ---------------------  i
!                 i
!       ---------------------

trndir(i+1) = trndir(i)*trnlay(i)            ! solar beam transmission from top
refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))  ! interface multiple scattering for i-1
tdrrdir     = trndir(i)*rdir(i)              ! direct tran times layer direct ref
tdndif      = trntdr(i) - trndir(i)          ! total down diffuse = tot tran - direct tran
trntdr(i+1) = trndir(i)*tdir(i) + &          ! total transmission to direct beam for layers above
             (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
! Eq. B4; Briegleb and Light 2007
rdndif(i+1) = rdif_b(i) + &                  ! reflectivity to diffuse radiation for layers above
             (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
trndif(i+1) = trndif(i)*refkm1*tdif_a(i)     ! diffuse transmission to diffuse beam for layers above

enddo       ! end i main level loop

! compute reflectivity to direct and diffuse radiation for layers
! below by adding succesive layers starting from the underlying
! ground and working upwards:
!
!              layers       interface
!
!       ---------------------  i
!                 i
!       ---------------------  i+1
!                i+1
!       ---------------------

! set the underlying ground albedo == albedo of near-IR
! unless bnd_idx < nir_bnd_bgn, for visible
rupdir(snl_btm_itf) = albsfc(c_idx,2)
rupdif(snl_btm_itf) = albsfc(c_idx,2)
if (bnd_idx < nir_bnd_bgn) then
rupdir(snl_btm_itf) = albsfc(c_idx,1)
rupdif(snl_btm_itf) = albsfc(c_idx,1)
endif

do i=snl_btm,snl_top,-1
! interface scattering Eq. B5; Briegleb and Light 2007
refkp1 = c1/( c1 - rdif_b(i)*rupdif(i+1))
! dir from top layer plus exp tran ref from lower layer, interface
! scattered and tran thru top layer from below, plus diff tran ref
! from lower layer with interface scattering tran thru top from below
rupdir(i) = rdir(i) &
           + (        trnlay(i)  *rupdir(i+1) &
           +  (tdir(i)-trnlay(i))*rupdif(i+1) ) * refkp1 * tdif_b(i)
! dif from top layer from above, plus dif tran upwards reflected and
! interface scattered which tran top from below
rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
enddo       ! i

! net flux (down-up) at each layer interface from the
! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
! the interface reflectivities and transmissivities required
! to evaluate interface fluxes are returned from solution_dEdd;
! now compute up and down fluxes for each interface, using the
! combined layer properties at each interface:
!
!              layers       interface
!
!       ---------------------  i
!                 i
!       ---------------------

do i = snl_top, snl_btm_itf
! interface scattering, Eq. 52; Briegleb and Light 2007
refk = c1/(c1 - rdndif(i)*rupdif(i))
! dir tran ref from below times interface scattering, plus diff
! tran and ref from below times interface scattering
! fdirup(i) = (trndir(i)*rupdir(i) + &
!                 (trntdr(i)-trndir(i))  &
!                 *rupdif(i))*refk
! dir tran plus total diff trans times interface scattering plus
! dir tran with up dir ref and down dif ref times interface scattering
! fdirdn(i) = trndir(i) + (trntdr(i) &
!               - trndir(i) + trndir(i)  &
!               *rupdir(i)*rdndif(i))*refk
! diffuse tran ref from below times interface scattering
! fdifup(i) = trndif(i)*rupdif(i)*refk
! diffuse tran times interface scattering
! fdifdn(i) = trndif(i)*refk

! netflux, down - up
! dfdir = fdirdn - fdirup
dfdir(i) = trndir(i) &
         + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
         -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
if (dfdir(i) < puny) dfdir(i) = c0
! dfdif = fdifdn - fdifup
dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
if (dfdif(i) < puny) dfdif(i) = c0
enddo  ! k

! SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
! direct incident
if (flg_slr_in == 1) then
albedo = rupdir(snl_top)
dftmp  = dfdir
refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
F_sfc_pls = (trndir(snl_top)*rupdir(snl_top) + &
           (trntdr(snl_top)-trndir(snl_top))  &
           *rupdif(snl_top))*refk
!diffuse incident
else
albedo = rupdif(snl_top)
dftmp  = dfdif
refk   = c1/(c1 - rdndif(snl_top)*rupdif(snl_top))
F_sfc_pls = trndif(snl_top)*rupdif(snl_top)*refk
endif

! Absorbed flux in each layer
do i=snl_top,snl_btm,1
F_abs(i) = dftmp(i)-dftmp(i+1)
flx_abs_lcl(i,bnd_idx) = F_abs(i)

! ERROR check: negative absorption
if (flx_abs_lcl(i,bnd_idx) < -0.00001) then
 write (*,"(a,e13.6,a,i6,a,i6)") "SNICAR ERROR: negative absoption : ", &
       flx_abs_lcl(i,bnd_idx), " at timestep: ", nstep, " at column: ", c_idx
 write(*,*) "SNICAR_AD STATS: snw_rds(0)= ", snw_rds(c_idx,0)
 write(*,*) "SNICAR_AD STATS: L_snw(0)= ", L_snw(0)
 write(*,*) "SNICAR_AD STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
 write(*,*) "SNICAR_AD STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
 write(*,*) "SNICAR_AD STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
 write(*,*) "SNICAR_AD STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
 write(*,*) "SNICAR_AD STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
 write(*,*) "SNICAR_AD STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
 write(*,*) "SNICAR_AD STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
!  call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
 call land_error_message("SNICAR_RT_HE in snicar_mod: Negative absorption!", severity=FATAL)
endif
enddo

! absobed flux by the underlying ground
F_btm_net = dftmp(snl_btm_itf)

! note here, snl_btm_itf = 1 by snow column set up in CLM
flx_abs_lcl(1,bnd_idx) = F_btm_net

if (flg_nosnl == 1) then
! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
!flx_abs_lcl(:,bnd_idx) = 0._r8
!flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

! changed on 20070408:
! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
! handles the case of no snow layers. Then, if a snow layer is addded between now and
! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
flx_abs_lcl(0,bnd_idx) = F_abs(0)
flx_abs_lcl(1,bnd_idx) = F_btm_net
endif

!Underflow check (we've already tripped the error condition above)
do i=snl_top,1,1
flx_abs_lcl(i,bnd_idx) = max(0.0, flx_abs_lcl(i,bnd_idx))
enddo

F_abs_sum = 0.0
do i=snl_top,snl_btm,1
F_abs_sum = F_abs_sum + F_abs(i)
enddo

! no need to repeat calculations for adding-doubling solver
flg_dover = 0

!--------------------------- End of Adding-doubling RT solver  --------------------------------

enddo !enddo while (flg_dover > 0)

! Energy conservation check:
! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
energy_sum = (mu_not*PI*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
if (abs(energy_sum) > 0.00001) then
write (*,"(a,e12.6,a,i6,a,i6)") "SNICAR ERROR: Energy conservation error of : ", energy_sum, &
  " at timestep: ", nstep, " at column: ", c_idx
write(*,*) "F_abs_sum: ",F_abs_sum
write(*,*) "F_btm_net: ",F_btm_net
write(*,*) "F_sfc_pls: ",F_sfc_pls
write(*,*) "mu_not*PI*flx_slrd_lcl(bnd_idx): ", mu_not*PI*flx_slrd_lcl(bnd_idx)
write(*,*) "flx_slri_lcl(bnd_idx)", flx_slri_lcl(bnd_idx)
write(*,*) "bnd_idx", bnd_idx
write(*,*) "F_abs", F_abs
write(*,*) "albedo", albedo
call land_error_message("SNICAR_RT_HE in snicar_mod: Energy conservation error!", severity=FATAL)
! call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
endif

albout_lcl(bnd_idx) = albedo

! Fail if albedo > 1
if (albout_lcl(bnd_idx) > 1.0) then

write (*,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx, " NSTEP= ",nstep
write (*,*) "SNICAR STATS: bnd_idx= ",bnd_idx
write (*,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), &
  " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
write (*,*) "SNICAR STATS: landtype= ", sfctype
write (*,*) "SNICAR STATS: h2osno_total= ", h2osno_lcl, " snl= ", snl_lcl
write (*,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in

write (*,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
write (*,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
write (*,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
write (*,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
write (*,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)

write (*,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
write (*,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
write (*,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
write (*,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
write (*,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

write (*,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
write (*,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
write (*,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
write (*,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
write (*,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)

call land_error_message("SNICAR_RT_HE in snicar_mod: Albedo out of bounds!", severity=FATAL)

! call endrun(subgrid_index=c_idx, subgrid_level=subgrid_level_column, msg=errmsg(sourcefile, __LINE__))
endif

enddo   ! loop over wvl bands


! Weight output NIR albedo appropriately
select case (numrad_snw)
case (5)  ! 5-band case
! VIS band
albout(c_idx,1) = albout_lcl(1)
case (480)  ! 480-band case
! average for VIS band
flx_sum = 0.0
do bnd_idx= 1, (nir_bnd_bgn-1)
flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
end do
albout(c_idx,1) = flx_sum / sum(flx_wgt(1:(nir_bnd_bgn-1)))
end select

! average for NIR band (5 or 480-band case)
flx_sum = 0.0
do bnd_idx = nir_bnd_bgn, nir_bnd_end
flx_sum = flx_sum + flx_wgt(bnd_idx) * albout_lcl(bnd_idx)
end do
albout(c_idx,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
select case (numrad_snw)
case (5)  ! 5-band case
! VIS band
flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
case (480)  ! 480-band case
! average for VIS band
do i=snl_top,1,1
flx_sum = 0.0
do bnd_idx= 1,(nir_bnd_bgn-1)
flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
enddo
flx_abs(c_idx,i,1) = flx_sum / sum(flx_wgt(1:(nir_bnd_bgn-1)))
end do
end select

! average for NIR band (5 or 480-band case)
do i = snl_top, 1, 1
flx_sum = 0.0
do bnd_idx = nir_bnd_bgn, nir_bnd_end
flx_sum = flx_sum + flx_wgt(bnd_idx) * flx_abs_lcl(i,bnd_idx)
end do
flx_abs(c_idx,i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
end do

! high solar zenith angle adjustment for Adding-doubling solver results
! near-IR direct albedo/absorption adjustment for high solar zenith angles
! solar zenith angle parameterization
! calculate the scaling factor for NIR direct albedo if SZA>75 degree
if ((mu_not < mu_75) .and. (flg_slr_in == 1)) then
sza_c1 = sza_a0 + sza_a1 * mu_not + sza_a2 * (mu_not * mu_not)
sza_c0 = sza_b0 + sza_b1 * mu_not + sza_b2 * (mu_not * mu_not)
sza_factor = sza_c1 * (log10(snw_rds_lcl(snl_top) * c1) - c6) + sza_c0
flx_sza_adjust  = albout(c_idx,2) * (sza_factor-c1) * sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
albout(c_idx,2) = albout(c_idx,2) * sza_factor
flx_abs(c_idx,snl_top,2) = flx_abs(c_idx,snl_top,2) - flx_sza_adjust
endif


! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
elseif ( (coszen(c_idx) > 0.0) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0.0) ) then
albout(c_idx,1) = albsfc(c_idx,1)
albout(c_idx,2) = albsfc(c_idx,2)

! There is either zero snow, or no sun
else
albout(c_idx,1) = 0.0
albout(c_idx,2) = 0.0
endif    ! if column has snow and coszen > 0

enddo    ! loop over all columns

! end associate

end subroutine SNICAR_RT_HE


   !-----------------------------------------------------------------------
subroutine piecewise_linear_interp1d(nd, xd, yd, xi, yi)

   ! piecewise linear interpolation method for 1-dimensional data
   ! original author: John Burkardt, Florida State University, 09/22/2012
   ! Added and modified by Cenlin He (NCAR), 01/27/2022

   implicit none

   integer , intent(in)   :: nd         ! number of data points of (xd)
   real, intent(in)   :: xd(1:nd)   ! x-value of data points
   real, intent(in)   :: yd(1:nd)   ! y-value of data points
   real, intent(in)   :: xi         ! x-value for to-be-interpolated point
   real, intent(out)  :: yi         ! the interpolated value at xi

   ! local variables
   integer  :: i, k    ! loop index
   real :: t

   yi = 0.0

   ! if only one data point
   if ( nd == 1 ) then
      yi = yd(1)
      return
   endif

   ! if multiple data points
   if ( xi < xd(1) ) then ! extrapolate
      t  = ( xi - xd(1) ) / ( xd(2) - xd(1) )
      yi = (1.0 - t) * yd(1) + t * yd(2)
   elseif ( xi > xd(nd) ) then ! extrapolate
      t  = ( xi - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi = (1.0 - t) * yd(nd-1) + t * yd(nd)
   else  ! piecsewise interpolate
      do k = 2, nd
         if ( (xd(k-1) <= xi) .and. (xi <= xd(k)) ) then
            t  = ( xi - xd(k-1) ) / ( xd(k) - xd(k-1) )
            yi = (1.0 - t) * yd(k-1) + t * yd(k)
            exit
         endif
      enddo
   endif

   return

 end subroutine piecewise_linear_interp1d

  

 end module snicar_mod
