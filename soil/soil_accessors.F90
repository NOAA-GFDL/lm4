!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Land Model 4 (LM4).
!*
!* LM4 is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* LM4 is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with LM4.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module soil_accessors_mod

use land_tile_mod, only: land_tile_type

implicit none
public

contains

subroutine soil_T_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%T(i)
    endif
end subroutine

subroutine soil_wl_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%wl(i)
    endif
end subroutine

subroutine soil_ws_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p;p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%ws(i)
    endif
end subroutine

subroutine soil_groundwater_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%groundwater(i)
    endif
end subroutine

subroutine soil_groundwater_T_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%groundwater_T(i)
    endif
end subroutine

subroutine soil_frozen_freq_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%frozen_freq(i)
    endif
end subroutine

subroutine soil_w_fc_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%w_fc(i)
    endif
end subroutine

subroutine soil_alpha_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%alpha(i)
    endif
end subroutine

subroutine soil_uptake_T_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%uptake_T
    endif
end subroutine

subroutine soil_tag_ptr(t,p)
    type(land_tile_type),pointer::t
    integer,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%tag
        endif
end subroutine

subroutine soil_fast_soil_C_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%fast_soil_C(i)
    endif
end subroutine

subroutine soil_slow_soil_C_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%slow_soil_C(i)
    endif
end subroutine

subroutine soil_asoil_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%asoil_in(i)
    endif
end subroutine

subroutine soil_fsc_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%fsc_in(i)
    endif
end subroutine

subroutine soil_ssc_in_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%ssc_in(i)
    endif
end subroutine

subroutine soil_is_peat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    integer,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%is_peat(i)
    endif
end subroutine

subroutine soil_tau_groundwater_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%tau_groundwater
    endif
end subroutine

subroutine soil_hillslope_length_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_length
    endif
end subroutine

subroutine soil_hillslope_relief_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_relief
    endif
end subroutine

subroutine soil_hillslope_a_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_a
    endif
end subroutine

subroutine soil_hillslope_n_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_n
    endif
end subroutine

subroutine soil_hillslope_zeta_bar_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%hillslope_zeta_bar
    endif
end subroutine

subroutine soil_soil_e_depth_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%soil_e_depth
    endif
end subroutine

subroutine soil_zeta_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%zeta
    endif
end subroutine

subroutine soil_tau_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%tau
    endif
end subroutine

subroutine soil_k_sat_gw_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%k_sat_gw
    endif
end subroutine

subroutine soil_vwc_wilt_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_wilt
    endif
end subroutine

subroutine soil_vwc_fc_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_fc
    endif
end subroutine

subroutine soil_vwc_sat_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%vwc_sat
    endif
end subroutine

subroutine soil_k_sat_ref_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%k_sat_ref
    endif
end subroutine

subroutine soil_Qmax_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%Qmax
    endif
end subroutine

subroutine soil_refl_dry_dir_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_dry_dir(i)
    endif
end subroutine

subroutine soil_refl_dry_dif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_dry_dif(i)
    endif
end subroutine

subroutine soil_refl_sat_dir_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_sat_dir(i)
    endif
end subroutine

subroutine soil_refl_sat_dif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%refl_sat_dif(i)
    endif
end subroutine

subroutine soil_f_iso_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_iso_dry(i)
    endif
end subroutine

subroutine soil_f_vol_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_vol_dry(i)
    endif
end subroutine

subroutine soil_f_geo_dry_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_geo_dry(i)
    endif
end subroutine

subroutine soil_f_iso_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_iso_sat(i)
    endif
end subroutine

subroutine soil_f_vol_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_vol_sat(i)
    endif
end subroutine

subroutine soil_f_geo_sat_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL();
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%pars%f_geo_sat(i)
    endif
end subroutine

subroutine soil_gross_nitrogen_flux_into_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_into_tile
    endif
end subroutine

subroutine soil_gross_nitrogen_flux_out_of_tile_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%gross_nitrogen_flux_out_of_tile
    endif
end subroutine


! stuff below is for CORPSE
subroutine sc_soil_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%litterC(k)
  endif
end subroutine

subroutine sc_negative_litter_C_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(associated(t))then
     if(associated(t%soil))p=>t%soil%neg_litt_C(i)
  endif
end subroutine

subroutine sc_negative_litter_N_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(associated(t))then
     if(associated(t%soil))p=>t%soil%neg_litt_N(i)
  endif
end subroutine

subroutine sc_soil_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%litterN(k)
  endif
end subroutine

subroutine sc_protected_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%protectedC(k)
  endif
end subroutine

subroutine sc_protected_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%protectedN(k)
  endif
end subroutine

subroutine sc_DOC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%dissolved_carbon(j)
  endif
end subroutine

subroutine sc_DON_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%dissolved_nitrogen(j)
  endif
end subroutine

subroutine sc_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%C_in(j)
  endif
end subroutine

subroutine sc_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%N_in(j)
  endif
end subroutine

subroutine sc_protected_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%protected_C_in(j)
  endif
end subroutine

subroutine sc_protected_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%protected_N_in(j)
  endif
end subroutine

subroutine sc_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%C_turnover(j)
  endif
end subroutine

subroutine sc_protected_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%protected_C_turnover(j)
  endif
end subroutine

subroutine sc_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%N_turnover(j)
  endif
end subroutine

subroutine sc_protected_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%protected_N_turnover(j)
  endif
end subroutine

subroutine sc_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%nitrate
  endif
end subroutine

subroutine sc_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%ammonium
  endif
end subroutine

subroutine sc_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%nitrif
  endif
end subroutine

subroutine sc_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(associated(t)) then
     if(associated(t%soil))p=>t%soil%org_matter(i)%denitrif
  endif
end subroutine

subroutine sc_litter_nitrate_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%nitrate
    endif
end subroutine

subroutine sc_litter_ammonium_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%ammonium
    endif
end subroutine

subroutine sc_litter_nitrif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%nitrif
    endif
end subroutine

subroutine sc_litter_denitrif_ptr(t,i,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(i)%denitrif
    endif
end subroutine

subroutine sc_litter_C_in_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%C_in(i)
    endif
end subroutine

subroutine sc_litter_N_in_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%N_in(i)
    endif
end subroutine

subroutine sc_litter_C_turnover_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%C_turnover(i)
    endif
end subroutine

subroutine sc_litter_N_turnover_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%N_turnover(i)
    endif
end subroutine

subroutine sc_litter_dissolved_carbon_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%dissolved_carbon(i)
    endif
end subroutine

subroutine sc_litter_dissolved_nitrogen_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%dissolved_nitrogen(i)
    endif
end subroutine

subroutine sc_litter_livingMicrobeC_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%livingMicrobeC
    endif
end subroutine

subroutine sc_litter_livingMicrobeN_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%livingMicrobeN
    endif
end subroutine

subroutine sc_litter_CO2_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j
    real,pointer::p
    p=>NULL()
    if(associated(t)) then
        if(associated(t%soil))p=>t%soil%litter(j)%litterCohorts(i)%CO2
    endif
end subroutine

subroutine sc_litter_litterC_ptr(t,i,j,k,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j,k
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(k)%litterCohorts(i)%litterC(j)
    endif
end subroutine

subroutine sc_litter_protectedC_ptr(t,i,j,k,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j,k
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(k)%litterCohorts(i)%protectedC(j)
    endif
end subroutine

subroutine sc_litter_litterN_ptr(t,i,j,k,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j,k
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(k)%litterCohorts(i)%litterN(j)
    endif
end subroutine

subroutine sc_litter_protectedN_ptr(t,i,j,k,p)
    type(land_tile_type),pointer::t
    integer,intent(in)::i,j,k
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%litter(k)%litterCohorts(i)%protectedN(j)
    endif
end subroutine

subroutine sc_livingMicrobeC_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    integer,intent(in)::i,j
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%livingMicrobeC
    endif
end subroutine

subroutine sc_CO2_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    integer,intent(in)::i,j
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%CO2
    endif
end subroutine

subroutine sc_livingMicrobeN_ptr(t,i,j,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    integer,intent(in)::i,j
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%org_matter(i)%litterCohorts(j)%livingMicrobeN
    endif
end subroutine

end module soil_accessors_mod
