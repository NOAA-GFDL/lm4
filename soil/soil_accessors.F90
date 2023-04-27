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
        if(associated(t%soil))p=>t%soil%w_wilt(1)
    endif
end subroutine

subroutine soil_vwc_fc_ptr(t,p)
    type(land_tile_type),pointer::t
    real,pointer::p
    p=>NULL()
    if(associated(t))then
        if(associated(t%soil))p=>t%soil%w_fc(1)
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

end module soil_accessors_mod