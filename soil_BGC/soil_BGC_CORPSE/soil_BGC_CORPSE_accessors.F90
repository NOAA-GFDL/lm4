module soilc_CORPSE_accessors_mod

use land_tile_mod, only: land_tile_type
use soil_BGC_CORPSE_type_mod, only: soil_BGC_CORPSE_t

implicit none
public

contains

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

subroutine soil_is_peat_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; integer,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%is_peat(i)
  end select
end subroutine

subroutine sc_soil_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%litterC(k)
  end select
end subroutine

subroutine sc_soil_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%litterN(k)
  end select
end subroutine

subroutine sc_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%nitrate
  end select
end subroutine

subroutine sc_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%ammonium
  end select
end subroutine

subroutine sc_livingMicrobeC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%livingMicrobeC
  end select
end subroutine

subroutine sc_livingMicrobeN_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%livingMicrobeN
  end select
end subroutine

subroutine sc_CO2_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; real,pointer::p; integer,intent(in)::i,j
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%CO2
  end select
end subroutine

subroutine sc_litter_nitrate_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%nitrate
  end select
end subroutine

subroutine sc_litter_ammonium_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%ammonium
  end select
end subroutine

subroutine sc_litter_CO2_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%CO2
  end select
end subroutine

subroutine sc_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%nitrif
  end select
end subroutine

subroutine sc_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%denitrif
  end select
end subroutine

subroutine sc_litter_nitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%nitrif
  end select
end subroutine

subroutine sc_litter_denitrif_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(i)%denitrif
  end select
end subroutine


subroutine sc_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%C_in(j)
  end select
end subroutine

subroutine sc_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%N_in(j)
  end select
end subroutine

subroutine sc_litter_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%C_in(i)
  end select
end subroutine

subroutine sc_litter_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%N_in(i)
  end select
end subroutine

subroutine sc_litter_dissolved_carbon_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%dissolved_carbon(i)
  end select
end subroutine

subroutine sc_litter_dissolved_nitrogen_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%dissolved_nitrogen(i)
  end select
end subroutine

subroutine sc_litter_livingMicrobeC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%livingMicrobeC
  end select
end subroutine

subroutine sc_litter_livingMicrobeN_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%litterCohorts(i)%livingMicrobeN
  end select
end subroutine

subroutine sc_protected_C_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_C_in(j)
  end select
end subroutine

subroutine sc_protected_N_in_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_N_in(j)
  end select
end subroutine

subroutine sc_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%C_turnover(j)
  end select
end subroutine

subroutine sc_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%N_turnover(j)
  end select
end subroutine

subroutine sc_protected_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_C_turnover(j)
  end select
end subroutine

subroutine sc_protected_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%protected_N_turnover(j)
  end select
end subroutine

subroutine sc_litter_C_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%C_turnover(i)
  end select
end subroutine

subroutine sc_litter_N_turnover_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j; real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(j)%N_turnover(i)
  end select
end subroutine

subroutine sc_protected_C_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%protectedC(k)
  end select
end subroutine

subroutine sc_protected_N_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j,k;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%litterCohorts(j)%protectedN(k)
  end select
end subroutine

subroutine sc_litter_litterC_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%litterC(j)
  end select
end subroutine

subroutine sc_litter_litterN_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%litterN(j)
  end select
end subroutine

subroutine sc_litter_protectedC_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%protectedC(j)
  end select
end subroutine

subroutine sc_litter_protectedN_ptr(t,i,j,k,p)
  type(land_tile_type),pointer::t
  integer,intent(in)::i,j,k
  real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%litter_corpse(k)%litterCohorts(i)%protectedN(j)
  end select
end subroutine

subroutine sc_DOC_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%dissolved_carbon(j)
  end select
end subroutine

subroutine sc_DON_ptr(t,i,j,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i,j;real,pointer::p
  p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
     p=>s%org_matter(i)%dissolved_nitrogen(j)
  end select
end subroutine

subroutine sc_negative_litter_C_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
    p=>s%neg_litt_C(i)
  end select
end subroutine

subroutine sc_negative_litter_N_ptr(t,i,p)
  type(land_tile_type),pointer::t; integer,intent(in)::i; real,pointer::p; p=>NULL()
  if(.not.associated(t))       return
  if(.not.associated(t%soilc)) return
  select type(s=>t%soilc); class is (soil_BGC_CORPSE_t)
    p=>s%neg_litt_N(i)
  end select
end subroutine


end module