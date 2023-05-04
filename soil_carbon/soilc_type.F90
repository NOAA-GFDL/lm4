module soilc_type_mod

implicit none; private

public :: soilc_t

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
end interface

end module
