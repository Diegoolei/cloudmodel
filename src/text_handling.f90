module io
   !! I/O related procedures
   implicit none

   interface str
      module procedure :: str_gen
   end interface

contains
   function str_gen(int_in) result(str_out)
      class(*), intent(in) :: int_in
      character(len=99) :: str_mid
      character(len=:), allocatable :: str_out

      select type(int_in)
       type is (real)
         write (str_mid, *) int_in
       type is (integer)
         write (str_mid, *) int_in
      end select
      str_out = trim(adjustl(str_mid))
   end function str_gen
end module io
