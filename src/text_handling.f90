module io
   !! I/O related procedures
   implicit none

   interface str
      module procedure :: str_gen
   end interface

contains
   function str_gen(int_in) result(str_out)
      integer, intent(in) :: int_in
      character(len=3) :: str_out_mid
      character(len=4) :: str_out
      str_out_mid = str_gen_aux(int_in)
      if (int_in < 10) then
         str_out = '0'//str_out_mid
      else
         str_out = str_out_mid
      end if
   end function str_gen

   function str_gen_aux(int_in) result(str_out)
      class(*), intent(in) :: int_in
      character(len=99) :: str_mid
      character(len=3) :: str_out

      select type (int_in)
      type is (real)
         write (str_mid, *) int_in
      type is (integer)
         write (str_mid, *) int_in
      end select
      str_out = trim(adjustl(str_mid))
   end function str_gen_aux
end module io
