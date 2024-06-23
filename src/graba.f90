subroutine graba120(Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV,&
   u_original, u_perturbed, v_original, v_perturbed, w_original, w_perturbed,&
   thermal_property_1, thermal_property_2, pressure_original, pressure_perturbed,&
   vapor_amt, perturbed_vapor_amt, drop_amt, perturbed_drop_amt, rain_amt,&
   perturbed_rain_amt, crystal_amt, perturbed_crystal_amt, snow_amt, perturbed_snow_amt,&
   hail_amt, perturbed_hail_amt, spray_amt, perturbed_spray_amt, heat_force, Tvis, Tlvl,&
   Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, Vtgra0, Qvaprel, aerrel, Eautcn, Eacrcn)


   use dimensions
   use config
   implicit none
   real, dimension(-3:nz1+3), intent(in) :: Den0, Temp0, Tita0, Pres00, Qvap0,&
      cc2, aer0, UU, VV
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_original,&
      u_perturbed, v_original, v_perturbed, w_original, w_perturbed,&
      thermal_property_1, thermal_property_2, pressure_original,&
      pressure_perturbed, heat_force, vapor_amt, perturbed_vapor_amt, drop_amt,&
      perturbed_drop_amt, spray_amt, perturbed_spray_amt, rain_amt,&
      perturbed_rain_amt, crystal_amt, perturbed_crystal_amt, snow_amt,&
      perturbed_snow_amt, hail_amt, perturbed_hail_amt

   real, dimension(210:320), intent(in) :: Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs,&
      Eautcn, Eacrcn

   real, dimension(-3:2*nz1+5), intent(in) :: Av,Vtnie,Vtgra0
   real, dimension(nz1), intent(in)  :: Qvaprel,aerrel
   integer :: unit_number

   open(newunit = unit_number, file = output_directory//"inis.da", status =&
      'unknown', form = 'unformatted')
   write(unit_number) Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV
   close(unit_number)

   open(newunit = unit_number, file = output_directory//"velos.da", status =&
      'unknown', form = 'unformatted')
   rewind unit_number
   write(unit_number) u_original, u_perturbed, v_original, v_perturbed,&
      w_original, w_perturbed, thermal_property_1, thermal_property_2,&
      pressure_original, pressure_perturbed, vapor_amt, perturbed_vapor_amt,&
      drop_amt, perturbed_drop_amt, rain_amt, perturbed_rain_amt, crystal_amt,&
      perturbed_crystal_amt, snow_amt, perturbed_snow_amt, hail_amt,&
      perturbed_hail_amt, spray_amt, perturbed_spray_amt, heat_force
   close(unit_number)

   open(newunit = unit_number, file = output_directory//"varconz.da", status = &
      'unknown', form = 'unformatted')
   rewind unit_number
   write(unit_number)  Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, Vtgra0,&
      Qvaprel, aerrel, Eautcn, Eacrcn
   close(unit_number)

end subroutine graba120

subroutine graba231(k, w_perturbed, thermal_property_1, vapor_amt, rain_amt,&
   hail_amt, spray_amt, Qvap0, aer0, file_number)
   !### Grabacion 2D ###
   use dimensions
   use config
   implicit none
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: w_perturbed,&
      thermal_property_1, vapor_amt, rain_amt, hail_amt, spray_amt
   real, dimension(-3:nz1+3), intent(in) :: Qvap0, aer0
   integer, intent(in) :: k
   character(len=3), intent(in) :: file_number
   character(len=100) :: nombre
   integer :: i,j
   integer :: unit_number

   nombre = output_directory//'w_perturbed'//file_number//'.m'
   open(newunit = unit_number, file = nombre)
   do i=1,nx1
      write(unit_number,2000) (w_perturbed(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Tita'//file_number//'.m'
   open(newunit = unit_number, file = nombre)
   do i=1,nx1
      write(unit_number,2000) (thermal_property_1(i,j,0),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qvap'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) ((vapor_amt(i,j,k)+Qvap0(k)),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qllu'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (rain_amt(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Aero'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (spray_amt(i,j,0)+aer0(0),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qgra'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (hail_amt(i,j,1),j=1,nx1)
   end do
   close(unit_number)


2000 format(50E11.3)

end subroutine graba231

subroutine graba320(u_original, v_original, w_original, thermal_property_1,&
   pressure_original, vapor_amt, drop_amt, rain_amt, crystal_amt, snow_amt,&
   hail_amt, spray_amt, file_number)
   !### Grabacion 3D ###
   use dimensions
   use config
   implicit none

   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_original,&
      v_original, w_original, thermal_property_1, pressure_original, vapor_amt,&
      drop_amt, rain_amt, crystal_amt, snow_amt, hail_amt, spray_amt
   character(len=3), intent(in) :: file_number
   character(len=30) :: file_name
   integer :: unit_number

   file_name = output_directory//"nube"//trim(file_number)//'.sal'
   open(newunit = unit_number, file = file_name, status = 'unknown', form =&
      'unformatted')

   write(unit_number) u_original, v_original, w_original, thermal_property_1,&
      pressure_original, vapor_amt, drop_amt, rain_amt, crystal_amt, snow_amt,&
      hail_amt, spray_amt
   close(unit_number)

end subroutine graba320
