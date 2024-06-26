subroutine graba120(Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV,&
   u_perturbed_base, u_perturbed_new, v_perturbed_base, v_perturbed_new, w_perturbed_base, w_perturbed_new,&
   potential_temperature_base, potential_temperature_new, pressure_base, pressure_new,&
   vapor_base, vapor_new, drop_base, drop_new, rain_base,&
   rain_new, crystal_base, crystal_new, snow_base, snow_new,&
   hail_base, hail_new, aerosol_base, aerosol_new, heat_force, Tvis, Tlvl,&
   Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, Vtgra0, Qvaprel, aerrel, Eautcn, Eacrcn)


   use dimensions
   use config
   implicit none
   real, dimension(-3:nz1+3), intent(in) :: Den0, Temp0, Tita0, Pres00, Qvap0,&
      cc2, aer0, UU, VV
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_perturbed_base,&
      u_perturbed_new, v_perturbed_base, v_perturbed_new, w_perturbed_base, w_perturbed_new,&
      potential_temperature_base, potential_temperature_new, pressure_base,&
      pressure_new, heat_force, vapor_base, vapor_new, drop_base,&
      drop_new, aerosol_base, aerosol_new, rain_base,&
      rain_new, crystal_base, crystal_new, snow_base,&
      snow_new, hail_base, hail_new

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
   write(unit_number) u_perturbed_base, u_perturbed_new, v_perturbed_base, v_perturbed_new,&
      w_perturbed_base, w_perturbed_new, potential_temperature_base, potential_temperature_new,&
      pressure_base, pressure_new, vapor_base, vapor_new,&
      drop_base, drop_new, rain_base, rain_new, crystal_base,&
      crystal_new, snow_base, snow_new, hail_base,&
      hail_new, aerosol_base, aerosol_new, heat_force
   close(unit_number)

   open(newunit = unit_number, file = output_directory//"varconz.da", status = &
      'unknown', form = 'unformatted')
   rewind unit_number
   write(unit_number)  Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, Vtgra0,&
      Qvaprel, aerrel, Eautcn, Eacrcn
   close(unit_number)

end subroutine graba120

subroutine graba231(k, w_perturbed_new, potential_temperature_base, vapor_base, rain_base,&
   hail_base, aerosol_base, Qvap0, aer0, file_number)
   !### Grabacion 2D ###
   use dimensions
   use config
   implicit none
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: w_perturbed_new,&
      potential_temperature_base, vapor_base, rain_base, hail_base, aerosol_base
   real, dimension(-3:nz1+3), intent(in) :: Qvap0, aer0
   integer, intent(in) :: k
   character(len=3), intent(in) :: file_number
   character(len=100) :: nombre
   integer :: i,j
   integer :: unit_number

   nombre = output_directory//'w_perturbed_new'//file_number//'.m'
   open(newunit = unit_number, file = nombre)
   do i=1,nx1
      write(unit_number,2000) (w_perturbed_new(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Tita'//file_number//'.m'
   open(newunit = unit_number, file = nombre)
   do i=1,nx1
      write(unit_number,2000) (potential_temperature_base(i,j,0),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qvap'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) ((vapor_base(i,j,k)+Qvap0(k)),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qllu'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (rain_base(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Aero'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (aerosol_base(i,j,0)+aer0(0),j=1,nx1)
   end do
   close(unit_number)

   nombre = output_directory//'Qgra'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (hail_base(i,j,1),j=1,nx1)
   end do
   close(unit_number)


2000 format(50E11.3)

end subroutine graba231

subroutine graba320(u_perturbed_base, v_perturbed_base, w_perturbed_base, potential_temperature_base,&
   pressure_base, vapor_base, drop_base, rain_base, crystal_base, snow_base,&
   hail_base, aerosol_base, file_number)
   !### Grabacion 3D ###
   use dimensions
   use config
   implicit none

   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_perturbed_base,&
      v_perturbed_base, w_perturbed_base, potential_temperature_base, pressure_base, vapor_base,&
      drop_base, rain_base, crystal_base, snow_base, hail_base, aerosol_base
   character(len=3), intent(in) :: file_number
   character(len=30) :: file_name
   integer :: unit_number

   file_name = output_directory//"nube"//trim(file_number)//'.sal'
   open(newunit = unit_number, file = file_name, status = 'unknown', form =&
      'unformatted')

   write(unit_number) u_perturbed_base, v_perturbed_base, w_perturbed_base, potential_temperature_base,&
      pressure_base, vapor_base, drop_base, rain_base, crystal_base, snow_base,&
      hail_base, aerosol_base
   close(unit_number)

end subroutine graba320
