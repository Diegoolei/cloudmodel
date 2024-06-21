subroutine graba120(Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV,&
   u_original,u_perturbed,v_original,v_perturbed,w_original,w_perturbed,&
   thermal_property_1,thermal_property_2,pressure_original,pressure_perturbed,&
   Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,aer1,&
   aer2,heat_force,Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,&
   Eautcn,Eacrcn)

   use dimensions
   use config
   implicit none
   real, dimension(-3:nz1+3), intent(in) :: Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_original,&
      u_perturbed, v_original, v_perturbed, w_original, w_perturbed,&
      thermal_property_1, thermal_property_2, pressure_original,&
      pressure_perturbed, heat_force, Qvap1, Qvap2, Qgot1, Qgot2, aer1, aer2,&
      Qllu1, Qllu2, Qcri1, Qcri2, Qnie1, Qnie2, Qgra1, Qgra2

   real, dimension(210:320), intent(in) :: Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Eautcn, Eacrcn

   real, dimension(-3:2*nz1+5), intent(in) :: Av,Vtnie,Vtgra0
   real, dimension(nz1), intent(in)  :: Qvaprel,aerrel
   integer :: unit_number

   open(newunit=unit_number,file=output_directory//"inis.da",status='unknown',form='unformatted')
   write(unit_number) Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV
   close(unit_number)

   open(newunit=unit_number,file=output_directory//"velos.da",status='unknown',form='unformatted')
   rewind unit_number
   write(unit_number) u_original, u_perturbed, v_original, v_perturbed, w_original,&
      w_perturbed, thermal_property_1, thermal_property_2, pressure_original,&
      pressure_perturbed, Qvap1, Qvap2, Qgot1, Qgot2, Qllu1, Qllu2, Qcri1, Qcri2,&
      Qnie1, Qnie2, Qgra1, Qgra2, aer1, aer2, heat_force
   close(unit_number)

   open(newunit=unit_number,file=output_directory//"varconz.da",status='unknown',form='unformatted')
   rewind unit_number
   write(unit_number)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn
   close(unit_number)

end subroutine graba120

subroutine graba231(k, w_perturbed, thermal_property_1, Qvap1, Qllu1, Qgra1, aer1, Qvap0, aer0, file_number)
   !### Grabacion 2D ###
   use dimensions
   use config
   implicit none
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: w_perturbed, thermal_property_1, Qvap1, Qllu1, Qgra1, aer1
   real, dimension(-3:nz1+3), intent(in) :: Qvap0, aer0
   integer, intent(in) :: k
   character(len=3), intent(in) :: file_number
   character(len=100) :: nombre
   integer :: i,j
   integer :: unit_number

   nombre=output_directory//'w_perturbed'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (w_perturbed(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre=output_directory//'Tita'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (thermal_property_1(i,j,0),j=1,nx1)
   end do
   close(unit_number)

   nombre=output_directory//'Qvap'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) ((Qvap1(i,j,k)+Qvap0(k)),j=1,nx1)
   end do
   close(unit_number)

   nombre=output_directory//'Qllu'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (Qllu1(i,j,1),j=1,nx1)
   end do
   close(unit_number)

   nombre=output_directory//'Aero'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (aer1(i,j,0)+aer0(0),j=1,nx1)
   end do
   close(unit_number)

   nombre=output_directory//'Qgra'//file_number//'.m'
   open(newunit=unit_number,file=nombre)
   do i=1,nx1
      write(unit_number,2000) (Qgra1(i,j,1),j=1,nx1)
   end do
   close(unit_number)


2000 format(50E11.3)

end subroutine graba231

subroutine graba320(u_original, v_original, w_original, thermal_property_1,&
   pressure_original, Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1, file_number)
   !### Grabacion 3D ###
   use dimensions
   use config
   implicit none

   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: u_original,&
      v_original, w_original, thermal_property_1, pressure_original, Qvap1,&
      Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1
   character(len=3), intent(in) :: file_number
   character(len=30) :: file_name
   integer :: unit_number

   file_name = output_directory//"nube"//trim(file_number)//'.sal'
   open(newunit=unit_number,file= file_name,status='unknown',form='unformatted')

   write(unit_number) u_original, v_original, w_original, thermal_property_1,&
      pressure_original, Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1
   close(unit_number)

end subroutine graba320
