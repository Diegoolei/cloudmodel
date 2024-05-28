subroutine graba120(Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV,&
   U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
   Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,aer1,aer2,Fcalo,&
   Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn,output_directory)
   
   USE dimen
   implicit none
   character(len=80), intent(in) :: output_directory
   real, dimension(-3:nz1+3), intent(in) :: Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: U1, U2, V1, V2, W1, W2, Titaa1,Titaa2
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: Pres1,Pres2,Fcalo
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: Qvap1,Qvap2,Qgot1,Qgot2,aer1,aer2,Qllu1,Qllu2,Qcri1,Qcri2,&
      Qnie1,Qnie2,Qgra1,Qgra2
   
   real, dimension(210:320), intent(in) :: Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Eautcn, Eacrcn
   
   real, dimension(-3:2*nz1+5), intent(in) :: Av,Vtnie,Vtgra0
   real, dimension(nz1), intent(in)  :: Qvaprel,aerrel

   open(unit=40,file=output_directory//"inis.da",status='unknown',form='unformatted')
   open(unit=41,file=output_directory//"velos.da",status='unknown',form='unformatted')
   rewind 41
   open(unit=42,file=output_directory//"varconz.da",status='unknown',form='unformatted')
   rewind 42

   write(40) Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV
!$$
   write(41) U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
      Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,aer1,aer2,Fcalo
   write(42)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn

   close(42)
   close(41)
   close(40)
end subroutine graba120

subroutine graba231(k, t1, W2, Titaa1, Qvap1, Qllu1, Qgra1, aer1, Qvap0, aer0, file_number,output_directory)
   USE dimen
   implicit none
   character(len=80), intent(in) :: output_directory
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: W2, Titaa1, Qvap1, Qllu1, Qgra1, aer1
   real, dimension(-3:nz1+3), intent(in) :: Qvap0, aer0
   integer, intent(in) :: k, t1
   character(len=3), intent(in) :: file_number
   character(len=100) :: nombre
   integer :: i,j
   nombre=output_directory//'W2'//file_number//'.m'
   open(unit=30,file=nombre)
   do 360 i=1,nx1
      write(30,2000) (W2(i,j,1),j=1,nx1)
360 continue
   close(30)
   nombre=output_directory//'Tita'//file_number//'.m'
   open(unit=30,file=nombre)
   do 370 i=1,nx1
      write(30,2000) (Titaa1(i,j,0),j=1,nx1)
370 continue
   close(30)
   nombre=output_directory//'Qvap'//file_number//'.m'
   open(unit=30,file=nombre)
   do 380 i=1,nx1
      write(30,2000) ((Qvap1(i,j,k)+Qvap0(k)),j=1,nx1)
380 continue
   close(30)
   nombre=output_directory//'Qllu'//file_number//'.m'
   open(unit=30,file=nombre)
   do 395 i=1,nx1
      write(30,2000) (Qllu1(i,j,1),j=1,nx1)
395 continue
   close(30)
   nombre=output_directory//'Aero'//file_number//'.m'
   open(unit=30,file=nombre)
   do 385 i=1,nx1
      write(30,2000) (aer1(i,j,0)+aer0(0),j=1,nx1)
385 continue
   close(30)

   nombre=output_directory//'Qgra'//file_number//'.m'
   open(unit=30,file=nombre)
   do 525 i=1,nx1
      write(30,2000) (Qgra1(i,j,1),j=1,nx1)
525 continue
   close(30)


2000 format(50E11.3)

end subroutine graba231

subroutine graba320(U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1,file_number)
   USE dimen
   USE config
   implicit none

   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1
   character(len=3), intent(in) :: file_number
   character(len=:), allocatable :: file_name

   file_name = output_directory//"nube"//trim(file_number)//'.sal'
   open(unit=60,file= file_name,status='unknown',form='unformatted')
   ! U1 - Pres1 : perdim.i
   ! Qvap1 - aer1 : permic.i
   write(60) U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1
   close(60)

end subroutine graba320
