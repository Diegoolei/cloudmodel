subroutine graba320(U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1,t1,tie,bre)
   USE constants
   implicit none

   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2), intent(in) :: U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1
   integer, intent(in) :: t1
   character(len=2), intent(in) :: tie
   character(len=50), intent(in) :: bre
   character(len=80) :: file_name

   file_name = "outputdata/nube"//tie//bre(2*t1-1:2*t1)//'.sal'
   open(unit=60,file= file_name,status='unknown',form='unformatted')
   ! U1 - Pres1 : perdim.i
   ! Qvap1 - aer1 : permic.i
   write(60) U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1
   close(60)

end subroutine graba320
