!     Este programa lee la salida del la nube en 3-D y calcula
!     la perturbacion de aerosoles y vapor en la estela

program aeroana
   USE cant01
   implicit none

   include 'dimen.i'
   include 'cant01.i'
   include 'permic.i'
   include 'perdim.i'
   include 'aeroana.i'


   character*2 paso
   character*2 nombre
   character*36 arch

   integer imax(-3:nz1+3),jmax(-3:nz1+3)
   real aermax(-3:nz1+3)
   real aerm(-3:nz1+3)

   integer lmax(-3:nz1+3),mmax(-3:nz1+3)
   real vapmax(-3:nz1+3)
   real vapm(-3:nz1+3)

!**   datos generales
   nombre='31'
   paso='10'

!*    lectura de los datos generales y de base


!*        lectura de la nube

   arch='outputdata/nube'//nombre//paso//'.sal'

   write(*,*) arch

   open(unit=60,file=arch,status='unknown',form='unformatted')
   read(60) U2,V2,W2,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,Qgra1,aer1
   close(60)

!*        calculo y grabacion

   open(unit=15,file='outputdata/aerest'//nombre//'.'//paso)
   open(unit=16,file='outputdata/vapest'//nombre//'.'//paso)

   do 100 k=0,nz1
      aerm(k)=0.
      imax(k)=0
      jmax(k)=0
      aermax(k)=0.
      vapm(k)=0.
      lmax(k)=0
      mmax(k)=0
      vapmax(k)=0.

      do 110 i=1,nx1
         do 110 j=1,nx1
            if (aer1(i,j,k).gt.aermax(k)) then
               aermax(k)=aer1(i,j,k)
               imax(k)=i
               jmax(k)=j
            endif
            if (Qvap1(i,j,k).gt.vapmax(k)) then
               vapmax(k)=Qvap1(i,j,k)
               lmax(k)=i
               mmax(k)=j
            endif


110   continue

      do 120 i=imax(k)-2,imax(k)+2
         do 120 j=jmax(k)-2,jmax(k)+2
            aerm(k)=aerm(k)+aer1(i,j,k)
            vapm(k)=vapm(k)+aer1(i,j,k)
120   continue

      aerm(k)=aerm(k)/5**2.
      vapm(k)=vapm(k)/5**2.

      write(15,*) k,aerm(k),imax(k),jmax(k),aermax(k)
      write(16,*) k,vapm(k),lmax(k),mmax(k),vapmax(k)

100 continue

   close(15)
   close(16)

end
