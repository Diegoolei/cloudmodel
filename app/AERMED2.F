!     Este programa lee la salida del la nube en 3-D y calcula
!     la perturbacion promedio de los aerosoles alrededor del maximo

      program aeroana
      implicit none

      include 'dimen.i'
      include 'cant01.i'
      include 'const.i'
      include 'permic.i'
      include 'perdim.i'
      include 'aeroana.i'


      character*2 paso
      character*2 nombre
      character*25 arch

      integer imax(-3:nz1+3),jmax(-3:nz1+3)
      real aermax(-3:nz1+3)
      real aerm(-3:nz1+3)

!**   datos generales
      nombre='31'
      paso='15'

!*    lectura de los datos generales y de base


!*        lectura de la nube

          arch='nube'//nombre
     &          //paso//'.sal'

          write(*,*) arch

          open(unit=6,file=arch,status='unknown',form='unformatted')
          read(6) U2,V2,W2,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1
     &          ,Qgra1,aer1
          close(6)

!*        calculo y grabacion

          open(15,file='aermed.s'//paso)

          do 100 k=0,nz1
            aerm(k)=0.
            imax(k)=0
            jmax(k)=0
            aermax(k)=0.
  	    do 110 i=1,nx1
	    do 110 j=1,nx1
              if (aer1(i,j,k).gt.aermax(k)) then
                aermax(k)=aer1(i,j,k)
                imax(k)=i
                jmax(k)=j
              endif
  110       continue

            do 120 i=imax(k)-5,imax(k)+5
            do 120 j=jmax(k)-5,jmax(k)+5
              aerm(k)=aerm(k)+aer1(i,j,k)

  120       continue

            aerm(k)=aerm(k)/11**2.

            write(15,*) k,aerm(k),imax(k),jmax(k),aermax(k)

  100     continue

          close(15)

      end
