*     Este programa lee la salida de aerdif y realiza cortes
*     planos para que puedan ser levantados por el Matlab

      program coraer
      implicit none

      include 'dimen.i'
      include 'const.i'
      include 'perdim.i'
      include 'permic.i'
      include 'estbas.i'


      integer ii,i,j,k,l,m,n
      integer pasoaux(25),ycort(50),lvar,nvar
      real*8 aux
      character*100 paso
      character*2 nombre
      character*25 arch,archaux
      character*2 var(13)
      integer*4 variable(nx1,nz1)
      real totpos,totneg

***   datos generales
      nombre='25'

      paso='01020304050607080910111213141516171819202122232425'
     &     //'26272829303132333435363738394041424344454647484950'
      pasoaux(1)=0
      pasoaux(2)=0
      pasoaux(3)=0
      pasoaux(4)=0
      pasoaux(5)=0
      pasoaux(6)=0
      pasoaux(7)=0
      pasoaux(8)=0
      pasoaux(9)=1
      pasoaux(10)=0
      pasoaux(11)=1
      pasoaux(12)=0
      pasoaux(13)=1
      pasoaux(14)=0
      pasoaux(15)=0
      pasoaux(16)=0
      pasoaux(17)=0
      pasoaux(18)=0
      pasoaux(19)=0
      pasoaux(20)=0
      pasoaux(21)=0
      pasoaux(22)=0
      pasoaux(23)=0
      pasoaux(24)=0
      pasoaux(25)=0

      ycort(1)=0
      ycort(2)=0
      ycort(3)=0
      ycort(4)=0
      ycort(5)=0
      ycort(6)=0
      ycort(7)=0
      ycort(8)=0
      ycort(9)=0
      ycort(10)=0
      ycort(11)=0
      ycort(12)=0
      ycort(13)=0
      ycort(14)=0
      ycort(15)=0
      ycort(16)=0
      ycort(17)=0
      ycort(18)=0
      ycort(19)=0
      ycort(20)=0
      ycort(21)=0
      ycort(22)=0
      ycort(23)=0
      ycort(24)=0
      ycort(25)=0
      ycort(26)=1
      ycort(27)=0
      ycort(28)=0
      ycort(29)=0
      ycort(30)=0
      ycort(31)=0
      ycort(32)=0
      ycort(33)=0
      ycort(34)=0
      ycort(35)=0
      ycort(36)=0
      ycort(37)=0
      ycort(38)=0
      ycort(39)=0
      ycort(40)=0
      ycort(41)=0
      ycort(42)=0
      ycort(43)=0
      ycort(44)=0
      ycort(45)=0
      ycort(46)=0
      ycort(47)=0
      ycort(48)=0
      ycort(49)=0
      ycort(50)=0

      write(*,*) 'aqui' 

      do 10 ii=1,25
        if(pasoaux(ii).eq.1) then

**        lectura de la nube

          arch='aerdif'//paso(ii*2-1:ii*2)

          write(*,*) arch

          open(unit=6,file=arch,status='unknown',form='unformatted')
          read(6) aer1
          close(6)

**        cuenta de positivos y negativos

          totpos=0
          totneg=0
          do 300 i=1,nx1           
          do 300 j=1,nx1           
          do 300 k=1,nz1           
	    if (aer1(i,j,k).gt.0) totpos=totpos+aer1(i,j,k)
	    if (aer1(i,j,k).lt.0) totneg=totneg+aer1(i,j,k)
  300     continue

**       cortes (no incluyen el piso)
          do 50 j=1,nx1

            if (ycort(j).eq.1) then
                write(*,*) j
                do 60 k=1,nz1
                  do 70 i=1,nx1
                    variable(i,k)=aer1(i,j,k)*10.
   70             continue
   60         continue

              archaux='cortes\caer'//paso(ii*2-1:ii*2)

              write(*,*) archaux

              open(10,file=archaux)
                do 90 k=1,nz1
                  write(10,1000) (variable(i,k),i=1,nx1)
 90            continue
                close(10)
 80          continue

            endif
   50     continue

        write(*,*) ii,totpos,totneg,totpos+totneg

        endif
   10 continue

 1000 format(50i7)
      end
