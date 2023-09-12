!     Este programa lee la salida de aerdif y realiza sumas
!     en Y para que el planos X-Z puedan ser levantados por el Matlab

      program coraer
         implicit none

         include 'dimen.i'
         include 'const.i'
         include 'perdim.i'
         include 'permic.i'
         include 'estbas.i'


         integer ii,i,j,k,l,m,n
         integer pasoaux(25),lvar,nvar
         real*8 aux
         character*100 paso
         character*2 nombre
         character*25 arch,archaux
         character*2 var(13)
         integer*4 variable(nx1,nz1)
         real totpos,totneg

!**   datos generales
         nombre='25'

         paso='01020304050607080910111213141516171819202122232425'
     &        //'26272829303132333435363738394041424344454647484950'
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
         pasoaux(15)=1
         pasoaux(16)=0
         pasoaux(17)=1
         pasoaux(18)=0
         pasoaux(19)=1
         pasoaux(20)=0
         pasoaux(21)=1
         pasoaux(22)=0
         pasoaux(23)=1
         pasoaux(24)=0
         pasoaux(25)=1


         write(*,*) 'aqui'

         do 10 ii=1,25
            if(pasoaux(ii).eq.1) then

!*        lectura de la nube

               arch='aerdif'//paso(ii*2-1:ii*2)

               write(*,*) arch

               open(unit=60,file=arch,status='unknown',form='unformatted')
               read(60) aer1
               close(60)

!*        cuenta de positivos y negativos

               totpos=0
               totneg=0
               do 300 i=1,nx1
                  do 300 k=1,nz1
                     aux=0
                     do 310 j=1,nx1

                        aux=aux+aer1(i,j,k)

                        if (aer1(i,j,k).gt.0) totpos=totpos+aer1(i,j,k)
                        if (aer1(i,j,k).lt.0) totneg=totneg+aer1(i,j,k)
  310                continue
                     variable(i,k)=aux*10.
  300          continue


               archaux='cortes\paer'//paso(ii*2-1:ii*2)

               write(*,*) archaux

               open(unit=10,file=archaux)
               do 90 k=1,nz1
                  write(10,1000) (variable(i,k),i=1,nx1)
   90          continue
               close(10)
   80          continue

   50          continue

               write(*,*) ii,totpos,totneg,totpos+totneg

            endif
   10    continue

 1000    format(50i7)
      end
