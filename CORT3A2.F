*     Este programa lee la salida del la nube en 3-D y realiza cortes
*     planos para que puedan ser levantados por el Matlab

      program cort3a2
      implicit none

      include 'dimen.i'
      include 'const.i'
      include 'perdim.i'
      include 'permic.i'
      include 'estbas.i'


      integer ii,i,j,k,l,m,n
      integer pasoaux(25),ycort(50),lvar,nvar
      parameter (nvar=13)
      real*8 aux
      character*100 paso
      character*2 nombre
      character*25 arch,archaux
      character*2 var(13)
      integer*4 variable(nx1,nz1,nvar)


***   datos generales
      nombre='11'
      var(1)='UU'
      var(2)='VV'
      var(3)='WW'
      var(4)='TT'
      var(5)='GO'
      var(6)='LL'
      var(7)='CR'
      var(8)='NI'
      var(9)='GR'
      var(10)='VA'
      var(11)='AE'


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
      ycort(23)=1
      ycort(24)=1
      ycort(25)=1
      ycort(26)=1
      ycort(27)=1
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

**    lectura de los datos generales y de base

      open(5,file='inis.da')
      open(unit=3,file='varconz.da',status='unknown'
     &        ,form='unformatted')

      read(5,*) Den0,Temp0,Tita0,Pres00,Qvap0
     &                  ,cc2,aer0,UU,VV
      read(3)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,
     &             Qvaprel,aerrel,Eautcn,Eacrcn

      close(3)
      close(5)

      write(*,*) 'aqui' 

      do 10 ii=1,25
        if(pasoaux(ii).eq.1) then

**        lectura de la nube

          arch='d:\nube'//nombre//'\nube'//nombre
     &          //paso(ii*2-1:ii*2)//'.sal'

          write(*,*) arch

          open(unit=6,file=arch,status='unknown',form='unformatted')
          read(6) U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1
     &          ,Qgra1,aer1
          close(6)

        write(*,*) W1(nx1/2,nx1/2,15),W1(nx1/2,nx1/2,20)

**       cortes (no incluyen el piso)
          do 50 j=1,nx1

            if (ycort(j).eq.1) then
                write(*,*) j
                do 60 k=1,nz1
                  do 70 i=1,nx1
*                   calculo de la temperatua
                    aux=Pres00(k)+Pres1(i,j,k)
                    variable(i,k,4)=(Tita0(k)+Titaa1(i,j,k))*aux*1e2

*                   velocidades
                    variable(i,k,1)=(UU(k)+U1(i,j,k))*1e3
                    variable(i,k,2)=(VV(k)+V1(i,j,k))*1e3
                    variable(i,k,3)=W1(i,j,k)*1e3

*                   filtrado de gotitas
                      aux=0.
                      do 200 l=i-1,i+1
                      do 200 m=j-1,j+1
                      do 200 n=k-1,k+1
                        aux=aux+Qgot1(l,m,n)*1e6
  200                 continue
                      variable(i,k,5)=aux/27.
*                   filtrado de lluvia
                      aux=0.
                      do 210 l=i-1,i+1
                      do 210 m=j-1,j+1
                      do 210 n=k-1,k+1
                        aux=aux+Qllu1(l,m,n)*1e6
  210                 continue
                      variable(i,k,6)=aux/27.
*                   filtrado de cristales
                      aux=0.
                      do 220 l=i-1,i+1
                      do 220 m=j-1,j+1
                      do 220 n=k-1,k+1
                        aux=aux+Qcri1(l,m,n)*1e6
  220                 continue
                      variable(i,k,7)=aux/27.
*                   filtrado de nieve
                      aux=0.
                      do 230 l=i-1,i+1
                      do 230 m=j-1,j+1
                      do 230 n=k-1,k+1
                        aux=aux+Qnie1(l,m,n)*1e6
  230                 continue
                      variable(i,k,8)=aux/27.
*                   filtrado de granizos
                      aux=0.
                      do 240 l=i-1,i+1
                      do 240 m=j-1,j+1
                      do 240 n=k-1,k+1
                        aux=aux+Qgra1(l,m,n)*1e6
  240                 continue
                      variable(i,k,9)=aux/27.

                    variable(i,k,10)=(Qvap0(k)+Qvap1(i,j,k))*1e6
                    variable(i,k,11)=(aer0(k)+aer1(i,j,k))*1e1
                    variable(i,k,12)=Qvap1(i,j,k)*1e6
                    variable(i,k,13)=aer1(i,j,k)*1e1
   70             continue
   60         continue

              do 80 lvar=1,11
                archaux=var(lvar)//paso(j*2-1:j*2)//
     &                  'y'//paso(ii*2-1:ii*2)
                write(*,*) archaux
                open(10+lvar,file='cortes\'//archaux)
                do 90 k=1,nz1
                  write(10+lvar,1000) (variable(i,k,lvar),i=1,nx1)
 90            continue
                close(10+lvar)
 80          continue

            endif
   50     continue

        endif
   10 continue

 1000 format(50i7)
      end
