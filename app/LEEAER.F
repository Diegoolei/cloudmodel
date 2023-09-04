!***** PROGRAMA DE LECTURA
!      Lee aerdif y nuestra una seleccion de puntos en la pantalla

      program leeaer

         implicit none

!     Parametros del modelo

         integer nx1,ny1,nz1
         integer x1,x2,y1,y2,z1,z2
         integer i,j,k
         parameter (nx1=50,ny1=50,nz1=45)

         real aerdif(-3:nx1+3,-3:nx1+3,-3:nz1+3)
         integer aer1(-3:nx1+3,-3:nx1+3,-3:nz1+3)
         character*2 nombre

         write(*,*) 'Ingrese nombre a leer '
         read(*,*) nombre

         write(*,*) 'Ingrese intervalo en x a mostrar '
         read(*,*) x1,x2
         write(*,*) 'Ingrese intervalo en y a mostrar '
         read(*,*) y1,y2
         write(*,*) 'Ingrese intervalo en z a mostrar '
         read(*,*) z1,z2

         open(unit=15, file='saleaer')

         open(unit=10,file='aerdif'//nombre,
     &        status='unknown',form='unformatted')
         read(10) aerdif
         close(10)

         do 100 i=1,nx1
            do 100 j=1,nx1
               do 100 k=1,nz1
                  aer1(i,j,k)=aerdif(i,j,k)*1000
  100    continue


         do 10 j=y1,y2
            write(*,*) j
            write(15,*) j
            do 20 k=z1,z2
               write(*,*) k,(aer1(i,j,k),i=x1,x2)
               write(15,*) k,(aer1(i,j,k),i=x1,x2)
   20       continue
   10    continue

         close(15)

      end








