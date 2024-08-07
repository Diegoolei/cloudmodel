****** PROGRAMA DE LECTURA y escritura para ser leido por foo_to_v5d
*      Este programa lee los archivos generados por el modelo de
*      Hielo10 (model01)
*      no tomamos en cuenta a la presion
*      sigue la nube solamente
*      calculamos el perfil visible de la nube (Qtot>=1e-3 g/m^3)

       program leelluvi

       implicit none

*      Parametros del modelo 

       integer nx1,ny1,nz1,nvar
       parameter (nx1=50,ny1=50,nz1=45,nvar=12)

      real dt1

      real U1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real V1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real W1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Titaa1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Pres1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qvap1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qgot1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qllu1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qcri1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qnie1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real Qgra1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      real aer1(-3:nx1+3,-3:nx1+3,-2:nz1+2)

      real vars(nx1,ny1,nz1,nvar)

      real Temp0(-3:nz1+3)
      real Tita0(-3:nz1+3)
      real Presi0(-3:nz1+3)
      real Pres00(-3:nz1+3)
      real Den0(-3:nz1+3)      
      real UU(-3:nz1+3)
      real VV(-3:nz1+3)
      real aer0(-3:nz1+3)
      real Qvap0(-3:nz1+3)

      character*12 httfla
      character*6 hthead

      character*6 httima
      character*15 entrav5d
      character*400 bre
      character*57 bre1,bre2,bre3,bre4
      character* 3 ent

       integer lx,ly,lz,ls,paso,inpasos,fipasos,delpaso,iaux
       integer haux,maux,saux,hora
       integer i,j,k,tt,nt1,nt
       real aux,T

      hthead='nube11'
      bre='01020304050607080910111213141516171819202122232425'

      ent='ent'
      
***********************************************************

       write(*,*) '     PROGRAMA DE LECTURA DE HIELO1 '
       write(*,*) '          Y ESCRITURA PARA VIS5D'
       write(*,*)
       write(*,*)
       write(*,*)
       write (*,*) 'Ingrese el paso de tiempo inicial '
       read(*,*) inpasos

       write (*,*) 'Ingrese el paso de tiempo final '
       read(*,*) fipasos

       write(*,*) 'Ingrese delta de tiempo (en segundos)  '
       read(*,*) delpaso

       write(*,*) 'Ingrese dt1 (en segundos)  '
       read(*,*) dt1

       nt1=delpaso/dt1
       nt=nt1*(fipasos-inpasos+1)

*********************************************************************

*     Lectura del estado base
      open (7,file='inic03.sa')
      do 100 k=0,nz1

	read(7,210) i,Temp0(k),Tita0(k),Presi0(k),Pres00(k),Den0(k),
     &               aer0(k),Qvap0(k),UU(k),VV(k)

  100 continue
      close(7)
  210 format(I3,9E12.4)

*     inicializacion todas las variables en su estado base
      do 105 i=1,nx1
      do 105 j=1,ny1
      do 105 k=1,nz1
       vars(i,j,k,1)=UU(k)
       vars(i,j,k,2)=VV(k)
       vars(i,j,k,3)=0.
       vars(i,j,k,4)=Temp0(k)
       vars(i,j,k,5)=Qvap0(k)*1000.
       vars(i,j,k,6)=0.
       vars(i,j,k,7)=0.
       vars(i,j,k,8)=0.
       vars(i,j,k,9)=0.
       vars(i,j,k,10)=0.
       vars(i,j,k,11)=0.
       vars(i,j,k,12)=0.
  105 continue
     
      write(*,*)
      write(*,*)

       write(*,*) '   Comienza la escritura de los datos para V5D'

      entrav5d='entra30'

      open(unit=3,file=entrav5d)

      write(3,*) fipasos+1-inpasos
      write(3,*) nvar
      write(3,*) ny1
      write(3,*) nx1-10
      write(3,*) nz1

      write(3,111) 'U  '
      write(3,111) 'V  '
      write(3,111) 'W  '
      write(3,111) 'T '
      write(3,111) 'Qv '
      write(3,111) 'Qg '
      write(3,111) 'Ql '
      write(3,111) 'Qc '
      write(3,111) 'Qn '
      write(3,111) 'Qh '
      write(3,111) 'Ae '
      write(3,111) 'Con'

      do 1 paso=inpasos,fipasos
        haux=paso*delpaso/3600
        maux=paso*delpaso/60-haux*60
        saux=paso*delpaso-maux*60-haux*3600
        hora=(haux+12)*10000+maux*100+saux
        write(3,*) hora
        write(*,*) hora,paso,delpaso,haux,maux,saux
    1 continue

      do 2 paso=inpasos,fipasos
        write(3,113) '98100'
    2 continue

      write(3,*) '0.'
      write(3,*) '.01'
      write(3,*) '0.'
      write(3,*) '.01'
      write(3,*) '0.'
      write(3,*) '1.'
 111    format (A3)
 112    format (A5)
 113    format (A6)


       write(*,*)
       write(*,*)
       write(*,*)
       write (*,*) 'Comienzo de la lectura'
       write(*,*) 
       write(*,*) 



       write(*,*)
       write(*,*) 'Dado que la matriz G se ordena en primero filas (y)'
       write(*,*) ', luego columnas (x) invierto el orden entre x e y'
       write(*,*) 'Como (1,1) significa (N-O) la variable ly va desde '
       write(*,*) '       nya a 1 (sentido inverso)'


      do 3 paso=inpasos,fipasos

        httima=bre(paso*2-1:paso*2)//'.sal'
        httfla = hthead//httima

        write(*,*)
        write(*,*) 'leyendo el archivo: '
        write(*,*) httfla
*           pause

      open(unit=6,file=httfla,status='unknown',form='unformatted')
        read(6) U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1,
     &          Qgra1,aer1
      close(6)

      write(*,*) 'archivo leido OK'

      write(*,*) paso

      tt=paso*nt1

      do 10 ly=1,ny1
      do 10 lx=1,nx1
      do 10 lz=1,nz1

        i=lx
        j=ly

*      Escritura de las velocidades

       vars(i,j,lz,1)=UU(lz)+U1(lx,ly,lz)
       vars(i,j,lz,2)=VV(lz)+V1(lx,ly,lz)
       vars(i,j,lz,3)=W1(lx,ly,lz)

*      Escritura de la variable Temperatura
	  aux=Pres00(lz)+Pres1(i,j,lz)
	  T=(Tita0(lz)+Titaa1(i,j,lz))*aux
       vars(i,j,lz,4)=T

*      Escritura de la variable Qvap

       vars(i,j,lz,5)=(Qvap0(lz)+Qvap1(lx,ly,lz))*1000

*      Escritura de la variable Qgot

       vars(i,j,lz,6)=Qgot1(lx,ly,lz)*1000

*      Escritura de la variable Qllu

       vars(i,j,lz,7)=Qllu1(lx,ly,lz)*1000

*      Escritura de la variable Qcri

       vars(i,j,lz,8)=Qcri1(lx,ly,lz)*1000

*      Escritura de la variable Qnie

       vars(i,j,lz,9)=Qnie1(lx,ly,lz)*1000

*      Escritura de la variable Qgra

       vars(i,j,lz,10)=Qgra1(lx,ly,lz)*1000

**      Escritura de la variable Aer

       vars(i,j,lz,11)=aer1(lx,ly,lz)

**      Escritura de la variable Contorno visible
       aux=0.
       do 201 k=6,10
         aux=aux+vars(i,j,lz,k)
  201  continue
       if (aux.ge.1e-6) then
         vars(i,j,lz,12)=1.
       else
         vars(i,j,lz,12)=0.
       endif

   10 continue


*      Escritura de las variables
      write(*,*) 'Escritura de las variables'
      do 20 ls=1,nvar
         do 30 ly=1,ny1
         do 30 lx=6,nx1-5
         do 30 lz=1,nz1
           write(3,*) vars(lx,ly,lz,ls) 
   30   continue
   20 continue

    3 continue
     
       close(unit=3)


      write(*,*)
      write(*,*) 'Las variables han sido grabadas en el archivo :'
     *            , entrav5d
      write(*,*) 
      write(*,*)
      write(*,*) '     Fin del programa'



*       do 555 k=1,nz1
*         write(*,*) k,Qgot1(30,25,k),vars(30,25,k,12),vars(32,25,k,12)
*         pause
*  555  continue

       end








