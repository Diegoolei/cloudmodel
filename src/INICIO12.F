!23456789*123456789*123456789*123456789*123456789*123456789*123456789*
!     revisado  3/4/99
!     Condiciones iniciales para las variables dinamicas
!     Corresponde a una nube con hielo
!     Calcula la presion inicial en forma iterativa, primero supone aire
!      seco y luego integra considerando la densidad total, incluyendo el 
!      vapor.
!     Solo hay una perturbacion en temperatura para iniciar la conveccion
!     Incluye viento de corte, tipo frente


      subroutine condi

      implicit none
      include 'dimen.i'
      include 'permic.i'
      include 'perdim.i'
      include 'const.i'
      include 'estbas.i'
      include 'cant01.i'
      include 'inicio11.i'
   
!**** comienzo de la rutina

      call PP(G,Rd,dx1,nz,Presi0,P00)

!**   viento de corte
      vc=1.5/3000**2.
      vh=5000./pi

      do 5 k=0,nz1
        zeta=k*dx1
        if (zeta .le. 500.) then
          UU(k)=0.
	  VV(k)=0.

!
       write(*,*) '0',zeta,UU(k),VV(k)
!

        elseif (zeta .le. 2000.) then
  	  zeta1=zeta-500.
          aux=4.*(zeta1/1500.)**2.
          UU(k)=aux
          VV(k)=0.
!
       write(*,*) '1',zeta,UU(k),VV(k),vb

        elseif (zeta .le. 9000.) then
          zeta1=zeta-2000.          
          vb=zeta1/7000.
          UU(k)=4.-10.*vb**2.
          VV(k)=3.*vb**.5

       write(*,*) '3',zeta,UU(k),VV(k),vb
!
        else
          zeta1=zeta-9000.
          UU(k)=-6.+4.*(zeta1/9000.)**2.
          VV(k)=3.-5.*(zeta1/9000.)**.5
!
       write(*,*) '4',zeta,UU(k),VV(k),vb
!
        endif

       UU(k)=UU(k)*.7
       VV(k)=VV(k)*0.


!##
!      write(*,*) 'vel0',zeta,UU(k),VV(k),zeta1,aux,vb
!##

    5 continue


!**   calculo de 'constantes' que dependen de T

      open(3,file='ccc')

      do 400 k=313,210,-1
        
        Tk=k-T0

        Tvis(k)=4.9e-8*Tk+Vis0
	if (k.lt.273.15) Tvis(k)=Tvis(k)-1.2e-10*Tk**2.

!       calores latentes de evaporacion, fusion y sublimacion
        gam=.167+3.67e-4*k
        Tlvl(k)=Lvl0*(T0/k)**gam
        Tlsl(k)=(Lsl0+0.485*Tk-2.5e-3*Tk**2.)*4180.
        Tlvs(k)=Tlvl(k)+Tlsl(k)

!       tension de vapor de saturacion liquido y solido
        aux=a3+Tk*(a4+Tk*(a5+Tk*a6))
        aux=a0+Tk*(a1+Tk*(a2+Tk*aux))
        Telvs(k)=aux*100.

        aux=b3+Tk*(b4+Tk*(b5+Tk*b6))
        aux=b0+Tk*(b1+Tk*(b2+Tk*aux))
        Tesvs(k)=aux*100.        
        if (k.lt.220) then
          aux=Tlvl(220)/Rv*(1./220.-1./k)
          Telvs(k)=Telvs(220)*exp(aux)
          aux=Tlvs(220)/Rv*(1./220.-1./k)
          Tesvs(k)=Tesvs(220)*exp(aux)
        endif
!        Eautcn(k)=.3*exp(.35*Tk)*(1.+4.*exp(-.4*(k-259.)**2.))
!        Eacrcn(k)=.2*exp(.35*Tk)*(1.+4.*exp(-.4*(k-259.)**2.))
!      cambio por las expresiones de Straka(17/01/99)
        Eautcn(k)=10.**(.035*(Tk)-.7)
        Eacrcn(k)=exp(.09*Tk)


        write(3,*) k,Tvis(k),Tlvl(k),Tlsl(k),Tlvs(k),Telvs(k),Tesvs(k)
     &              ,Eautcn(k),Eacrcn(k)

  400 continue    	
      close(3)

!**   condiciones de tiempo bueno

      do 10 k=-1,nz1+2
      do 15 i=-1,nx1+2      
      do 15 j=-1,nx1+1

!     cantidades primas
	U1(i,j,k)=0.
	U2(i,j,k)=0. 
	V1(i,j,k)=0.
	V2(i,j,k)=0.     
	W1(i,j,k)=0.
	W2(i,j,k)=0.
	Pres1(i,j,k)=0.
	Pres2(i,j,k)=0.
	Titaa2(i,j,k)=0.
	Titaa1(i,j,k)=0.
	Qvap1(i,j,k)=0.
	Qvap2(i,j,k)=0.
	Qgot1(i,j,k)=0.
	Qgot2(i,j,k)=0.
	Qllu1(i,j,k)=0.
	Qllu2(i,j,k)=0.
!$$
	Qcri1(i,j,k)=0.
	Qcri2(i,j,k)=0.
	Qnie1(i,j,k)=0.
	Qnie2(i,j,k)=0.
        Qgra1(i,j,k)=0.
        Qgra2(i,j,k)=0.
	aer1(i,j,k)=0.
	aer2(i,j,k)=0.
   15 continue

!     cantidades base
	zeta=k*dx1
	Temp0(k)=TT(zeta)
	Den0(k)=Presi0(k)/Rd/Temp0(k)   
	Tita0(k)=Temp0(k)*(P00/Presi0(k))**Kapa
	Pres00(k)=Temp0(k)/Tita0(k)
	aer0(k)=10000.*exp(-zeta/2500.)
   10 continue

	Temp0(-1)=Temp0(0)
	Den0(-1)=Den0(0)
	Tita0(-1)=Tita0(0)
	Pres00(-1)=Pres00(0)
	aer0(-1)=-aer0(0)

      do 30 k=0,nz1
      do 35 i=1,nx1
      do 35 j=1,nx1

!**   perturbaciones iniciales en la temperatura y en los aerosoles

	
	zeta=k*dx1
	equis=i*dx1
	ygrie=j*dx1

	G1=exp(-((centx-equis)**2.+(centy-ygrie)**2.)*.5
     &     /radiomed**2.)

       Titaa1(i,j,k)=temper*exp(-(zeta-centz)**2./sigmat)*G1

	if (Titaa1(i,j,k).lt.1e-5) Titaa1(i,j,k)=0.

	G1=exp(-((cenaerx-equis)**2.+(cenaery-ygrie)**2.)*.5
     &     /radiomed**2.)


      aer1(i,j,k)=aerper*exp(-zeta**2./sigmaa)*G1


   35 continue


!**   vapor base

      Tem1=Temp0(k)

      if (zeta.le.500) then
	rel1=.55+.05*zeta/500.
      else if (zeta.le. 1500.) then
	  rel1=.6
      else if (zeta.le.4000) then
            rel1=.6-(zeta-1500)/2500.*.25
      else if (zeta.le.7000) then
	    rel1=.35-(zeta-4000.)/3000.*.25
      else if(zeta.gt.7000) then
	    rel1=.1-(zeta-7000)/3000.*.02
      endif


      n=int(Tem1)
      aux=Tem1-n
      elv1=Telvs(n)*(1-aux)+Telvs(n+1)*aux

       write(*,*) k,Rv,Tem1


      Qvap0(k)=rel1*elv1/Rv/Tem1      

!     recalculo de la densidad
      Den0(k)=Den0(k)+Qvap0(k)

   30 continue  

!**   Velocidad terminal para gota de lluvia, cte que depende de P
      do 450 k=1,nz1+1
        Av(2*k-1)=Av0*((P00/Presi0(k-1))**.286+    !puntos intermedios
     &            (P00/Presi0(k))**.286)/2.
        Av(2*k)=Av0*(P00/Presi0(k))**.286
  450 continue

!**   Velocidad terminal para la nieve, cte que depende de P
      do 460 k=1,nz1+1
        Vtnie(2*k-1)=Vtnie0*((P00/Presi0(k-1))**.3+  !puntos intermedio
     &            (P00/Presi0(k))**.3)/2.
        Vtnie(2*k)=Vtnie0*(P00/Presi0(k))**.3
  460 continue

!**   Velocidad terminal para el granizo, cte que depende de z
      do 470 k=0,nz1+1
        aux=2.754*rhogra**.605
        Vtgra0(2*k)=aux/Tvis(Temp0(k))**.21/Den0(k)**.395
  470 continue

      do 475 k=1,nz1+1
        Vtgra0(2*k-1)=(Vtgra0(2*k-2)+Vtgra0(2*k))/2.  ! punto intermedio
       write(*,*) 'vt',k,Vtgra0(2*k-2)*5e-3**.8,Vtgra0(2*k-1)*5e-3**.8
     &           ,Tvis(Temp0(k-1)),Den0(k-1)
  475 continue

!**************************************************************
!    Recalculo de la Presion y de Tita
   
!    Recalculo de la Presion a partir de la densidad

      call PP2(G,dx1,Den0,Presi0,P00)

      open (7,file='inic03.sa')
      do 100 k=0,nz1
	Tita0(k)=Temp0(k)*(P00/Presi0(k))**Kapa
	Pres00(k)=Temp0(k)/Tita0(k)
	cc2(k)=Cp*Rd*Tita0(k)*Pres00(k)/Cv            

	write(7,210) k,Temp0(k),Tita0(k),Presi0(k),Pres00(k),Den0(k),
     &               aer0(k),Qvap0(k),UU(k),VV(k)

  100 continue
	close(7)


      Tita0(-1)=Tita0(0)
      Pres00(-1)=Pres00(0)
      Den0(-1)=Den0(0)
      Qvap0(-1)=0

      do 300 i=1,nx1
      do 300 j=1,nx1
	Pres1(i,j,0)=Pres1(i,j,1)
	Pres1(i,j,-1)=Pres1(i,j,1)
	Pres2(i,j,0)=Pres1(i,j,1)
	Pres2(i,j,-1)=Pres1(i,j,1)
	Titaa1(i,j,0)=Titaa1(i,j,1)
	Titaa1(i,j,-1)=Titaa1(i,j,1)
	Qvap1(i,j,0)=Qvap1(i,j,1)
	Qvap1(i,j,-1)=Qvap1(i,j,1)
  300 continue

!     calculo del Qvaprel
      Qvaptot=0.
      do 350 k=1,nz1
        Qvaptot=Qvaptot+Qvap0(k)
  350 continue
      do 355 k=1,nz1
        Qvaprel(k)=Qvap0(k)/Qvaptot
  355 continue

!     calculo del aerrel
      aertot=0.
      do 360 k=1,nz1
        aertot=aertot+aer0(k)
  360 continue
      do 365 k=1,nz1
        aerrel(k)=aer0(k)/aertot
  365 continue

  200 format(I3,4E11.3)   
  210 format(I3,9E12.4)
 
      return
      
      end       

!*********************************************************

      real function TT (zeta)
	real a,zeta,xx
	a= 298.15
	if (zeta.le.2000) then
	  TT=a-9.e-3*zeta
	elseif (zeta.le.5500) then
	  xx=zeta-2000.
	  TT=a-18.-xx*(9.e-3-2e-3*xx/3500./2.)
        elseif (zeta.le.9000) then
	  xx=zeta-5500.
	  TT=a-46.-7e-3*xx
        elseif (zeta.le.11000) then
          xx=zeta-9000
          TT=a-70.5-7e-3*xx+1.75e-6*xx**2.
        elseif (zeta.le.12000) then
          TT=a-77.5
        else
          xx=zeta-12000
          TT=a-77.5+50.*(xx/9000.)**2.
        endif

      end
!*****************************************************
      subroutine PP(G,Rd,dx,nz1,Pres,Pres0)
	integer k,nz1,nx4
	parameter (nx4=500)
	real TT
	real Pres(-3:nz1+3)
	real integ(-2:nx4+2)
	real Pres0
	real G,Rd,dx,dx4
	real zetaa,zetam,zetad
        dx4=dx/4.
	
	integ(0)=0         
      do 10 k=1,nx4
	  zetaa=(2*k-2)*dx4
	  zetam=(2*k-1)*dx4
	  zetad=(2*k)*dx4
	  ya=1/TT(zetaa)
	  ym=1/TT(zetam)
	  yd=1/TT(zetad)
	  integ(k)=integ(k-1)+ya+4*ym+yd
   10 continue      

      do 20 k=1,nz1+2
	 Pres(k)=Pres0*exp(-G/Rd*(integ(2*k)*dx4/3))
   20 continue
      Pres(0)=Pres0
       Pres(-1)=Pres0

       return        
       end

!*****************************************************
      subroutine PP2(G,dx,Den0,Pres00,Pres0)

	include 'dimen.i'

	integer k
	real Pres00(-3:nz1+3)
	real Den0(-3:nz1+3)
	real Den00(-3:3*nz1+3)
	real integ(-3:3*nz1+3)
	real Pres0
	real G,dx

      do 5 k=0,nz1-1
	  Den00(2*k)=Den0(k)
	  Den00(2*k+1)=(Den0(k)+Den0(k+1))/2.
    5 continue
	  Den00(2*nz1)=Den0(nz1)
	  Den00(2*nz1+1)=2.*Den0(nz1)-Den00(2*nz1-1)

	integ(0)=0         
      do 10 k=1,nz1
	  ya=Den00(2*k-1)
	  ym=Den00(2*k)
	  yd=Den00(2*k+1)
	  integ(k)=integ(k-1)+ya+4*ym+yd
   10 continue      
      do 20 k=1,nz1
	Pres00(k)=Pres0-G*integ(k)*dx/6.

   20 continue
      Pres00(0)=Pres0
      Pres00(-1)=Pres0

      return
      end
