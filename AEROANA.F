*     Este programa lee la salida del la nube en 3-D y realiza calculos
*     para indicar la eliminacion y liberacion de aerosoles por las
*     particulas

      program aeroana
      implicit none

      include 'dimen.i'
      include 'const.i'
      include 'perdim.i'
      include 'permic.i'
      include 'estbas.i'
      include 'cant01.i'
      include 'advecs.i'
      include 'lmngot.i'
      include 'lmnllu.i'
      include 'lmncri.i'
      include 'lmnnie.i'
      include 'lmngra.i'
      include 'aeroana.i'


      integer ii,ymax
      integer pasoaux(25),lvar,nvar
      parameter (nvar=13)
      character*100 paso
      character*2 nombre
      character*25 arch
      character*2 var(13)
      integer*4 variable(nx1,nz1,nvar)


***   datos generales
      nombre='21'
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

**    lectura de los datos generales y de base

      open(15,file='aerores')
      open(13,file='inis.da')
      open(unit=14,file='velos.da',status='unknown',form='unformatted')
      rewind 14
      open(unit=3,file='varconz.da',status='unknown',form='unformatted')
      rewind 3

	  read(13,*) Den0,Temp0,Tita0,Pres00,Qvap0
     &                  ,cc2,aer0,UU,VV

*$$
	  read(3)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,
     &             Qvaprel,aerrel,Eautcn,Eacrcn

	close(3)
	close(13)

***   definicion de calores y presiones de vapor a 0 K

      Lsl00=Lsl0*4180.

**********************************************************

      do 10 ii=1,25
        if(pasoaux(ii).eq.1) then

**        lectura de la nube

          arch='d:\nube'//nombre//'\nube'//nombre
     &          //paso(ii*2-1:ii*2)//'.sal'

          write(*,*) arch

          open(unit=6,file=arch,status='unknown',form='unformatted')
          read(6) U2,V2,W2,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1
     &          ,Qgra1,aer1
          close(6)


**    inicializado para cada paso

	s=0
	Qvapneg=0.
	lvapneg=0
	aerneg=0.
	laerneg=0
	llluneg=0
*$$
        lcrineg=0.
        lnieneg=0.
        lgraneg=0.
	ener=0.
	ener1=0.        
	ener2=0.     
	ener3=0.     
	ener4=0.     
	ener5=0.     
	qv=0.
	qg=0.
	daitot=0.     
	lgot(1)=nx1
	lgot(2)=0
	mgot(1)=nx1
	mgot(2)=0
	ngot(1)=nz1
	ngot(2)=0
	lllu(1)=nx1
	lllu(2)=0
	mllu(1)=nx1
	mllu(2)=0
	nllu(1)=nz1
	nllu(2)=0
*$$
        lcri(1)=nx1
        lcri(2)=0
        mcri(1)=nx1
        mcri(2)=0
        ncri(1)=nz1
        ncri(2)=0
        lnie(1)=nx1
        lnie(2)=0
        mnie(1)=nx1
        mnie(2)=0
        nnie(1)=nz1
        nnie(2)=0
        lgra(1)=nx1
        lgra(2)=0
        mgra(1)=nx1
        mgra(2)=0
        ngra(1)=nz1
        ngra(2)=0
	vapt1=0.
	vapt2=0.
	vapt3=0.
	vapt4=0.
	gott1=0.
	gott2=0.
	gott3=0.
	gott4=0.
	aert1=0.
	aert2=0.
	aert3=0.
	aert4=0.
	totnuc=0.
	totmic=0.

	do 15 i=0,nx1+1
	do 15 j=0,nx1+1
	  advvap1(i,j)=W2(i,j,1)*(Qvap1(i,j,1)+Qvap1(i,j,0))/4.
	  advgot1(i,j)=0.
	  advllu1(i,j)=W2(i,j,1)*Qllu1(i,j,1)

*     modifique i,j,k por i,j,1 (11/01/99)

	  if (W2(i,j,1).gt.0) advllu1(i,j)=0.
	  advaer1(i,j)=W2(i,j,1)*(aer1(i,j,1)+aer1(i,j,0))/4.
	  if(W2(i,j,1).lt.0) advaer1(i,j)=advaer1(i,j)*1.5
*$$
	  advcri1(i,j)=0.
	  advnie1(i,j)=W2(i,j,1)*Qnie1(i,j,1)
	  if (W2(i,j,1).gt.0) advnie1(i,j)=0.
          advgra1(i,j)=W2(i,j,1)*Qgra1(i,j,1)
          if (W2(i,j,1).gt.0) advgra1(i,j)=0.
   15   continue

*********************************************************************
**    calculo de la dinamica y de la termodinamica

	do 20 k=1,nz1-1
	  n=k
	  dden0z=(Den0(k+1)-Den0(k-1))/Den0(k)    
	  call turbu1(n)
	do 20 i=1,nx1
	  l=i
	do 20 j=1,nx1
	  m=j

*    calculo del coeficiente de turbulencia y derivadas
	  call turbu2(l,m,n)

*    calculo de las inhomogeneidades para las velocidades
*$$
	  call inomo(l,m,n,dden0z)    

***     calculo de la energia cinetica
	ener1=.5*Den0(k)*(U2(i,j,k)**2.+V2(i,j,k)**2.+W2(i,j,k)**2.)
     &       +ener1

*     calculo de la temperatura potencial
	  call tempot(l,m,n,dden0z,Fcalo(i,j,k))
	  Fcalo(i,j,k)=0.
*     dinamica del vapor y de las gotitas

	  call dvapor(l,m,n)
	    advvap1(i,j)=advvap2(i,j)

	  call dgotit(l,m,n)
	    advgot1(i,j)=advgot2(i,j)

	  call dlluvi(l,m,n)

	    advllu1(i,j)=advllu2(i,j)
	  

*&&

	  call dcrist(l,m,n)
	    advcri1(i,j)=advcri2(i,j)
          call dnieve(l,m,n)
	    advnie1(i,j)=advnie2(i,j)
          call dgrani(l,m,n)
            advgra1(i,j)=advgra2(i,j)
	  
	  call daeros(l,m,n)
	    advaer1(i,j)=advaer2(i,j)

*     limites de la nube
	  if(Qgot2(i,j,k).ne.0) then
	    if (i.lt.lgot(1)) lgot(1)=i
	    if (i.gt.lgot(2)) lgot(2)=i
	    if (j.lt.mgot(1)) mgot(1)=j
	    if (j.gt.mgot(2)) mgot(2)=j
	    if (k.lt.ngot(1)) ngot(1)=k
	    if (k.gt.ngot(2)) ngot(2)=k
	    s=s+1
	  endif

*     limites de la lluvia
	  if(Qllu2(i,j,k).ne.0) then
	    if (i.lt.lllu(1)) lllu(1)=i
	    if (i.gt.lllu(2)) lllu(2)=i
	    if (j.lt.mllu(1)) mllu(1)=j
	    if (j.gt.mllu(2)) mllu(2)=j
	    if (k.lt.nllu(1)) nllu(1)=k
	    if (k.gt.nllu(2)) nllu(2)=k
	    llluneg=1
	  endif

*     limites de los cristales
*$$
          if(Qcri2(i,j,k).ne.0) then
            if (i.lt.lcri(1)) lcri(1)=i
            if (i.gt.lcri(2)) lcri(2)=i
            if (j.lt.mcri(1)) mcri(1)=j
            if (j.gt.mcri(2)) mcri(2)=j
            if (k.lt.ncri(1)) ncri(1)=k
            if (k.gt.ncri(2)) ncri(2)=k
            lcrineg=1
	  endif

*     limites de la nieve
          if(Qnie2(i,j,k).ne.0) then
            if (i.lt.lnie(1)) lnie(1)=i
            if (i.gt.lnie(2)) lnie(2)=i
            if (j.lt.mnie(1)) mnie(1)=j
            if (j.gt.mnie(2)) mnie(2)=j
            if (k.lt.nnie(1)) nnie(1)=k
            if (k.gt.nnie(2)) nnie(2)=k
            lnieneg=1
	  endif

*     limites del granizo
          if(Qgra2(i,j,k).ne.0) then
            if (i.lt.lgra(1)) lgra(1)=i
            if (i.gt.lgra(2)) lgra(2)=i
            if (j.lt.mgra(1)) mgra(1)=j
            if (j.gt.mgra(2)) mgra(2)=j
            if (k.lt.ngra(1)) ngra(1)=k
            if (k.gt.ngra(2)) ngra(2)=k
            lgraneg=1
	  endif

	  if(Qvap0(k)+Qvap2(i,j,k).lt.0) then
	    write(*,*) 'moco con el vapor',l,m,n,Qvap0(k)+Qvap2(i,j,k)
	    write(*,*) Qvap2(i,j,k),Qvap0(k),P,T,W2(i,j,k)
	    Qvapneg=Qvapneg+Qvap0(k)+Qvap2(i,j,k)
	    lvapneg=1
*            pause
	  endif

	  if(aer0(k)+aer2(i,j,k).lt.0) then
	    write(*,*) 'moco con aerosoles',l,m,n,aer0(k)+aer2(i,j,k)
	    write(*,*) aer2(i,j,k),aer0(k),P,T,W2(i,j,k)
	    aerneg=aerneg+aer0(k)+aer2(i,j,k)
	    laerneg=1
*            pause
	  endif

   20   continue

**    correccion de negativos
	if(s.ge.1) call corgot
	if (llluneg.eq.1) call corllu
*$$
        if (lcrineg.eq.1) call corcri
        if (lnieneg.eq.1) call cornie
        if (lgraneg.eq.1) call corgra
        if (lvapneg.eq.1) call corvap(Qvapneg)
	if (laerneg.eq.1) call coraer(aerneg)

*#### primer calculo de agua (sin laterales)
      do 23 i=1,nx1
      do 23 j=1,nx1
      do 23 k=1,nz1-1
	vapt1=vapt1+Qvap2(i,j,k)
	gott1=gott1+Qgot2(i,j,k)
	aert1=aert1+aer2(i,j,k)

      if(Qvap2(i,j,k)+Qvap0(k).lt.0) then
        write(*,*) 'vapomoco',i,j,k,Qvap2(i,j,k),Qvap0(k),Qvap1(i,j,k)
        stop
      endif

   23 continue
*###

      write(*,*) 'antes de microfis'

*********************************************************************

*********************************************************************
**    sublazo microfisico       

	do 30 k=1,nz1-1
	  n=k
	do 30 i=1,nx1
	  l=i
	do 30 j=1,nx1
	  m=j

*     calculo de T,P,Densi,Dv,Vis

	  aux=Pres00(k)+Pres2(i,j,k)
	  P=aux**ikapa*P00
	  T=(Tita0(k)+Titaa2(i,j,k))*aux
	  Tempa1(i,j,k)=T-Temp0(k)
	  Qvap=Qvap0(k)+Qvap2(i,j,k)
          Naer=aer0(k)+aer2(i,j,k)
	  densi=P/T/Rd-AA*Qvap
	  Dv=Dv0*(T/273.15)**1.94*(P00/P)

*     calculo de Vis, Lvl, Lsl, Lvs, elvs y  esvs

	  iT=int(T)
	  aux2=T-iT
	  Vis=Tvis(iT)*(1-aux2)+Tvis(iT+1)*aux2
	  elvs=Telvs(iT)*(1-aux2)+Telvs(iT+1)*aux2
	  esvs=Tesvs(iT)*(1-aux2)+Tesvs(iT+1)*aux2
	  Lvl=Tlvl(iT)*(1-aux2)+Tlvl(iT+1)*aux2
	  Lsl=Tlsl(iT)*(1-aux2)+Tlsl(iT+1)*aux2
	  Lvs=Tlvs(iT)*(1-aux2)+Tlvs(iT+1)*aux2
	  Eaccn=Eacrcn(iT)*(1-aux2)+Eacrcn(iT+1)*aux2
	  Eaucn=Eautcn(iT)*(1-aux2)+Eautcn(iT+1)*aux2
          if (T.ge.T0) then
            Eacng=1.
          else
            Eacng=exp(.08*(T-T0)) !mod 6/2/2000
          endif

	  nu=Vis/densi

**    nucleacion (de ser necesario tiene otro paso de tiempo)
      lll=tt

	  Qliq=Qgot2(i,j,k)
	  e1=Qvap*Rv*T
          rl=(e1-elvs)/elvs
          rs=(e1-esvs)/esvs
          yy=0

        if ((rl.gt.1e-3 .or. rs.gt.1e-3).and. Naer.gt.0) then
       
          call nuclea(Qvap,Qliq,Naer,T,densi,e1,elvs,esvs,
     &                rl,rs,Lvl,Lvs,l,m,n,daer,dqgot,dqcri)


	  Taux=T-Temp0(k)-Tempa1(i,j,k)
	  Titaa2(i,j,k)=T/aux-Tita0(k)
          if (dqgot.gt.0) yy=1
	else    
          Taux=0.
	  dqgot=0.
          dqcri=0.
	  daer=0.
        endif   

	totnuc=totnuc+daer

*#### segundo calculo de agua (sin laterales)
	vapt2=vapt2+Qvap2(i,j,k)
	gott2=gott2+Qgot2(i,j,k)
	aert2=aert2+aer2(i,j,k)

*$$
	if (Qgot2(i,j,k).gt.0 .or. dqgot.gt.0 .or. 
     &      Qllu2(i,j,k).gt.0 .or. Qcri2(i,j,k).gt.0 .or.
     &      Qnie2(i,j,k).gt.0) then

	  qgotaux=Qgot2(i,j,k)
          if (Qgot2(i,j,k).eq.0) qgotaux=0d0
	  qvapaux=Qvap2(i,j,k)+Qvap0(k)
	  qlluaux=Qllu2(i,j,k)
          if (Qllu2(i,j,k).eq.0) qlluaux=0d0
	  qcriaux=Qcri2(i,j,k)
          if (Qcri2(i,j,k).eq.0) then
            qcriaux=0d0
          endif

          qnieaux=Qnie2(i,j,k)
          if (Qnie2(i,j,k).eq.0) qnieaux=0d0
          qgraaux=Qgra2(i,j,k)
          if (Qgra2(i,j,k).eq.0) qgraaux=0d0
          Naer=aer2(i,j,k)+aer0(k)
	  T=Tempa1(i,j,k)+Temp0(k)

	  do 35 t2=1,lt2              
	    qgotaux=qgotaux+dqgot/float(lt2)
            qcriaux=qcriaux+dqcri/float(lt2)
            qvapaux=qvapaux-(dqgot+dqcri)/float(lt2)
	    Naer=Naer+daer/float(lt2)
	    T=T+Taux/float(lt2)

*     calculo de elvs y esvs

	    iT=int(T)
	    aux2=T-iT
	    elvs=Telvs(iT)*(1-aux2)+Telvs(iT+1)*aux2
	    esvs=Tesvs(iT)*(1-aux2)+Tesvs(iT+1)*aux2

	    call microfis(elvs,esvs,Lvl,Lvs,Lsl,T,Dv,Eaccn,Eaucn,
     &                    Eacng,Lsl00,Fcal,l,m,n,
     &                    qvapaux,qgotaux,qlluaux,qcriaux,qnieaux,
     &                    qgraaux,Naer,daer2,nu,yy)

	    Fcalo(l,m,n)=Fcalo(l,m,n)+Fcal/dt1/densi
	    Naer=Naer+daer2

	    totmic=totmic+daer2

 35       continue

*$$
	  Qgot2(i,j,k)=qgotaux
	  Qllu2(i,j,k)=qlluaux
	  Qcri2(i,j,k)=qcriaux
	  Qnie2(i,j,k)=qnieaux
          Qgra2(i,j,k)=qgraaux
	  Qvap2(i,j,k)=qvapaux-Qvap0(k)
	  aer2(i,j,k)=Naer-aer0(k) !mod 14/9/99
	  Tempa1(i,j,k)=T-Temp0(k)

	endif


       if (Tita0(k).lt.abs(Titaa2(i,j,k))+200 .or. Temp0(k).lt.
     &   abs(Tempa1(i,j,k))+200) then
       write(*,*) 'problemas con la temperatura 2'
       write(*,*) i,j,k,Titaa2(i,j,k),Tita0(k),Tempa1(i,j,k),Temp0(k)
       write(*,*) T,Titaa2(i,j,k),aux
       stop
       endif

       if(aer2(i,j,k)+aer0(k).le.0) then

         if (aer2(i,j,k)+aer0(k).lt.-aer0(k)*.05) then

	 write(*,*) 'problemas con los aerosoles'
	 write(*,*) i,j,k,aer1(i,j,k),aer2(i,j,k),aer0(k),daer

           write(*,*) 'sin aerosoles',W2(i,j,k)
           stop
         endif

	 aer2(i,j,k)=-aer0(k)
       endif

*#### tercer calculo de agua (sin laterales)
	vapt3=vapt3+Qvap2(i,j,k)
	gott3=gott3+Qgot2(i,j,k)
	aert3=aert3+aer2(i,j,k)
*###


***   calculo de la energia
      ener2=densi*G*k*dx1+ener2
      ener3=densi*(Cp-Rd)*T+ener3
      ener4=P+ener4
      ener5=(Qvap2(i,j,k)+Qgot2(i,j,k))*G*k*dx1+ener5
      qv=Qvap2(i,j,k)+qv
      qg=Qgot2(i,j,k)+qg
      daitot=densi+daitot
***
         aerdif(i,j,k)=aer2(i,j,k)-aer1(i,j,k)

   30   continue


******************************************************
*                 grabacion

        open(unit=10,file='aerdif'//paso(ii*2-1:ii*2),
     &       status='unknown',form='unformatted')
        write(10) aerdif
        close(10)

        write(15,*) ii,totnuc,totmic

      endif
  10  continue

      end
