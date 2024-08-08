!23456789*123456789*123456789*123456789*123456789*123456789*123456789*
!     Revision 8/04/00
!     Incluye microfisica con vapor, gotitas, lluvia, cristales, nieve
!      y granizos (por unidad de volumen)
!     Con viento de corte.
!     Este modelo simula una nube tridimensional, con diferencias
!      finitas adelantadas en el tiempo y centradas en el espacio
!     Este modelo (asi como sus variantes) sirve de test
!     Las grillas son las mismas para las cantidades dinamicas y
!      microfisicas, como asi tambien los intervalos de tiempo.
!     Todas las variables son reales*4
!     Condiciones de contorno homogeneas para las variables microfisicas
!     Graba el valor de todas las variables cada ltb segundos.
!     Graba el valor de las variables para analisis cada ltg segundos.
!     Condicion de contorno nula para el vapor
!     Contempla el desplazamiento de la nube.
!     Mejora la condicion en el piso para los aerosoles cuando hay agua (815)


program modelo

   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE lmngot
   USE lmnllu
   USE lmncri
   USE lmnnie
   USE lmngra
   USE turbvar
   USE mein
   USE mode20
   USE estad03
   USE posnub02
   USE corrinu2
   implicit none

   call mode20_init()

!########### Condiciones iniciales ############
   if (ini.eq.0) then
      !Si ini=0 el calculo empieza por primera###
      call condi
   else
      !si ini=1 el calculo recomienza desde algun paso
      open(unit=13,file='outputdata/inis.da',status='unknown',form='unformatted')
      open(unit=14,file='outputdata/velos.da',status='unknown',form='unformatted')
      rewind 14
      open(unit=30,file='outputdata/varconz.da',status='unknown',form='unformatted')
      rewind 30

      read(13,*) Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV
      read(14) U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,&
         Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
         Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,&
         aer1,aer2,Fcalo
      read(30)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn

      close(30)
      close(13)
      close(14)

      do 5001 i=1,nx1
         do 5001 j=1,nx1
            do 5001 k=1,nz1
               aux=.8*exp(-((i-35)**2.+(j-25.5)**2.+(k-8)**2.)/50.)
               Titaa1(i,j,k)=Titaa1(i,j,k)+aux
               Titaa2(i,j,k)=Titaa2(i,j,k)+aux
5001  continue
   endif

   !definicion de calores y presiones de vapor a 0 K
   Lsl00=Lsl0*4180.
   write(*,*) 'Qgot',Qgot1(32,34,0),Qgot1(32,34,1),Qgot1(32,34,2)
   write(*,*) 'Qcri',Qcri1(32,34,0),Qcri1(32,34,1),Qcri1(32,34,2)
   write(*,*) 'Qnie',Qnie1(32,34,0),Qnie1(32,34,1),Qnie1(32,34,2)
   write(*,*) 'Qllu',Qllu1(32,34,0),Qllu1(32,34,1),Qllu1(32,34,2)
   write(*,*) Qllu1(16,12,27),Qllu1(12,17,27)
   write(*,*) 'des cond',Qvap1(20,20,8),Qvap0(8)
   write(*,*) Tita0(2),Tita0(3),Tita0(4)
   write(*,*) Tita0(9),Tita0(10),Tita0(11)
   write(*,*) Tita0(20),Tita0(17),Tita0(18)

   open(unit=13,file='outputdata/u2'//tie//bre(2*t1+1:2*t1+2))
   open(unit=14,file='outputdata/v2'//tie//bre(2*t1+1:2*t1+2))
   open(unit=15,file='outputdata/w2'//tie//bre(2*t1+1:2*t1+2))
   open(unit=16,file='outputdata/va'//tie//bre(2*t1+1:2*t1+2))
   open(unit=17,file='outputdata/go'//tie//bre(2*t1+1:2*t1+2))
   open(unit=18,file='outputdata/ae'//tie//bre(2*t1+1:2*t1+2))
   open(unit=19,file='outputdata/pr'//tie//bre(2*t1+1:2*t1+2))
   open(unit=20,file='outputdata/ti'//tie//bre(2*t1+1:2*t1+2))

   open(unit=30,file='outputdata/esta'//tie//bre(2*t1+1:2*t1+2))
   open(unit=31,file='outputdata/vara'//tie//bre(2*t1+1:2*t1+2))
   open(unit=32,file='outputdata/posnub'//tie//'.sa')
   open(unit=33,file='outputdata/est'//tie//bre(2*t1+1:2*t1+2))

   write(*,*) 'inis',Titaa1(16,16,4),Qvap1(16,16,4),Tita0(4),Qvap0(4)
   write(*,*) 'tem',Fcalo(21,16,5)

!#####################################################################
!############### comienzo de la evolucion temporal ###################
   !lazo temporal principal
   do 1 tt=1,lt1
      !inicializado para cada paso
      s=0
      Qvapneg=0.
      lvapneg=0
      aerneg=0.
      laerneg=0
      llluneg=0
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

!######################## Advección de vapores #######################
      
      do 15 i=0,nx1+1
         do 15 j=0,nx1+1
            advvap1(i,j)=W2(i,j,1)*(Qvap1(i,j,1)+Qvap1(i,j,0))/4.
            advgot1(i,j)=0.
            advllu1(i,j)=W2(i,j,1)*Qllu1(i,j,1)

            if (W2(i,j,1).gt.0) advllu1(i,j)=0.
            advaer1(i,j)=W2(i,j,1)*(aer1(i,j,1)+aer1(i,j,0))/4.
            if(W2(i,j,1).lt.0) advaer1(i,j)=advaer1(i,j)*1.5
            advcri1(i,j)=0.
            advnie1(i,j)=W2(i,j,1)*Qnie1(i,j,1)
            if (W2(i,j,1).gt.0) advnie1(i,j)=0.
            advgra1(i,j)=W2(i,j,1)*Qgra1(i,j,1)
            if (W2(i,j,1).gt.0) advgra1(i,j)=0.
15    continue

!########### calculo de la dinamica y de la termodinamica ############

      do 20 k=1,nz1-1
         n=k
         dden0z=(Den0(k+1)-Den0(k-1))/Den0(k)
         call turbu1(n)
         do 20 i=1,nx1
            l=i
            do 20 j=1,nx1
               m=j
               !calculo del coeficiente de turbulencia y derivadas
               call turbu2(l,m,n)

               !calculo de las inhomogeneidades para las velocidades
               call inomo(l,m,n,dden0z)

               !calculo de la energia cinetica
               ener1=.5*Den0(k)*(U2(i,j,k)**2.+V2(i,j,k)**2.+W2(i,j,k)**2.)+ener1

               !calculo de la temperatura potencial
               call tempot(l,m,n,dden0z,Fcalo(i,j,k))
               Fcalo(i,j,k)=0.

               !dinamica del vapor y de las gotitas
               call dvapor(l,m,n)
               advvap1(i,j)=advvap2(i,j)
               call dgotit(l,m,n)
               advgot1(i,j)=advgot2(i,j)
               call dlluvi(l,m,n)
               advllu1(i,j)=advllu2(i,j)
               call dcrist(l,m,n)
               advcri1(i,j)=advcri2(i,j)
               call dnieve(l,m,n)
               advnie1(i,j)=advnie2(i,j)
               call dgrani(l,m,n)
               advgra1(i,j)=advgra2(i,j)

               if  (i.eq.20.and.m.eq.22.and.(k.ge.27.and. k.le.28)) then
                  write(*,*) 'cri0',k,Qcri1(i,j,k),Qcri2(i,j,k)
               endif

               call daeros(l,m,n)
               advaer1(i,j)=advaer2(i,j)

               !limites de la nube
               if(Qgot2(i,j,k).ne.0) then
                  if (i.lt.lgot(1)) lgot(1)=i
                  if (i.gt.lgot(2)) lgot(2)=i
                  if (j.lt.mgot(1)) mgot(1)=j
                  if (j.gt.mgot(2)) mgot(2)=j
                  if (k.lt.ngot(1)) ngot(1)=k
                  if (k.gt.ngot(2)) ngot(2)=k
                  s=s+1
               endif

               !limites de la lluvia
               if(Qllu2(i,j,k).ne.0) then
                  if (i.lt.lllu(1)) lllu(1)=i
                  if (i.gt.lllu(2)) lllu(2)=i
                  if (j.lt.mllu(1)) mllu(1)=j
                  if (j.gt.mllu(2)) mllu(2)=j
                  if (k.lt.nllu(1)) nllu(1)=k
                  if (k.gt.nllu(2)) nllu(2)=k
                  llluneg=1
               endif

               !limites de los cristales
               if(Qcri2(i,j,k).ne.0) then
                  if (i.lt.lcri(1)) lcri(1)=i
                  if (i.gt.lcri(2)) lcri(2)=i
                  if (j.lt.mcri(1)) mcri(1)=j
                  if (j.gt.mcri(2)) mcri(2)=j
                  if (k.lt.ncri(1)) ncri(1)=k
                  if (k.gt.ncri(2)) ncri(2)=k
                  lcrineg=1
               endif

               !limites de la nieve
               if(Qnie2(i,j,k).ne.0) then
                  if (i.lt.lnie(1)) lnie(1)=i
                  if (i.gt.lnie(2)) lnie(2)=i
                  if (j.lt.mnie(1)) mnie(1)=j
                  if (j.gt.mnie(2)) mnie(2)=j
                  if (k.lt.nnie(1)) nnie(1)=k
                  if (k.gt.nnie(2)) nnie(2)=k
                  lnieneg=1
               endif

               !limites del granizo
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
               endif

               if(aer0(k)+aer2(i,j,k).lt.0) then
                  write(*,*) 'moco con aerosoles',l,m,n,aer0(k)+aer2(i,j,k)
                  write(*,*) aer2(i,j,k),aer0(k),P,T,W2(i,j,k)
                  aerneg=aerneg+aer0(k)+aer2(i,j,k)
                  laerneg=1
               endif
20    continue

      write(*,3000) 'limites',lnie(1),lnie(2),mnie(1),mnie(2),nnie(1),nnie(2)

      !correccion de negativos
      if(s.ge.1) call corgot
      if (llluneg.eq.1) call corllu
      if (lcrineg.eq.1) call corcri
      if (lnieneg.eq.1) call cornie
      if (lgraneg.eq.1) call corgra
      if (lvapneg.eq.1) call corvap(Qvapneg)
      if (laerneg.eq.1) call coraer(aerneg)
      write(*,*) 'vap0',Qvap2(16,16,25),Qvap0(25)
      write(*,*) 'vap0',Qvap2(16,16,26),Qvap0(26)
      write(*,*) 'vap0',Qvap2(16,16,27),Qvap0(27)

      !primer calculo de agua (sin laterales)
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
23    continue

!####################### Sublazo Microfísico #########################
      write(*,*) 'antes de microfis'
      do 30 k=1,nz1-1
         n=k
         do 30 i=1,nx1
            l=i
            do 30 j=1,nx1
               m=j

               !calculo de T,P,Densi,Dv,Vis
               aux=Pres00(k)+Pres2(i,j,k)
               P=aux**ikapa*P00
               T=(Tita0(k)+Titaa2(i,j,k))*aux
               Tempa1(i,j,k)=T-Temp0(k)
               Qvap=Qvap0(k)+Qvap2(i,j,k)
               Naer=aer0(k)+aer2(i,j,k)
               densi=P/T/Rd-AA*Qvap
               Dv=Dv0*(T/273.15)**1.94*(P00/P)

               !calculo de Vis, Lvl, Lsl, Lvs, elvs y  esvs
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
                  Eacng=exp(.08*(T-T0))
               endif

               nu=Vis/densi

               !nucleacion (de ser necesario tiene otro paso de tiempo)
               lll=tt

               Qliq=Qgot2(i,j,k)
               e1=Qvap*Rv*T
               rl=(e1-elvs)/elvs
               rs=(e1-esvs)/esvs
               yy=0

               if ((rl.gt.1e-3 .or. rs.gt.1e-3).and.Naer.gt.0) then
                  call nuclea(Qvap,Qliq,Naer,T,densi,e1,elvs,esvs,rl,rs,Lvl,Lvs,l,m,n,daer,dqgot,dqcri)
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

               !segundo calculo de agua (sin laterales)
               vapt2=vapt2+Qvap2(i,j,k)
               gott2=gott2+Qgot2(i,j,k)
               aert2=aert2+aer2(i,j,k)
               if (Qgot2(i,j,k).gt.0 .or. dqgot.gt.0 .or.Qllu2(i,j,k).gt.0 .or. Qcri2(i,j,k).gt.0 .or.Qnie2(i,j,k).gt.0) then
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

                     !calculo de elvs y esvs
                     iT=int(T)
                     aux2=T-iT
                     elvs=Telvs(iT)*(1-aux2)+Telvs(iT+1)*aux2
                     esvs=Tesvs(iT)*(1-aux2)+Tesvs(iT+1)*aux2
                     call microfis(elvs,esvs,Lvl,Lvs,Lsl,T,Dv,Eaccn,Eaucn,Eacng,Lsl00,&
                        Fcal,l,m,n,qvapaux,qgotaux,qlluaux,qcriaux,qnieaux,&
                        qgraaux,Naer,daer2,nu,yy)
                     Fcalo(l,m,n)=Fcalo(l,m,n)+Fcal/dt1/densi
                     Naer=Naer+daer2
                     totmic=totmic+daer2
                     if (i.eq.20.and.(j.eq.22.or.j.eq.22).and.k.eq.28) then
                        write(*,*) qgotaux,qlluaux,qcriaux
                     endif
35                continue

                  Qgot2(i,j,k)=qgotaux
                  Qllu2(i,j,k)=qlluaux
                  Qcri2(i,j,k)=qcriaux
                  Qnie2(i,j,k)=qnieaux
                  Qgra2(i,j,k)=qgraaux
                  Qvap2(i,j,k)=qvapaux-Qvap0(k)
                  aer2(i,j,k)=Naer-aer0(k)
                  Tempa1(i,j,k)=T-Temp0(k)

               endif

               if (Tita0(k).lt.abs(Titaa2(i,j,k))+200.or.Temp0(k).lt.abs(Tempa1(i,j,k))+200) then
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

               !tercer calculo de agua (sin laterales)
               vapt3=vapt3+Qvap2(i,j,k)
               gott3=gott3+Qgot2(i,j,k)
               aert3=aert3+aer2(i,j,k)

               !calculo de la energia
               ener2=densi*G*k*dx1+ener2
               ener3=densi*(Cp-Rd)*T+ener3
               ener4=P+ener4
               ener5=(Qvap2(i,j,k)+Qgot2(i,j,k))*G*k*dx1+ener5
               qv=Qvap2(i,j,k)+qv
               qg=Qgot2(i,j,k)+qg
               daitot=densi+daitot
30    continue

!####################### Contornos ##############################
      Qvap=(qv+qg)/nx1**2.*nz1*.1
      write(*,*) 'antes de las redefiniciones'

      !contornos en el piso y en el techo
      do 400 i=1,nx1
         do 400 j=1,nx1
            Titaa2(i,j,0)=-W2(i,j,1)*(Tita0(0)+Tita0(1))*dt1/dx2+Titaa1(i,j,0)
            Titaa2(i,j,nz1)=Titaa2(i,j,nz1-1)

            !suponemos que las velocidades horizontales a nivel de piso son
            !iguales a 1/4 de la correspondiente en el nivel 1
            auxx=((U2(i+1,j,1)+U2(i,j,1))*(Qvap1(i+1,j,0)+Qvap1(i,j,0))&
               -(U2(i-1,j,1)+U2(i,j,1))*(Qvap1(i-1,j,0)+Qvap1(i,j,0)))&
               /4.*.25
            auxy=((V2(i,j+1,1)+V2(i,j,1))*(Qvap1(i,j+1,0)+Qvap1(i,j,0))&
               -(V2(i,j-1,1)+V2(i,j,1))*(Qvap1(i,j-1,0)+Qvap1(i,j,0)))&
               /4.*.25
            auxz=W2(i,j,1)*((Qvap1(i,j,1)+Qvap1(i,j,0))+Qvap0(1)+Qvap0(0))/2.*.5

            Qvap2(i,j,nz1)=Qvap2(i,j,nz1-1)

            Qgot2(i,j,0)=Qgot2(i,j,1)
            Qgot2(i,j,nz1)=Qgot2(i,j,nz1-1)

            Qcri2(i,j,0)=Qcri2(i,j,1)
            Qcri2(i,j,nz1)=Qcri2(i,j,nz1-1)

            !Para que no se acumulen el piso
            Qllu2(i,j,0)=Qllu2(i,j,1)/2.
            Qllu2(i,j,nz1)=Qllu2(i,j,nz1-1)

            Qnie2(i,j,0)=Qnie2(i,j,1)
            Qnie2(i,j,nz1)=Qnie2(i,j,nz1-1)

            !Para que no se acumulen el piso
            Qgra2(i,j,0)=Qgra2(i,j,1)/2.
            Qgra2(i,j,nz1)=Qgra2(i,j,nz1-1)

            !suponemos que las velocidades horizontales a nivel de piso son
            !iguales a 1/4 de la correspondiente en el nivel 1

            auxx=((U2(i+1,j,1)+U2(i,j,1))*(aer1(i+1,j,0)+aer1(i,j,0))&
               -(U2(i-1,j,1)+U2(i,j,1))*(aer1(i-1,j,0)+aer1(i,j,0)))&
               /4.*.25
            auxy=((V2(i,j+1,1)+V2(i,j,1))*(aer1(i,j+1,0)+aer1(i,j,0))&
               -(V2(i,j-1,1)+V2(i,j,1))*(aer1(i,j-1,0)+aer1(i,j,0)))&
               /4.*.25
            auxz=W2(i,j,1)*((aer1(i,j,1)+aer1(i,j,0))+aer0(1)+aer0(0))/2.*.5

            if (W2(i,j,0).gt.0) then
               aeraux=-((auxx+auxy)+2.*auxz)*dt1/dx1
            else
               !se refleja un 25 % de los aerosoles que caen
               aeraux=-((auxx+auxy)+.25*2.*auxz)*dt1/dx1
            endif

            !agregamos un termino de turbulencia para los aerosoles
            !a nivel de piso
            turbu=cks/dx1*.25*(abs(U2(i,j,1))+abs(V2(i,j,1))+2.*abs(W2(i,j,1)))
            lapla=((aer1(i+1,j,0)+aer1(i-1,j,0))+(aer1(i,j+1,0)+aer1(i,j-1,0)+aer1(i,j,1)))-5.*aer1(i,j,0)
            lapla=lapla+(aer0(1)-aer0(0))

            aer2(i,j,0)=aeraux+aer1(i,j,0)+turbu*lapla

            aer2(i,j,nz1)=aer2(i,j,nz1-1)
400   continue

      !contornos laterales
      do 410 k=1,nz1-1
         do 410 j=1,nx1

            Titaa2(0,j,k)=Titaa2(1,j,k)
            Titaa2(nx1+1,j,k)=Titaa2(nx1,j,k)
            Titaa2(j,0,k)=Titaa2(j,1,k)
            Titaa2(j,nx1+1,k)=Titaa2(j,nx1,k)
            Qvap2(0,j,k)=0.
            Qvap2(nx1+1,j,k)=0.
            Qvap2(j,0,k)=0.
            Qvap2(j,nx1+1,k)=0.
            Qgot2(0,j,k)=Qgot2(1,j,k)
            Qgot2(nx1+1,j,k)=Qgot2(nx1,j,k)
            Qgot2(j,0,k)=Qgot2(j,1,k)
            Qgot2(j,nx1+1,k)=Qgot2(j,nx1,k)
            Qllu2(0,j,k)=0.
            Qllu2(nx1+1,j,k)=0.
            Qllu2(j,0,k)=0.
            Qllu2(j,nx1+1,k)=0.
            Qcri2(0,j,k)=Qcri2(1,j,k)
            Qcri2(nx1+1,j,k)=Qcri2(nx1,j,k)
            Qcri2(j,0,k)=Qcri2(j,1,k)
            Qcri2(j,nx1+1,k)=Qcri2(j,nx1,k)
            Qnie2(0,j,k)=Qnie2(1,j,k)
            Qnie2(nx1+1,j,k)=Qnie2(nx1,j,k)
            Qnie2(j,0,k)=Qnie2(j,1,k)
            Qnie2(j,nx1+1,k)=Qnie2(j,nx1,k)
            Qgra2(0,j,k)=0.
            Qgra2(nx1+1,j,k)=0.
            Qgra2(j,0,k)=0.
            Qgra2(j,nx1+1,k)=0.
            aer2(0,j,k)=aer2(1,j,k)
            aer2(nx1+1,j,k)=aer2(nx1,j,k)
            aer2(j,0,k)=aer2(j,1,k)
            aer2(j,nx1+1,k)=aer2(j,nx1,k)

410   continue

!############### calculo de la velocidad y la presion ################
      call velpre

!#####################################################################
!**  contornos, redefiniciones y filtros
!*** modificada las condiciones en el piso
      !Redefinicion
      do 155 i=1,nx1
         do 155 j=1,nx1

            k=0

            Titaa1(i,j,k)=pro3*Titaa2(i,j,k)+&
               pro4*((Titaa2(i+1,j,k)+Titaa2(i-1,j,k))+(Titaa2(i,j+1,k)+Titaa2(i,j-1,k)))

            if (abs(Titaa1(i,j,k)).lt.1e-10) Titaa1(i,j,k)=0

            Qvap1(i,j,k)=pro3*Qvap2(i,j,k)+&
               pro4*((Qvap2(i+1,j,k)+Qvap2(i-1,j,k))+(Qvap2(i,j+1,k)+Qvap2(i,j-1,k)))


            if (abs(Qvap1(i,j,k)).lt.1e-10) Qvap1(i,j,k)=0

            Qgot1(i,j,k)=pro3*Qgot2(i,j,k)+&
               pro4*((Qgot2(i+1,j,k)+Qgot2(i-1,j,k))+(Qgot2(i,j+1,k)+Qgot2(i,j-1,k)))

            if (Qgot1(i,j,k).lt.1e-10) Qgot1(i,j,k)=0

            Qllu1(i,j,k)=Qllu2(i,j,k)

            if (Qllu1(i,j,k).lt.1e-10) Qllu1(i,j,k)=0

            Qcri1(i,j,k)=pro3*Qcri2(i,j,k)+&
               pro4*((Qcri2(i+1,j,k)+Qcri2(i-1,j,k))+(Qcri2(i,j+1,k)+Qcri2(i,j-1,k)))

            if (Qcri1(i,j,k).lt.1e-10) Qcri1(i,j,k)=0

            Qnie1(i,j,k)=pro3*Qnie2(i,j,k)+pro4*((Qnie2(i+1,j,k)&
               +Qnie2(i-1,j,k))+(Qnie2(i,j+1,k)+Qnie2(i,j-1,k)))

            if (Qnie1(i,j,k).lt.1e-10) Qnie1(i,j,k)=0

            Qgra1(i,j,k)=Qgra2(i,j,k)

            if (Qgra1(i,j,k).lt.1e-10) Qgra1(i,j,k)=0

            aer1(i,j,k)=pro3*aer2(i,j,k)+pro4*((aer2(i+1,j,k)+aer2(i-1,j,k))+(aer2(i,j+1,k)+aer2(i,j-1,k)))

            !correccion cambiando la absorcion de aerosoles
            if (Qllu1(i,j,0).gt.1e-6 .and. W2(i,j,1).lt.-.1) then
               aeraux=(-W2(i,j,1)/1.)*(Qllu1(i,j,0)/1e-3)*.02*(dt1/5.)
               aer1(i,j,k)=aer1(i,j,k)-(aer1(i,j,k)+aer0(k))*aeraux
             endif

            if (abs(aer1(i,j,k)).lt.1e-10) aer1(i,j,k)=0


            do 155 k=1,nz1-1

               Titaa1(i,j,k)=pro1*Titaa2(i,j,k)+&
                  pro2*(&
                  (Titaa2(i+1,j,k)+Titaa2(i-1,j,k))+&
                  (Titaa2(i,j+1,k)+Titaa2(i,j-1,k))+&
                  Titaa2(i,j,k+1)+Titaa2(i,j,k-1))

               if (abs(Titaa1(i,j,k)).lt.1e-10) Titaa1(i,j,k)=0

               Qvap1(i,j,k)=pro1*Qvap2(i,j,k)+&
                  pro2*((Qvap2(i+1,j,k)+Qvap2(i-1,j,k))+&
                  (Qvap2(i,j+1,k)+Qvap2(i,j-1,k))+&
                  Qvap2(i,j,k+1)+Qvap2(i,j,k-1))

               if (abs(Qvap1(i,j,k)).lt.1e-10) Qvap1(i,j,k)=0

               Qgot1(i,j,k)=pro1*Qgot2(i,j,k)+&
                  pro2*((Qgot2(i+1,j,k)+Qgot2(i-1,j,k))+&
                  (Qgot2(i,j+1,k)+Qgot2(i,j-1,k))+&
                  Qgot2(i,j,k+1)+Qgot2(i,j,k-1))

               if (Qgot1(i,j,k).lt.1e-10) Qgot1(i,j,k)=0

               Qllu1(i,j,k)=pro1*Qllu2(i,j,k)+&
                  pro2*((Qllu2(i+1,j,k)+Qllu2(i-1,j,k))+&
                  (Qllu2(i,j+1,k)+Qllu2(i,j-1,k))+Qllu2(i,j,k+1)+Qllu2(i,j,k-1))

               if (Qllu1(i,j,k).lt.1e-10) Qllu1(i,j,k)=0

               Qcri1(i,j,k)=pro1*Qcri2(i,j,k)+&
                  pro2*((Qcri2(i+1,j,k)+Qcri2(i-1,j,k))+&
                  (Qcri2(i,j+1,k)+Qcri2(i,j-1,k))+Qcri2(i,j,k+1)+Qcri2(i,j,k-1))

               if (Qcri1(i,j,k).lt.1e-10) Qcri1(i,j,k)=0

               Qnie1(i,j,k)=pro1*Qnie2(i,j,k)+&
                  pro2*((Qnie2(i+1,j,k)+Qnie2(i-1,j,k))+&
                  (Qnie2(i,j+1,k)+Qnie2(i,j-1,k))+&
                  Qnie2(i,j,k+1)+Qnie2(i,j,k-1))

               if (Qnie1(i,j,k).lt.1e-10) Qnie1(i,j,k)=0

               Qgra1(i,j,k)=pro1*Qgra2(i,j,k)+&
                  pro2*((Qgra2(i+1,j,k)+Qgra2(i-1,j,k))+&
                  (Qgra2(i,j+1,k)+Qgra2(i,j-1,k))+&
                  Qgra2(i,j,k+1)+Qgra2(i,j,k-1))

               if (Qgra1(i,j,k).lt.1e-10) Qgra1(i,j,k)=0

               aer1(i,j,k)=pro1*aer2(i,j,k)+&
                  pro2*((aer2(i+1,j,k)+aer2(i-1,j,k))+&
                  (aer2(i,j+1,k)+aer2(i,j-1,k))+aer2(i,j,k+1)+aer2(i,j,k-1))


               if (abs(aer1(i,j,k)).lt.1e-10) aer1(i,j,k)=0
155   continue

      !contornos en el piso y en el techo
      do 420 i=1,nx1
         do 420 j=1,nx1
            Titaa1(i,j,0)=Titaa1(i,j,0)
            if (Titaa1(i,j,0).gt.0.5) Titaa1(i,j,0)=.5
            if (-Titaa1(i,j,0).gt.0.5) Titaa1(i,j,0)=-.5

            Titaa1(i,j,nz1)=Titaa1(i,j,nz1-1)

            !corregido para el vapor
            if (Qvap1(i,j,0).gt.Qvap0(0)*.5) then
               Qvap1(i,j,0)=.8*Qvap0(0)
            endif
            if (-Qvap1(i,j,0).gt.Qvap0(0)*.5) then
               Qvap1(i,j,0)=-.8*Qvap0(0)
            endif

            Qvap1(i,j,nz1)=Qvap1(i,j,nz1-1)
            Qgot1(i,j,0)=0.
            Qgot1(i,j,nz1)=Qgot1(i,j,nz1-1)
            Qllu1(i,j,0)=Qllu1(i,j,0)
            Qllu1(i,j,nz1)=Qllu1(i,j,nz1-1)
            Qcri1(i,j,0)=0.
            Qcri1(i,j,nz1)=Qcri1(i,j,nz1-1)

            Qnie1(i,j,0)=0.
            Qnie1(i,j,nz1)=Qnie1(i,j,nz1-1)

            Qgra1(i,j,0)=Qgra1(i,j,0)
            Qgra1(i,j,nz1)=Qgra1(i,j,nz1-1)

            !corregido para los aerosoles
            if (-aer1(i,j,0).gt.0.8*aer0(0)) then
               aer1(i,j,0)=-.8*aer0(0)
            endif
            aer1(i,j,nz1)=aer1(i,j,nz1-1)
420   continue

      !contornos laterales
      do 430 k=1,nz1-1
         do 430 j=1,nx1

            Titaa1(0,j,k)=Titaa1(1,j,k)
            Titaa1(nx1+1,j,k)=Titaa1(nx1,j,k)
            Titaa1(j,0,k)=Titaa1(j,1,k)
            Titaa1(j,nx1+1,k)=Titaa1(j,nx1,k)
            Qvap1(0,j,k)=0.
            Qvap1(nx1+1,j,k)=0.
            Qvap1(j,0,k)=0.
            Qvap1(j,nx1+1,k)=0.
            Qgot1(0,j,k)=Qgot1(1,j,k)
            Qgot1(nx1+1,j,k)=Qgot1(nx1,j,k)
            Qgot1(j,0,k)=Qgot1(j,1,k)
            Qgot1(j,nx1+1,k)=Qgot1(j,nx1,k)
            Qllu1(0,j,k)=0.
            Qllu1(nx1+1,j,k)=0.
            Qllu1(j,0,k)=0.
            Qllu1(j,nx1+1,k)=0.
            Qcri1(0,j,k)=Qcri1(1,j,k)
            Qcri1(nx1+1,j,k)=Qcri1(nx1,j,k)
            Qcri1(j,0,k)=Qcri1(j,1,k)
            Qcri1(j,nx1+1,k)=Qcri1(j,nx1,k)
            Qnie1(0,j,k)=Qnie1(1,j,k)
            Qnie1(nx1+1,j,k)=Qnie1(nx1,j,k)
            Qnie1(j,0,k)=Qnie1(j,1,k)
            Qnie1(j,nx1+1,k)=Qnie1(j,nx1,k)
            Qgra1(0,j,k)=0.
            Qgra1(nx1+1,j,k)=0.
            Qgra1(j,0,k)=0.
            Qgra1(j,nx1+1,k)=0.

            aer1(0,j,k)=aer1(1,j,k)
            aer1(nx1+1,j,k)=aer1(nx1,j,k)
            aer1(j,0,k)=aer1(j,1,k)
            aer1(j,nx1+1,k)=aer1(j,nx1,k)

430   continue

      !filtro para Titaa1 Qvap1
      call filtro(Titaa1,.01,.01,.02)

      !correccion de negativos para el vapor
      do 26 i=0,nx1+1
         do 26 j=0,nx1+1
            do 26 k=0,nz1
               if(Qvap1(i,j,k)+Qvap0(k).lt.0) then
                  Qvap1(i,j,k)=-Qvap0(k)
                  write(*,*) 'vapp',i,j,k,Qvap1(i,j,k)
               endif
26    continue

!############################# Prints ################################

      lll=tt
      ener=ener1+ener2+ener3+ener4+ener5
      write(*,*) 'fin',lll,Qcri1(16,16,14),Qcri2(16,16,14)
      write(*,*) W1(20,22,30),W1(19,22,30),W1(21,22,30),W1(20,21,30),W1(20,23,30)
      write(*,*) W1(19,21,30),W1(21,21,30),W1(21,23,30),W1(19,23,30),W1(21,22,30)
      write(*,*) Qcri1(20,22,27),Qnie1(20,22,27),Qvap1(20,22,27),Qvap0(27),W1(20,22,27)
      write(*,*) Qcri1(20,22,28),Qnie1(20,22,28),Qvap1(20,22,28),Qvap0(28),W1(20,22,28)
      write(*,*) Qcri1(20,22,29),Qnie1(20,22,29),Qvap1(20,22,29),Qvap0(29),W1(20,22,29)
      write(*,*) Qcri1(20,22,30),Qnie1(20,22,30),Qvap1(20,22,30),Qvap0(30),W1(20,22,30)
      write(*,*) Qcri1(20,22,31),Qnie1(20,22,31),Qvap1(20,22,31),Qvap0(31),W1(20,22,31)
      write(*,*) Qcri1(20,22,32),Qnie1(20,22,32),Qvap1(20,22,32),Qvap0(32),W1(20,22,32)
      write(*,*) Qcri1(20,22,33),Qnie1(20,22,33),Qvap1(20,22,33),Qvap0(33),W1(20,22,33)

!     Evolucion para puntos seleccionados
      write(13,44) U2(16,16,1),U2(17,17,1),U2(16,19,1),U2(17,19,1),U2(19,16,1),U2(19,17,1)
      write(14,44) V2(16,16,1),V2(17,17,1),V2(16,19,1),V2(17,19,1),V2(19,16,1),V2(19,17,1)
      write(15,44) W2(16,16,1),W2(17,17,1),W2(16,19,1),W2(17,19,1),W2(19,16,1),W2(19,17,1)
      write(16,44) Qvap1(16,16,1),Qvap1(17,17,1),Qvap1(16,19,1),Qvap1(17,19,1),Qvap1(19,16,1),Qvap1(19,17,1)
      write(17,44) Qgot1(16,16,1),Qgot1(17,17,1),Qgot1(16,19,1),Qgot1(17,19,1),Qgot1(19,16,1),Qgot1(19,17,1)
      write(18,44) aer1(16,16,1),aer1(17,17,1),aer1(16,19,1),aer1(17,19,1),aer1(19,16,1),aer1(19,17,1)
      write(19,44) Pres2(16,16,1),Pres2(17,17,1),Pres2(16,19,1),Pres2(17,19,1),Pres2(19,16,1),Pres2(19,17,1)
      write(20,44) Titaa1(16,16,1),Titaa1(17,17,1),Titaa1(16,19,1),Titaa1(17,19,1),Titaa1(19,16,1),Titaa1(19,17,1)

!     este es el unico que interesa
      aux1=0.
      aux2=0.
      aux3=0.
      aux4=0.
      do 810 j=nx1/2-1,nx1/2+2
         do 810 i=1,4
            aux1=aux1+aer1(i+nx1/2-10-posx(tt),j-posy(tt),0)
            aux2=aux2+aer1(i+nx1/2-2-posx(tt),j-posy(tt),0)
            aux3=aux3+aer1(i+nx1/2+6-posx(tt),j-posy(tt),0)
            aux4=aux4+aer1(j-posx(tt),i+nx1/2-10-posy(tt),0)
810   continue

      write(31,*)  tt,totnuc,totmic


!#####################################################################
!grabacion normal de las diferentes cantidades
      if (tt/nint(lte/dt1)*nint(lte/dt1).eq.tt) then
         call estad03_init()
         !desplazamiento de la nube
         tte=tte+1
         call posnub02_init()
         call corrinu2_init()
         write(32,*) tte,posx(tte),posy(tte),Xnub(tte),Ynub(tte),posxx,posyy
      endif

      if (tt/nint(ltg/dt1)*nint(ltg/dt1).eq.tt) then
         t1=t1+1
         
!############################ Grabación 2D ###########################
         !call graba231(k, t1, W2, Titaa1, Qvap1, Qllu1, Qgra1, aer1, Qvap0, aer0, tie, bre)
         
!############################ Grabación 3D ###########################
         write(*,*) 'graba320'
         call graba320(U1, V1, W1, Titaa1, Pres1, Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1, t1, tie, bre)
      endif

!######################### Grabación Backup ##########################
      if (tt/nint(ltb/dt1)*nint(ltb/dt1).eq.tt) then
         call graba120(Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV,&
         U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
         Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,aer1,aer2,Fcalo,&
         Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn)
      endif

      write(*,*) ""
      write(*,*) '----Tiempo transcurrido:',tt,'de',lt1,'----'
      write(*,*) ""
1  continue


!############### Fin de la evolucion temporal ########################
!#####################################################################

   close(13)
   close(14)
   close(15)
   close(16)
   close(17)
   close(18)
   close(19)
   close(20)
   close(30)
   close(31)
   close(32)
   close(33)

3000 format(a10,6i4)
4003 format(7g17.9)
4004 format(6g17.9)
4005 format(4g17.9)
4006 format(4g17.9)
4007 format(2g17.9)

44 format(6g16.8)

end