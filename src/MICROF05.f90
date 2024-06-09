!23456789*123456789*123456789*123456789*123456789*123456789*123456789*
!     Esta subrutina resuelve lo que le pasa a las gotitas,las gotas,
!     los cristales, la nieve, el granizo  y al vapor entre si.
!     Si hay nucleacion de gotitas no hay condensacion de agua liquida
!     Se ha modificado la parte de evaporacion teniendo en cuenta que
!     se debe tender a alcanzar el nivel de saturacion (25/6/99)
!     revisado 24/09/00

!$$
subroutine microfis(els,ess,Lvl,Lvs,Lsl,T,Dv,Eaccn,Eaucn,Eacng,Lsl00,Fcal,n,&
   qvapaux,qgotaux,qlluaux,qcriaux,qnieaux,qgraaux,Naer,daer2,nu,yy)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE microf05
   implicit none

   real, intent(in) :: els,ess,Lvl,Lvs,Lsl,T,Dv,Eaccn,Eaucn,Eacng,Lsl00,Naer,nu
   real, intent(inout) :: Fcal,daer2
   real(8), intent(inout) :: qvapaux,qgotaux,qlluaux,qcriaux,qnieaux,qgraaux
   integer, intent(in) :: n,yy

!*    Parametros comunes
   Nsc=nu/Dv

!*    parametros de las distribuciones y variables relacionadas
!     Gotitas
   Rgot=(qgotaux/cteqgot)**(1/6.)
   Ngot=2./27.*N0got*Rgot**3.

!     Lluvia
   Rllu=(qlluaux/cteqllu)**(1/4.)
   Nllu=N0llu*Rllu
   Vtm=Av(2*n)*Rllu**.8
   Nre=2.*Vtm*Rllu/nu
   A=Nre**(1/2.)*Nsc**(1/3.)
   if (A .le. 1.4) then
      fventl=1.+.108*A**2.
   else
      fventl=.78+.308*A
   endif

!     Cristales
   Rcri=5e-5
   if(T.lt.T0) Rcri=Rcri-4e-5*(T0-T)/40.
   if(Rcri.lt.1e-5) Rcri=1e-5
   Ncri=qcriaux/(pi*Rcri**3./10*rhocri)

!     Nieve
   Rnie=(qnieaux/cteqnie)**(1/4.)
   Nnie=N0nie*Rnie
   Nre=2.*Vtnie(2*n)*Rnie/nu
   if (T.le.T0) then
      fventn=.78+.308*Nsc**(1./3.)*Nre**.5
   else
      fventn=.86+.28*Nsc**(1./3.)*Nre**.5
   endif

!     Granizo
   Rgra=(qgraaux/cteqgra)**(1/4.)
   Ngra=N0gra*Rgra
   Vtgra=Vtgra0(2*n)*Rgra**.8
   Nre=2.*Vtgra*Rgra/nu
   fventgs=.78+.308*Nsc**(1./3.)*Nre**.5
   fventgl=.94+.308*Nsc**(1./3.)*Nre**.5

   s=0
   colilc=0.
   coliln=0.
   Fcal=0.
   Naux=Naer*1d6
!**********************************************************
!**    calculo de los distintos procesos para los hidrometeoros
!**********************************************************

!**   Evaporacion - Condensacion

   Qvls=els/Rv/T
   Qvls0=elvs0/Rv/T0
   Qvss=ess/Rv/T
   qauxl=qvapaux-Qvls
   qauxl0=qvapaux-Qvls0
   qauxs=qvapaux-Qvss

!*    Gotitas

   if(Rgot.gt.0 .and. yy.eq.0) then
      coevgot=-4.*pi*Dv*Rgot*Ngot*qauxl
!     para que no de cantidades negativas de gotas evaporando
      if (coevgot.gt. .95*qgotaux/dt2) then
         coevgot=qgotaux/dt2*.95
         s=1
      endif
!     para el vapor no se pase del de saturacion
      if (abs(coevgot).gt.abs(.95*qauxl/dt2)) coevgot=-.95*qauxl/dt2
   else
      coevgot=0.
   endif

!*    Lluvia

   if(Rllu.gt.0 .and. yy.eq.0) then
      coevllu=-4.*pi*Dv*qauxl*Rllu*Nllu*fventl
!     para que no de cantidades negativas de gotas evaporando
      if (coevllu.gt. .95*qlluaux/dt2) coevllu=qlluaux/dt2*.95
!     para el vapor no se pase del de saturacion
      if (abs(coevllu).gt.abs(.95*qauxl/dt2)) coevllu=-.95*qauxl/dt2
   else
      coevllu=0.
   endif

!*    Cristales

   If (qcriaux.gt.0 .and. T.lt.T0) then

      coevcri=-8.*Dv*qauxs*Rcri*Ncri

!     para que no de cantidades negativas de cristales evaporando
      if (coevcri.gt. .95*qcriaux/dt2) coevcri=qcriaux/dt2*.95
!     para el vapor no se pase del de saturacion
      if (abs(coevcri).gt.abs(.95*qauxs/dt2)) coevcri=-.95*qauxs/dt2
   else
      coevcri=0.
   endif

!$$
!*    Nieve

   if(Rnie.gt.0) then
      if (T.lt.T0) then

         coevnie=-8.*Dv*qauxs*Rnie*Nnie*fventn
      else
         coevnie=-8.*Dv*qauxl0*Rnie*Nnie*fventn
      endif

!     para que no de cantidades negativas de nieve evaporando
      if (coevnie.gt. .95*qnieaux/dt2) coevnie=qnieaux/dt2*.95
!     para el vapor no se pase del de saturacion
      if (abs(coevnie).gt.abs(.95*qauxs/dt2)) coevnie=-.95*qauxs/dt2
   else
      coevnie=0.
   endif

!*    para granizos (sup crecimiento humedo, se verifica en Tgra)
   if (Rgra.gt.0) then
      coevgra=-4.*pi*Dv*qauxl0*Rgra*Ngra*fventgl
!     para el vapor no se pase del de saturacion
      if (abs(coevgra).gt.abs(.95*qauxl0/dt2)) coevgra=-.95*qauxl0/dt2
   else
      coevgra=0.
   endif

!*    para aerosoles
   if (coevgot.gt.0 .and. yy.eq.0) then
      aux=max(Rgot,Rgotmin)
      libaer=coevgot/(4./3.*pi*rhow*aux**3.)
   else
      libaer=0.
   endif
   if (coevcri.gt.0) then
      libaer=libaer+coevcri/(pi*rhocri*Rcri**3./10.)*.5 !mod(2/4/00)
   endif
!**********************************************************

!**   Autoconversion

!*    gotitas a lluvia

   if (qgotaux.gt.2e-3) then !(mod 7/11/99)
      liqconv=-1e-3*qgotaux
   else
      liqconv=0.
   endif

!*    cristales a nieve

   if (T.lt.T0) then
      aux=20./3.*Eaucn*5e-3*.25*6.
      hieconv=-aux/Rcri/rhocri*qcriaux**2.

      if (qcriaux.gt.5e-4) then
         hieconv=hieconv-2.*qcriaux**2.
      endif
      if (Ncri.gt.100.) then
         aux=Acri*exp(Bcri*(T0-T))
         if (aux.lt.100.) aux=100.
         hieconv=hieconv-1e-6*qcriaux*Ncri/aux
      endif

      aux=(qcriaux-coevcri*dt2)*.9
      if (-hieconv.gt.aux) hieconv=-aux

   else
      hieconv=0.
   endif

   if (qauxs.lt.0 .and.T.le.T0 .and.qcriaux.lt.5e-4) then
      hieconv=hieconv+1e-2*qnieaux
   endif

   aux=exp(.09*(T-T0))
   if(aux.lt.0.2) aux=.2
   if (qnieaux .gt. 5e-4) then
      aux=aux*qnieaux/5e-4
      nieconv=-1e-2*qnieaux*aux
   else
      nieconv=0.
   endif


!**********************************************************

!**   Acrecion

!*    gotitas por lluvia
   if (Rgot.gt.0 .and. Rllu.gt.0 .and. s.eq.0) then
      acgollu=-gam3p8*qgotaux*Vtm*pi*Rllu**2.*Nllu*Efcol
!     para que no de cantidades negativas de gotas siendo colectadas
      if (acgollu.gt.qgotaux/dt2*.98) acgollu=qgotaux/dt2*.98
   else
      acgollu=0.
   endif

!*    gotitas por nieve
   if (Rgot.gt.0 .and. Rnie.gt.0 .and. s.eq.0) then
      acgonie=-2.*qgotaux*Vtnie(2*n)*pi*Rnie**2.*Nnie*Efcolgn
!     para que no de cantidades negativas de gotas siendo colectadas
      if (acgonie.gt.qgotaux/dt2*.98) acgonie=qgotaux/dt2*.98
   else
      acgonie=0.
   endif

!*    cristales por nieve
   if (Ncri.gt.0 .and. Rnie.gt.0 .and. T.lt.T0) then
      accrnie=-2.*qcriaux*Vtnie(2*n)*pi*(Rnie+Rcri)**2.*Nnie*Eaccn !mod 7/2/2000)
!     para que no de cantidades negativas de cristales siendo colectados
      if (accrnie.gt.qcriaux/dt2*.98) accrnie=qcriaux/dt2*.98
   else
      accrnie=0.
   endif

!*    gotitas por granizo
   if (Rgot.gt.0 .and. Rgra.gt.0 .and. s.eq.0) then
      acgogra=-gam3p8*qgotaux*Vtgra*pi*Rgra**2.*Ngra*Efcol
!     para que no de cantidades negativas de gotas siendo colectadas
      if (acgogra.gt.qgotaux/dt2*.98) acgogra=qgotaux/dt2*.98
   else
      acgogra=0.
   endif

!*    lluvia por granizo
   if (Rllu.gt.0 .and. Rgra.gt.0 .and. s.eq.0 .and.Vtgra.gt.Vtm) then
      acllgra=-pi*Efcol*Ngra*qlluaux*&
         (gam3p8*(Vtgra*Rgra**2.-Vtm*Rllu**2.)+&
         2.*gam2p8*Rgra*Rllu*(Vtgra-Vtm)+&
         2.*gam1p8*(Vtgra*Rllu**2.-Vtm*Rgra**2.))
!     para que no de cantidades negativas de lluvia siendo colectada
      if (acllgra.gt.qlluaux/dt2*.98) acllgra=qlluaux/dt2*.98
   else
      acllgra=0.
   endif

!*    cristales por granizo (sup aqui crecimiento humedo, si es seco se corrige en Tgra)
   if (Ncri.gt.0 .and. Rgra.gt.0 .and. s.eq.0 .and. T.lt.T0) then
      accrgra=-gam3p8*qcriaux*Vtgra*pi*Rgra**2.*Ngra*Efcol
!     para que no de cantidades negativas de critales siendo colectados
      if (accrgra.gt.qcriaux/dt2*.98) accrgra=qcriaux/dt2*.98
   else
      accrgra=0.
   endif

!*    nieve por granizo
   if (Rnie.gt.0 .and. Rgra.gt.0 .and. s.eq.0) then
      ccnigra=gam3p8*Vtgra*pi*Rgra**2.*Ngra*Nnie
      acnigra=-ccnigra*qnieaux/Nnie*Eacng
!     para que no de cantidades negativas de nieve siendo colectados
      if (acnigra.gt.qnieaux/dt2*.98) acnigra=qnieaux/dt2*.98
   else
      ccnigra=0.
      acnigra=0.
   endif

!**********************************************************
!**   Congelacion-Fusion de gotas y cristales, mod(14/4/99)

!*    congelacion homogenea de gotitas
   if (T.lt.238. .and. Rgot.gt.0) then
      aux=exp((238.-T)/8.)-1.
      if (T.lt.233) aux=exp(5./8.)-1.
      cfgotcri=-aux*qgotaux
   else
      cfgotcri=0.
   endif

!*    congelacion homogenea de lluvia
   if (T.lt.T0 .and. Rllu.gt.0) then
      aux=-6.7e-8*(exp(.66*(T0-T))-1.)
      if (aux.lt.-.9) aux=-.9
      cfllugra=aux*qlluaux
   else
      cfllugra=0.
   endif

   if (T.ge.T0 .and. Rcri.gt.0) then
      cfgotcri=.95*qcriaux
   else
      cfgotcri=0.
   endif

!*    fusion de la nieve
   if (Rnie.gt.0) then
      cfln1=-acgonie
      if (T.le.T0) then
         cfln2=-coevnie*Lvs/Lsl
         if(T.lt.T0-10.) cfln2=0. !mod 11/12/99
      else
         cfln2=-coevnie*Lvl0/Lsl00+&
            (8.*Kair*fventl*Rnie*Nnie+Cwl*cfln1)*(T-T0)/Lsl00
      endif
      if (cfln2.lt.0) cfln2=0.
      cfllunie=cfln1+cfln2

!     para que no de cantidades negativas nieve fundiendose
      if (cfln2.gt.qnieaux/dt2-coevnie) then
         cfln2=(qnieaux/dt2-coevnie)*.9
         cfllunie=cfln1+cfln2
      endif

   else
      cfllunie=0.
      cfln1=0.
      cfln2=0.
   endif

   if (cfllunie.lt.0) then
      stop
   endif

!**   colision entre lluvia y nieve o cristales que forman granizos (o nieve)

!*    lluvia con cristales     (mod 22/6/99)
   if(T.lt.T0.and.Rllu.gt.5e-5.and.Ncri.gt.0) then
      colilc=gam3p8*Vtm*pi*Rllu**2.*Nllu*Ncri*Efcol    !numero de colisiones
      colilc=min(colilc,Nllu*.3,Ncri*.3)
      congagua=colilc*(qcriaux/Ncri+qlluaux/Nllu)
   else
      colilc=0.
      congagua=0.
   endif

!*    lluvia con nieve     (mod 10/12/99)
   if(T.lt.T0.and.(Rllu.gt.5e-5.and.Rnie.gt.0))then
      coliln=pi*Efcol*Nllu*Nnie*(Vtm*&
         (gam3p8*Rllu**2.+gam1p8*2.*Rnie**2.+2.*gam2p8*Rnie*Rllu)&
         -Vtnie(2*n)*2.*(Rllu**2.+Rnie**2.+Rnie*Rllu))
      coliln=min(abs(coliln),Nllu*.3,Nnie*.3)
      congagua=congagua+coliln*(qnieaux/Nnie+qlluaux/Nllu)
   else
      coliln=0.
   endif

!**********************************************************
!**   Coleccion de aerosoles (7/2/2000)

!*    coleccion por lluvia
   coaerllu=-gam3p8*Vtm*pi*Rllu**2.*Nllu*Naux*Efcaer    !numero de colisiones

!*    coleccion por nieve
   coaernie=-2.*Vtnie(2*n)*pi*Rnie**2.*Nnie*Naux*Efcaer    !numero de colisiones

!*    coleccion por granizo
   coaergra=-gam3p8*Vtgra*pi*Rgra**2.*Ngra*Naux*Efcaer    !numero de colisiones

!*    difusion a gotitas

   coaergot=-4.*pi*Dfaer*Naux*Rgot*Ngot

!*    difusion a cristales

   coaercri=-8.*Dfaer*Naux*Rcri*Ncri


!********************************************************************
!     Calculo de la temperatura del granizo y determinacion del tipo
!     de crecimiento, se usan los terminos de intercambio anteriores
!     Primero suponemos crecimiento humedo y calculamos alfagra.
!     Si alfagra>1 el crecimiento es seco.
!     Si alfagra<0 o Tg>T0 hay melting

   if (Rgra.gt.0) then
      Qvaux=Qvls0
      Tg=T0
      alfagra=1.
      fugra=0.
      agual=-acgogra-acllgra
      hielo=-accrgra-acnigra
      A1=4.*pi*Rgra*Ngra*Kair*T*fventgl
      B1=4.*pi*Rgra*Ngra*Kair*fventgl
      A2=4.*pi*Rgra*Ngra*Lvl0*Dv*qvapaux*fventgl
      B2b=4.*pi*Rgra*Ngra*Lvl0*Dv*fventgl
      A3=agual*(Cwl*(T-T0)+alfagra*(Lsl00+Cwi*T0))
      A3b=agual*(Cwl*(T-T0))
      B3=agual*alfagra*Cwi
      A4=hielo*Cwi*T
      B4=hielo*Cwi

      Aalfa=A1+A2+A3b+A4
      Balfa=(B1+B4)*Tg+B2b*Qvls0

      if (agual.gt.0) then
         alfagra=-(Aalfa-Balfa)/agual/(Lsl00+Cwi*(T0-Tg))
         if (alfagra.ge.0 .and. alfagra.le.1) then ! crecimiento humedo
            crecigra=0
            A3=agual*(Cwl*(T-T0)+alfagra*(Lsl00+Cwi*T0))
            B3=agual*alfagra*Cwi
            Q1=A1-B1*Tg
            Q2=A2-B2b*Qvls0
            Q3=A3-B3*Tg
            Q4=A4-B4*Tg
            Qt=Q1+Q2+Q3+Q4
         endif
         if (alfagra.lt.0) then                          !fusion
            crecigra=1
            alfagra=0.
            A3=agual*(Cwl*(T-T0)+alfagra*(Lsl00+Cwi*T0))
            B3=agual*alfagra*Cwi
            Q1=A1-B1*Tg
            Q2=A2-B2b*Qvls0
            Q3=A3-B3*Tg
            Q4=A4-B4*Tg
            Qt=Q1+Q2+Q3+Q4
            fugra=-Qt/Lsl00
!       para que no den cantidades de granizo fundiendo
            if (-fugra.gt.qgraaux/dt2) fugra=-qgraaux/dt2*.9
         endif
      endif

      if (alfagra.gt.1. .or. agual.eq.0) then
         alfagra=1.
         accrgra=accrgra/8.
         hielo=-accrgra-acnigra
         A2=4.*pi*Rgra*Ngra*Lvs0*Dv*qvapaux*fventgs
         B2b=4.*pi*Rgra*Ngra*Lvs0*Dv*fventgs
         A4=hielo*Cwi*T
         B4=hielo*Cwi
         Q1=A1-B1*Tg
         Q2=A2-B2b*Qvls0
         Q3=A3-B3*Tg
         Q4=A4-B4*Tg
         Qt=Q1+Q2+Q3+Q4
         BB=B1+B3+B4
         CC2=B2b/Rv
         CC3=B2b*esvs0/Rv**2.*Lvs0
         esaux=esvs0*exp(Lvs0/Rv*(1./T0-1./Tg))
         dQt=-BB+CC2*esaux/Tg**2.-CC3/Tg**3.
         Taux=Tg-Qt/dQt
         do i=1,10
            if (abs(Taux-Tg) .gt. .05) then
               Taux=Tg
               Qvaux=esaux/Rv/Taux
               Q1=A1-B1*Taux
               Q2=A2-B2b*Qvaux
               Q3=A3-B3*Taux
               Q4=A4-B4*Taux
               Qt=Q1+Q2+Q3+Q4
               dQt=-BB+CC2*esaux/Tg**2.-CC3/Tg**3.
               Tg=Taux-Qt/dQt
               esaux=esvs0*exp(Lvs0/Rv*(1./T0-1./Tg))
            endif
         end do

         if (Tg .le. T0)  then
            crecigra=2          ! crecimiento seco
            Qvaux=esvs0/Rv/Tg*exp(Lvs0/Rv*(1./T0-1./Tg))
            coevgra=-4.*pi*Dv*Rgra*Ngra*fventgs*&
               (qvapaux-Qvaux)
            Q2=-coevgra*Lvs0
            Qt=Q1+Q2+Q3+Q4
         else
            crecigra=3          ! hay fusion
            Tg=T0
            Q1=A1-B1*Tg
            Q2=A2-B2b*Qvls0
            Q3=A3-B3*Tg
            Q4=A4-B4*Tg
            Qt=Q1+Q2+Q3+Q4
            fugra=-Qt/Lsl00
!       para que no den cantidades de granizo fundiendo
            if (-fugra.gt.qgraaux/dt2*.8) fugra=-qgraaux/dt2*.8
         endif
      endif
      fugrallu=-agual*(1.-alfagra)+fugra
      Fcalgra=(fugrallu-congagua)*Cwl*(Tg-T)
   else
      fugrallu=0.
      Fcalgra=0.
   endif

!**********************************************************
!**   multiplicacion de cristales por colisiones
!*     colisiones entre granizo con nieve
   mucrgrni=ccnigra*Ncrgrni*pi*Rcri**3./10.*rhocri
   if (T.gt.T0) mucrgrni=0.

!*     colisiones entre granizo con granizos
   if (Rgra.gt.0 .and. T.gt.0) then
      aux=4.*pi*Rgra**2.*Vtgra*.25*Ngra**2.
      mucrgrgr=aux*Ncrgrgr*pi*Rcri**3./10.*rhocri
   else
      mucrgrgr=0.
   endif

!**********************************************************
!**********************************************************

!**   terminos de intercambio
!$$
   invapgot=coevgot
   invapllu=coevllu
   invapcri=coevcri
   invapnie=coevnie
   invapgra=coevgra
   ingotllu=liqconv+acgollu
   ingotcri=cfgotcri
   ingotnie=acgonie
   ingotgra=acgogra
   incrinie=hieconv+accrnie+mucrgrni*.8
   incrigra=accrgra+mucrgrni*.2+mucrgrgr
   inllunie=cfllunie
   inllugra=cfllugra-fugrallu
   inniegra=nieconv

   if (Rllu.gt.1e-4 .and. Ncri.gt.0 .and. T.lt.T0) then !mod 10/12/99
      incrigra=incrigra-colilc*qcriaux/Ncri
      inllugra=inllugra-colilc*qlluaux/Nllu
      Fcal=colilc*qlluaux/Nllu*Lsl*dt2
   endif
   if (Rllu.gt.5e-5.and.Rllu.le.1e-4.and.Ncri.gt.0.and.T.lt.T0)then !mod 10/12/99
      incrinie=incrinie-colilc*qcriaux/Ncri
      inllunie=inllunie-colilc*qlluaux/Nllu
      Fcal=colilc*qlluaux/Nllu*Lsl*dt2
   endif

   if (Rllu.gt.5e-5 .and. Rnie.gt.0.) then !mod 10/12/99
      inllugra=inllugra-coliln*qlluaux/Nllu
      inniegra=inniegra-coliln*qnieaux/Nnie
      Fcal=Fcal+coliln*qlluaux/Nllu*Lsl*dt2
   endif

!     calculo del numero de particulas
!$$
   Intvap=invapgot+invapllu+invapcri+invapnie+invapgra
   Intgot=-invapgot+ingotllu+ingotcri+ingotnie+ingotgra
   Intllu=-invapllu-ingotllu+inllunie+inllugra
   Intcri=-invapcri-ingotcri+incrinie+incrigra
   Intnie=-invapnie-ingotnie-inllunie-incrinie+inniegra
   Intgra=-invapgra-ingotgra-inllugra-incrigra-inniegra
   Intaer=libaer+coaergot+coaerllu+coaercri+coaernie+coaergra

   qvapaux1=qvapaux+Intvap*dt2
   qgotaux1=qgotaux+Intgot*dt2
   qlluaux1=qlluaux+Intllu*dt2
   qcriaux1=qcriaux+Intcri*dt2
   qnieaux1=qnieaux+Intnie*dt2
   qgraaux1=qgraaux+Intgra*dt2

!     Correccion de lluvia negativa
   if (qlluaux1.lt.0 .and. T.lt.253) then
      qgraaux1=qgraaux1+qlluaux1
      qlluaux1=0.
   endif

!     Correccion de gotas negativas
!$$
   if (qgotaux1.lt.0) then
      invapgot=invapgot*(1.-qgotaux1/(Intgot*dt2))
      ingotllu=ingotllu*(1.-qgotaux1/(Intgot*dt2))
      ingotcri=ingotcri*(1.-qgotaux1/(Intgot*dt2))
      ingotnie=ingotnie*(1.-qgotaux1/(Intgot*dt2))
      ingotgra=ingotgra*(1.-qgotaux1/(Intgot*dt2))

      Intvap=invapgot+invapllu+invapcri+invapnie+invapgra
      Intgot=-invapgot+ingotllu+ingotcri+ingotnie+ingotgra
      Intllu=-invapllu-ingotllu+inllunie+inllugra
      Intcri=-invapcri-ingotcri+incrinie+incrigra
      Intnie=-invapnie-ingotnie-inllunie-incrinie+inniegra
      Intgra=-invapgra-ingotgra-inllugra-incrigra-inniegra

      qvapaux1=qvapaux+Intvap*dt2
      qgotaux1=qgotaux+Intgot*dt2
      qlluaux1=qlluaux+Intllu*dt2
      qcriaux1=qcriaux+Intcri*dt2
      qnieaux1=qnieaux+Intnie*dt2
      qgraaux1=qgraaux+Intgra*dt2

   endif

!     prevencion de negativos
!$$
   if (qgotaux1.lt.0 .or. qlluaux1.lt.0 .or. qcriaux1.lt.0 .or.&
      qnieaux1.lt.0 .or. qgraaux1.lt.0) then

      if (qgotaux1.lt.-1e-10 .or. qlluaux1.lt.-1e-10 .or. &
         qcriaux1.lt.-1e-10 .or.&
         qnieaux1.lt.-1e-10 .or. qgraaux1.lt.-1e-10) then
         stop
      endif
      if (qgotaux1.lt.0) qgotaux1=0.
      if (qlluaux1.lt.0) qlluaux1=0.
      if (qcriaux1.lt.0) qcriaux1=0.
      if (qnieaux1.lt.0) qnieaux1=0.
      if (qgraaux1.lt.0) qgraaux1=0.

!	stop
   endif

   aux=Intvap+Intgot+Intllu+Intcri+Intnie+Intgra
   if (aux.gt.1e-9) then
      stop
   endif

   qvapaux=qvapaux1
   qgotaux=qgotaux1
   qlluaux=qlluaux1
   qcriaux=qcriaux1
   qnieaux=qnieaux1
   qgraaux=qgraaux1
   daer2=Intaer*dt2/1e6


!     calculo del calor obtenido por cambio de fase (por m^-3)
!$$
   Fcal=Fcal+(-(invapgot+invapllu)*Lvl-&
      ingotcri*Lsl-invapcri*Lvs+Fcalgra)*dt2

   if (invapnie.gt.0 .or. T.lt.T0-10.) then !mod11/12/99
      Fcal=Fcal-invapnie*Lvs*dt2
   endif

   if (T.lt.T0) Fcal=Fcal-inllunie*Cwl*(T-T0)*dt2
   return
end
