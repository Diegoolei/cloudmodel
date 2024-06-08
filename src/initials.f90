module velpre01
   integer :: t, i, j, k
   real :: dvx, dvy, dvz, diver, dprex, dprey, dprez, vel0, vel1, vel2,&
      vel3, presi, presix, presiy, presiz, facx, facy, facz, prom1, prom,&
      kkk, presprom, nnn
contains
   subroutine velpre01_init()
      USE dimen
      implicit none
      facx=.05
      facy=.05
      facz=.05

      prom=.3/6.*(dt3/.2)
      prom1=1.-prom*6.
      kkk=.01
      nnn=(nx1+2)**2.*(nz1+1)

   end subroutine velpre01_init
end module velpre01

module aeroana
!> Variables de aeroana
   USE dimen
   integer tt,t1,t2,n,m,l,i,j,k,lll,s,iT,tte, lvapneg,llluneg,lcrineg,&
      laerneg,lnieneg,lgraneg,yy
   real T,P, Dv,Lvl,Lvs,Lsl,Vis,Qvap,Qliq,densi,nu, Lsl00, Eaucn,Eaccn,&
      Eacng, Naer,dqgot,dqcri,daer,daer2, Fcal, elvs,esvs,e1,rl,rs,dden0z,&
      aux,aux1,aux2,aux3,aux4, cks,turbu,lapla
   real(8) qgotaux,qvapaux,qlluaux,qcriaux,qnieaux,qgraaux,aeraux, auxx,auxy,&
      auxz, Taux,Qvapneg,aerneg, ener,ener1,ener2,ener3,ener4,ener5,qv,qg,&
      daitot,vapt1,vapt2,vapt3,vapt4, gott1,gott2,gott3,gott4, aert1,aert2,&
      aert3,aert4,totnuc,totmic
   real aerdif(-3:nx1+3,-3:nx1+3,-3:nz1+3)

contains
   subroutine aeroana_init()
      ctur=.5
      lt2=nint(dt1/dt2)
      lt3=2*nint(dt1/dt3)
      cteturb=ctur/2.**.5
      cks=cteturb*2.
      dx2=2.*dx1
      dx8=8.*dx1
      dx12=12.*dx1
      AA=1./Eps-1.
      ikapa=1./Kapa
      cteqgot=160./3**6.*pi*rhow*N0got
      cteqllu=8.*pi*rhow*N0llu
      cteqnie=.6*pi*rhonie*N0nie
      cteqgra=8.*pi*rhogra*N0gra
      tte=0
   end subroutine aeroana_init

end module aeroana

module mode20
   implicit none

   real :: T, P, Dv, Lvl, Lvs, Lsl, Vis, Qvap, Qliq, densi, nu, Lsl00, Eaucn,&
      Eaccn, Eacng, Naer, dqgot, dqcri, daer, daer2, Fcal, elvs, esvs, e1, rl,&
      rs, dden0z, aux, aux1, aux2, aux3, aux4, posxx, posyy, zmed,cks, turbu, lapla

   real(8) :: qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, aeraux,&
      auxx, auxy, auxz, Taux, Qvapneg, aerneg, ener, ener1, ener2, ener3,&
      ener4, ener5, qv, qg, daitot, vapt1, vapt2, vapt3, vapt4, gott1, gott2,&
      gott3, gott4, aert1, aert2, aert3, aert4, totnuc, totmic, impx, impy,&
      Qagua, Qaguat

   real(8) :: Xnub(5000), Ynub(5000)
   integer :: posx(-3:5000), posy(-3:5000), spos

   character(len=3) :: file_number
   character(len=12) :: nombre

   integer :: tt, t1, t2, n, m, l, i, j, k, lll, s, iT, tte, lvapneg,&
      llluneg, lcrineg, laerneg, lnieneg, lgraneg, yy, laux1, laux2, &
      maux1, maux2, naux2, umax, umin, vmax, vmin, wmax, wmin, titamax,&
      titamin, qvapmax, qvapmin, qgotmax, qllumax, qcrimax, qniemax,&
      qgramax, aermax, lumax, mumax, numax,lumin, mumin, numin, lvmax,&
      mvmax, nvmax, lvmin, mvmin, nvmin, lwmax, mwmax, nwmax, lwmin,&
      mwmin, nwmin, ltitamax, mtitamax, ntitamax, ltitamin, mtitamin,&
      ntitamin, lqvapmax, mqvapmax, nqvapmax, lqvapmin, mqvapmin, nqvapmin,&
      lqgotmax,mqgotmax, nqgotmax, lqllumax, mqllumax, nqllumax, laermax,&
      maermax, naermax, lqcrimax, mqcrimax, nqcrimax, lqniemax, mqniemax,&
      nqniemax, lqgramax, mqgramax, nqgramax, qgottot, qllutot, qcritot,&
      qnietot, qgratot


contains

   subroutine mode20_init()
      USE cant01
      USE dimen
      USE const
      USE config

      call init_config()
      ini = 0                                   !inicio por vez primera= 0
      t1 = 0                                    !paso a inicio (si ini=0->t1=0)
      ltt = sim_time_minutes * 60. * 2.         !tiempo total de simulacion
      ltg = save_lapse_minutes * 60. * 2.       !tiempo de grabacion
      lte = 3. * 60. * 2.                       !tiempo de grabacion estadistica
      ltb = 3. * 60. * 2.                       !tiempo de backup
      ctur = 0.5

      pro1 = 1. - 2e-2 * (dt1 / 5.)
      pro2 = (1. - pro1) / 6.
      pro3 = 1. - 2e-2 * (dt1 / 5.)
      pro4 = (1. - pro1) / 4.

      lt1 = nint(ltt / dt1)
      lt2 = nint(dt1 / dt2)                  ! Proporción Física/Microfísica
      lt3 = 2 * nint(dt1 / dt3)
      cteturb = ctur / 2.**.5

      cks = cteturb * 2.

      dx2 = 2. * dx1
      dx8 = 8. * dx1
      dx12 = 12. * dx1
      AA = 1. / Eps - 1.
      ikapa = 1. / Kapa

      cteqgot = 160. / 3**6. * pi * rhow * N0got
      cteqllu = 8. * pi * rhow * N0llu
      cteqnie = 0.6 * pi * rhonie * N0nie
      cteqgra = 8. * pi * rhogra * N0gra

      spos = 0
      posxx = 0.
      posyy = 0.
      posx(0) = 0
      posy(0) = 0
      tte = 0
   end subroutine mode20_init

end module mode20

module estad03
   !     calculo de algunas cantidades de interes (2/03/2000)
contains
   subroutine estad03_init()
      USE mode20
      USE permic
      USE perdim
      USE config
      integer :: unit_number
      umax=0.
      lumax=0
      mumax=0
      numax=0
      vmax=0.
      lvmax=0
      mvmax=0
      nvmax=0
      wmax=0.
      lwmax=0
      mwmax=0
      nwmax=0
      titamax=0.
      ltitamax=0
      mtitamax=0
      ntitamax=0
      qvapmax=0.
      lqvapmax=0
      mqvapmax=0
      nqvapmax=0
      qgotmax=0.
      lqgotmax=0
      mqgotmax=0
      nqgotmax=0
      qllumax=0.
      lqllumax=0
      mqllumax=0
      nqllumax=0
      qcrimax=0.
      lqcrimax=0
      mqcrimax=0
      nqcrimax=0
      qniemax=0.
      lqniemax=0
      mqniemax=0
      nqniemax=0
      qgramax=0.
      lqgramax=0
      mqgramax=0
      nqgramax=0
      aermax=0.
      laermax=0
      maermax=0
      naermax=0
      umin=0.
      lumin=0
      mumin=0
      numin=0
      vmin=0.
      lvmin=0
      mvmin=0
      nvmin=0
      wmin=0.
      lwmin=0
      mwmin=0
      nwmin=0
      titamin=0.
      ltitamin=0
      mtitamin=0
      ntitamin=0
      qvapmin=0.
      lqvapmin=0
      mqvapmin=0
      nqvapmin=0
      qgottot=0.
      qllutot=0.
      qcritot=0.
      qnietot=0.
      qgratot=0.

      do 700 k=1,nz1
         do 700 i=1,nx1
            do 700 j=1,nx1

               if (umax.lt.U1(i,j,k)*100) then
                  umax=U1(i,j,k)*100
                  lumax=i
                  mumax=j
                  numax=k
               endif

               if (umin.gt.U1(i,j,k)*100) then
                  umin=U1(i,j,k)*100
                  lumin=i
                  mumin=j
                  numin=k
               endif

               if (vmax.lt.V1(i,j,k)*100) then
                  vmax=V1(i,j,k)*100
                  lvmax=i
                  mvmax=j
                  nvmax=k
               endif

               if (vmin.gt.V1(i,j,k)*100) then
                  vmin=V1(i,j,k)*100
                  lvmin=i
                  mvmin=j
                  nvmin=k
               endif

               if (wmax.lt.W1(i,j,k)*100) then
                  wmax=W1(i,j,k)*100
                  lwmax=i
                  mwmax=j
                  nwmax=k
               endif

               if (wmin.gt.W1(i,j,k)*100) then
                  wmin=W1(i,j,k)*100
                  lwmin=i
                  mwmin=j
                  nwmin=k
               endif

               if (titamax.lt.Titaa1(i,j,k)*1000) then
                  titamax=Titaa1(i,j,k)*1000
                  ltitamax=i
                  mtitamax=j
                  ntitamax=k
               endif

               if (titamin.gt.Titaa1(i,j,k)*1000) then
                  titamin=Titaa1(i,j,k)*1000
                  ltitamin=i
                  mtitamin=j
                  ntitamin=k
               endif

               if (qvapmax.lt.Qvap1(i,j,k)*1e6) then
                  qvapmax=Qvap1(i,j,k)*1e6
                  lqvapmax=i
                  mqvapmax=j
                  nqvapmax=k
               endif

               if (qvapmin.gt.Qvap1(i,j,k)*1e6) then
                  qvapmin=Qvap1(i,j,k)*1e6
                  lqvapmin=i
                  mqvapmin=j
                  nqvapmin=k
               endif

               if (qgotmax.lt.Qgot1(i,j,k)*1e6) then
                  qgotmax=Qgot1(i,j,k)*1e6
                  lqgotmax=i
                  mqgotmax=j
                  nqgotmax=k
               endif
               qgottot=qgottot+Qgot1(i,j,k)*1e6

               if (qllumax.lt.Qllu1(i,j,k)*1e6) then
                  qllumax=Qllu1(i,j,k)*1e6
                  lqllumax=i
                  mqllumax=j
                  nqllumax=k
               endif
               qllutot=qllutot+Qllu1(i,j,k)*1e6

               if (qcrimax.lt.Qcri1(i,j,k)*1e6) then
                  qcrimax=Qcri1(i,j,k)*1e6
                  lqcrimax=i
                  mqcrimax=j
                  nqcrimax=k
               endif
               qcritot=qcritot+Qcri1(i,j,k)*1e6

!$$
               if (qniemax.lt.Qnie1(i,j,k)*1e6) then
                  qniemax=Qnie1(i,j,k)*1e6
                  lqniemax=i
                  mqniemax=j
                  nqniemax=k
               endif
               qnietot=qnietot+Qnie1(i,j,k)*1e6

               if (qgramax.lt.Qgra1(i,j,k)*1e6) then
                  qgramax=Qgra1(i,j,k)*1e6
                  lqgramax=i
                  mqgramax=j
                  nqgramax=k
               endif
               qgratot=qgratot+Qgra1(i,j,k)*1e6

               if (aermax.lt.aer1(i,j,k)/1000) then
                  aermax=aer1(i,j,k)/1000
                  laermax=i
                  maermax=j
                  naermax=k
               endif
700   continue


      qgotmax=0.
      qcrimax=0.
      qniemax=0.
      do 719 i=-1,1
         do 719 j=-1,1
            do 719 k=-1,1
               qgotmax=qgotmax+1e5*Qgot1(lqgotmax+i,mqgotmax+j,nqgotmax+k)
               qcrimax=qcrimax+1e5*Qcri1(lqcrimax+i,mqcrimax+j,nqcrimax+k)
               qniemax=qniemax+1e5*Qnie1(lqniemax+i,mqniemax+j,nqniemax+k)

719   continue

      qgotmax=qgotmax/27.
      qcrimax=qcrimax/27.
      qniemax=qniemax/27.
      umax=umax/10
      umin=umin/10
      vmax=vmax/10
      vmin=vmin/10
      wmax=wmax/10
      wmin=wmin/10
      titamax=titamax/10
      titamin=titamin/10
      qvapmax=qvapmax/10
      qvapmin=qvapmin/10
      qgottot=qgottot/1000
      qllutot=qllutot/1000
      qcritot=qcritot/1000
      qnietot=qnietot/1000
      qgratot=qgratot/1000

      open(newunit=unit_number,file=output_directory//"esta", ACCESS="append")
      write(unit_number,710) umax,umin,vmax,vmin,wmax,wmin,titamax,titamin&
         ,qvapmax,qvapmin,qgotmax,qllumax,qcrimax,qniemax&
         ,qgramax,aermax&
         ,lumax,mumax,numax,lumin,mumin,numin&
         ,lvmax,mvmax,nvmax,lvmin,mvmin,nvmin&
         ,lwmax,mwmax,nwmax,lwmin,mwmin,nwmin&
         ,ltitamax,mtitamax,ntitamax,ltitamin,mtitamin&
         ,ntitamin,lqvapmax,mqvapmax,nqvapmax,lqvapmin&
         ,mqvapmin,nqvapmin,lqgotmax,mqgotmax,nqgotmax&
         ,lqllumax,mqllumax,nqllumax&
         ,lqcrimax,mqcrimax,nqcrimax&
         ,lqniemax,mqniemax,nqniemax&
         ,lqgramax,mqgramax,nqgramax&
         ,laermax,maermax,naermax

      open(newunit=unit_number,file=output_directory//"est", ACCESS="append")
      write(unit_number,715) qgottot,qllutot,qcritot,qnietot,qgratot
      close(unit_number)
710   format(16i5,48i4)
715   format(5i9)
   end subroutine estad03_init
end module estad03


module posnub02
contains
   subroutine posnub02_init()
      USE mode20
      USE lmngot
      USE lmncri
      USE permic
      USE estbas
      USE cant01
!     desplazamientos horizontales a partir de la velocidad media de la nube
!     calculo la altura media de la nube, la velocidad media de la nube
!     es tomada como la velocidad del aire sin perturbar a esa altura
!                  (1/3/2000)

      if (ngot(2).ge.1 .or. ncri(2).gt.1) then

         impx=0.
         impy=0.
         zmed=0.
         Qaguat=0.
         spos=1
         Xnub(tte)=Xnub(tte-1)
         Ynub(tte)=Ynub(tte-1)

         laux1=min(lgot(1),lcri(1))
         laux2=max(lgot(2),lcri(2))
         maux1=min(mgot(1),mcri(1))
         maux2=max(mgot(2),mcri(2))
         naux2=max(ngot(2),ncri(2))

         do 1300 k=1,naux2
            do 1300 i=laux1,laux2
               do 1300 j=maux1,maux2

                  Qagua=Qgot1(i,j,k)+Qcri1(i,j,k)+Qllu1(i,j,k)+&
                     Qnie1(i,j,k)+Qgra1(i,j,k)
                  Qagua=Qgot1(i,j,k)+Qllu1(i,j,k)+&
                     (Qcri1(i,j,k)+Qnie1(i,j,k)+Qgra1(i,j,k))/1000
                  zmed=zmed+k*Qagua
                  Qaguat=Qaguat+Qagua

1300     continue

         if (Qaguat.gt.1e-3) then
            zmed=zmed/Qaguat
            Xnub(tte)=Xnub(tte)+UU(nint(zmed))*lte
            Ynub(tte)=Ynub(tte)+VV(nint(zmed))*lte
         endif

      endif

   end subroutine posnub02_init
end module posnub02

module corrinu2
contains
   subroutine corrinu2_init()
      USE mode20
      USE perdim
      USE permic
      !     movimiento de la nube (4/01/99)
!     Redefine el valor de todas las variables (deberian
!      ser las que se graban solamente)

!     calculo de la posicion media de la nube, es decir de las gotitas
!     el centro esta inicialmente en nx1/2+.5, nx1/2+.5
!     posxx y posyy son siempre en modulo menores que dx1
!     En posx y posy se guarda para cada tt la posicion en
!     puntos de red

      if (spos .eq. 1) then
         posxx=Xnub(tte)
         posyy=Ynub(tte)
      else
         posxx=0.
         posyy=0.
      endif

      posx(tte)=posx(tte-1)
      posy(tte)=posy(tte-1)

!*    corrimiento en x

      if (posxx.gt.dx1) then
         posx(tte)=posx(tte)+1
         Xnub(tte)=Xnub(tte)-dx1

         do 1500 k=0,nz1+1
            do 1501 j=0,nx1+1
               do 1502 i=1,nx1+1
                  U1(i-1,j,k)=U1(i,j,k)
                  V1(i-1,j,k)=V1(i,j,k)
                  W1(i-1,j,k)=W1(i,j,k)
                  Pres1(i-1,j,k)=Pres1(i,j,k)
                  U2(i-1,j,k)=U2(i,j,k)
                  V2(i-1,j,k)=V2(i,j,k)
                  W2(i-1,j,k)=W2(i,j,k)
                  Pres2(i-1,j,k)=Pres2(i,j,k)
                  Titaa1(i-1,j,k)=Titaa1(i,j,k)
                  Qvap1(i-1,j,k)=Qvap1(i,j,k)
                  Qgot1(i-1,j,k)=Qgot1(i,j,k)
                  Qllu1(i-1,j,k)=Qllu1(i,j,k)
                  Qcri1(i-1,j,k)=Qcri1(i,j,k)
                  Aer1(i-1,j,k)=Aer1(i,j,k)
                  Fcalo(i-1,j,k)=Fcalo(i,j,k)
1502           continue
               i=nx1+1
               U1(i,j,k)=U1(i-1,j,k)
               V1(i,j,k)=V1(i-1,j,k)
               W1(i,j,k)=W1(i-1,j,k)
               Pres1(i,j,k)=Pres1(i-1,j,k)
               U2(i,j,k)=U2(i-1,j,k)
               V2(i,j,k)=V2(i-1,j,k)
               W2(i,j,k)=W2(i-1,j,k)
               Pres2(i,j,k)=Pres2(i-1,j,k)
               Titaa1(i,j,k)=Titaa1(i-1,j,k)
               Qvap1(i,j,k)=Qvap1(i-1,j,k)
               Qgot1(i,j,k)=Qgot1(i-1,j,k)
               Qllu1(i,j,k)=Qllu1(i-1,j,k)
               Qcri1(i,j,k)=Qcri1(i-1,j,k)
               Aer1(i,j,k)=Aer1(i-1,j,k)
               Fcalo(i,j,k)=0.
1501        continue
            i=nx1
            j=0
            U1(i,j,k)=(U1(i-1,j,k)+U1(i,j+1,k))/2.
            V1(i,j,k)=(V1(i-1,j,k)+V1(i,j+1,k))/2.
            W1(i,j,k)=(W1(i-1,j,k)+W1(i,j+1,k))/2.
            Pres1(i,j,k)=(Pres1(i-1,j,k)+Pres1(i,j+1,k))/2.
            U2(i,j,k)=(U2(i-1,j,k)+U2(i,j+1,k))/2.
            V2(i,j,k)=(V2(i-1,j,k)+V2(i,j+1,k))/2.
            W2(i,j,k)=(W2(i-1,j,k)+W2(i,j+1,k))/2.
            Pres2(i,j,k)=(Pres2(i-1,j,k)+Pres2(i,j+1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i-1,j,k)+Titaa1(i,j+1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i-1,j,k)+Qvap1(i,j+1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i-1,j,k)+Qgot1(i,j+1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i-1,j,k)+Qllu1(i,j+1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i-1,j,k)+Qcri1(i,j+1,k))/2.
            Aer1(i,j,k)=(Aer1(i-1,j,k)+Aer1(i,j+1,k))/2.
            Fcalo(i,j,k)=0.
            j=nx1+1
            U1(i,j,k)=(U1(i-1,j,k)+U1(i,j-1,k))/2.
            V1(i,j,k)=(V1(i-1,j,k)+V1(i,j-1,k))/2.
            W1(i,j,k)=(W1(i-1,j,k)+W1(i,j-1,k))/2.
            Pres1(i,j,k)=(Pres1(i-1,j,k)+Pres1(i,j-1,k))/2.
            U2(i,j,k)=(U2(i-1,j,k)+U2(i,j-1,k))/2.
            V2(i,j,k)=(V2(i-1,j,k)+V2(i,j-1,k))/2.
            W2(i,j,k)=(W2(i-1,j,k)+W2(i,j-1,k))/2.
            Pres2(i,j,k)=(Pres2(i-1,j,k)+Pres2(i,j-1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i-1,j,k)+Titaa1(i,j-1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i-1,j,k)+Qvap1(i,j-1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i-1,j,k)+Qgot1(i,j-1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i-1,j,k)+Qllu1(i,j-1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i-1,j,k)+Qcri1(i,j-1,k))/2.
            Aer1(i,j,k)=(Aer1(i-1,j,k)+Aer1(i,j-1,k))/2.
            Fcalo(i,j,k)=0.
1500     continue
      endif

      if (posxx.lt.-dx1) then
         posx(tte)=posx(tte)-1
         Xnub(tte)=Xnub(tte)+dx1

         do 1510 k=0,nz1+1
            do 1511 j=0,nx1+1
               do 1512 i=nx1,0,-1
                  U1(i+1,j,k)=U1(i,j,k)
                  V1(i+1,j,k)=V1(i,j,k)
                  W1(i+1,j,k)=W1(i,j,k)
                  Pres1(i+1,j,k)=Pres1(i,j,k)
                  U2(i+1,j,k)=U2(i,j,k)
                  V2(i+1,j,k)=V2(i,j,k)
                  W2(i+1,j,k)=W2(i,j,k)
                  Pres2(i+1,j,k)=Pres2(i,j,k)
                  Titaa1(i+1,j,k)=Titaa1(i,j,k)
                  Qvap1(i+1,j,k)=Qvap1(i,j,k)
                  Qgot1(i+1,j,k)=Qgot1(i,j,k)
                  Qllu1(i+1,j,k)=Qllu1(i,j,k)
                  Qcri1(i+1,j,k)=Qcri1(i,j,k)
                  Aer1(i+1,j,k)=Aer1(i,j,k)
                  Fcalo(i+1,j,k)=Fcalo(i,j,k)
1512           continue
               i=0
               U1(i,j,k)=U1(i+1,j,k)
               V1(i,j,k)=V1(i+1,j,k)
               W1(i,j,k)=W1(i+1,j,k)
               Pres1(i,j,k)=Pres1(i+1,j,k)
               U2(i,j,k)=U2(i+1,j,k)
               V2(i,j,k)=V2(i+1,j,k)
               W2(i,j,k)=W2(i+1,j,k)
               Pres2(i,j,k)=Pres2(i+1,j,k)
               Titaa1(i,j,k)=Titaa1(i+1,j,k)
               Qvap1(i,j,k)=Qvap1(i+1,j,k)
               Qgot1(i,j,k)=Qgot1(i+1,j,k)
               Qllu1(i,j,k)=Qllu1(i+1,j,k)
               Qcri1(i,j,k)=Qcri1(i+1,j,k)
               Aer1(i,j,k)=Aer1(i+1,j,k)
               Fcalo(i,j,k)=0.
1511        continue
            i=1
            j=0
            U1(i,j,k)=(U1(i+1,j,k)+U1(i,j+1,k))/2.
            V1(i,j,k)=(V1(i+1,j,k)+V1(i,j+1,k))/2.
            W1(i,j,k)=(W1(i+1,j,k)+W1(i,j+1,k))/2.
            Pres1(i,j,k)=(Pres1(i+1,j,k)+Pres1(i,j+1,k))/2.
            U2(i,j,k)=(U2(i+1,j,k)+U2(i,j+1,k))/2.
            V2(i,j,k)=(V2(i+1,j,k)+V2(i,j+1,k))/2.
            W2(i,j,k)=(W2(i+1,j,k)+W2(i,j+1,k))/2.
            Pres2(i,j,k)=(Pres2(i+1,j,k)+Pres2(i,j+1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i+1,j,k)+Titaa1(i,j+1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i+1,j,k)+Qvap1(i,j+1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i+1,j,k)+Qgot1(i,j+1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i+1,j,k)+Qllu1(i,j+1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i+1,j,k)+Qcri1(i,j+1,k))/2.
            Aer1(i,j,k)=(Aer1(i+1,j,k)+Aer1(i,j+1,k))/2.
            Fcalo(i,j,k)=0.
            j=nx1+1
            U1(i,j,k)=(U1(i+1,j,k)+U1(i,j-1,k))/2.
            V1(i,j,k)=(V1(i+1,j,k)+V1(i,j-1,k))/2.
            W1(i,j,k)=(W1(i+1,j,k)+W1(i,j-1,k))/2.
            Pres1(i,j,k)=(Pres1(i+1,j,k)+Pres1(i,j-1,k))/2.
            U2(i,j,k)=(U2(i+1,j,k)+U2(i,j-1,k))/2.
            V2(i,j,k)=(V2(i+1,j,k)+V2(i,j-1,k))/2.
            W2(i,j,k)=(W2(i+1,j,k)+W2(i,j-1,k))/2.
            Pres2(i,j,k)=(Pres2(i+1,j,k)+Pres2(i,j-1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i+1,j,k)+Titaa1(i,j-1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i+1,j,k)+Qvap1(i,j-1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i+1,j,k)+Qgot1(i,j-1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i+1,j,k)+Qllu1(i,j-1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i+1,j,k)+Qcri1(i,j-1,k))/2.
            Aer1(i,j,k)=(Aer1(i+1,j,k)+Aer1(i,j-1,k))/2.
            Fcalo(i,j,k)=0.
1510     continue
      endif

!*    corrimiento en y

      if (posyy.gt.dx1) then
         posy(tte)=posy(tte)+1
         Ynub(tte)=Ynub(tte)-dx1


         do 1520 k=0,nz1+1
            do 1521 i=0,nx1+1
               do 1522 j=1,nx1+1
                  U1(i,j-1,k)=U1(i,j,k)
                  V1(i,j-1,k)=V1(i,j,k)
                  W1(i,j-1,k)=W1(i,j,k)
                  Pres1(i,j-1,k)=Pres1(i,j,k)
                  U2(i,j-1,k)=U2(i,j,k)
                  V2(i,j-1,k)=V2(i,j,k)
                  W2(i,j-1,k)=W2(i,j,k)
                  Pres2(i,j-1,k)=Pres2(i,j,k)
                  Titaa1(i,j-1,k)=Titaa1(i,j,k)
                  Qvap1(i,j-1,k)=Qvap1(i,j,k)
                  Qgot1(i,j-1,k)=Qgot1(i,j,k)
                  Qllu1(i,j-1,k)=Qllu1(i,j,k)
                  Qcri1(i,j-1,k)=Qcri1(i,j,k)
                  Aer1(i,j-1,k)=Aer1(i,j,k)
                  Fcalo(i,j-1,k)=Fcalo(i,j,k)
1522           continue
               j=nx1+1
               U1(i,j,k)=U1(i,j-1,k)
               V1(i,j,k)=V1(i,j-1,k)
               W1(i,j,k)=W1(i,j-1,k)
               Pres1(i,j,k)=Pres1(i,j-1,k)
               U2(i,j,k)=U2(i,j-1,k)
               V2(i,j,k)=V2(i,j-1,k)
               W2(i,j,k)=W2(i,j-1,k)
               Pres2(i,j,k)=Pres2(i,j-1,k)
               Titaa1(i,j,k)=Titaa1(i,j-1,k)
               Qvap1(i,j,k)=Qvap1(i,j-1,k)
               Qgot1(i,j,k)=Qgot1(i,j-1,k)
               Qllu1(i,j,k)=Qllu1(i,j-1,k)
               Qcri1(i,j,k)=Qcri1(i,j-1,k)
               Aer1(i,j,k)=Aer1(i,j-1,k)
               Fcalo(i,j,k)=0.
1521        continue
            j=nx1
            i=0
            U1(i,j,k)=(U1(i,j-1,k)+U1(i+1,j,k))/2.
            V1(i,j,k)=(V1(i,j-1,k)+V1(i+1,j,k))/2.
            W1(i,j,k)=(W1(i,j-1,k)+W1(i+1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j-1,k)+Pres1(i+1,j,k))/2.
            U2(i,j,k)=(U2(i,j-1,k)+U2(i+1,j,k))/2.
            V2(i,j,k)=(V2(i,j-1,k)+V2(i+1,j,k))/2.
            W2(i,j,k)=(W2(i,j-1,k)+W2(i+1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j-1,k)+Pres2(i+1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j-1,k)+Titaa1(i+1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j-1,k)+Qvap1(i+1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j-1,k)+Qgot1(i+1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j-1,k)+Qllu1(i+1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j-1,k)+Qcri1(i+1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j-1,k)+Aer1(i+1,j,k))/2.
            Fcalo(i,j,k)=0.
            i=nx1+1
            U1(i,j,k)=(U1(i,j-1,k)+U1(i-1,j,k))/2.
            V1(i,j,k)=(V1(i,j-1,k)+V1(i-1,j,k))/2.
            W1(i,j,k)=(W1(i,j-1,k)+W1(i-1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j-1,k)+Pres1(i-1,j,k))/2.
            U2(i,j,k)=(U2(i,j-1,k)+U2(i-1,j,k))/2.
            V2(i,j,k)=(V2(i,j-1,k)+V2(i-1,j,k))/2.
            W2(i,j,k)=(W2(i,j-1,k)+W2(i-1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j-1,k)+Pres2(i-1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j-1,k)+Titaa1(i-1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j-1,k)+Qvap1(i-1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j-1,k)+Qgot1(i-1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j-1,k)+Qllu1(i-1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j-1,k)+Qcri1(i-1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j-1,k)+Aer1(i-1,j,k))/2.
            Fcalo(i,j,k)=0.
1520     continue
      endif

      if (posyy.lt.-dx1) then
         posy(tte)=posy(tte)-1
         Xnub(tte)=Xnub(tte)+dx1

         do 1530 k=0,nz1+1
            do 1531 i=0,nx1+1
               do 1532 j=nx1,0,-1
                  U1(i,j+1,k)=U1(i,j,k)
                  V1(i,j+1,k)=V1(i,j,k)
                  W1(i,j+1,k)=W1(i,j,k)
                  Pres1(i,j+1,k)=Pres1(i,j,k)
                  U2(i,j+1,k)=U2(i,j,k)
                  V2(i,j+1,k)=V2(i,j,k)
                  W2(i,j+1,k)=W2(i,j,k)
                  Pres2(i,j+1,k)=Pres2(i,j,k)
                  Titaa1(i,j+1,k)=Titaa1(i,j,k)
                  Qvap1(i,j+1,k)=Qvap1(i,j,k)
                  Qgot1(i,j+1,k)=Qgot1(i,j,k)
                  Qllu1(i,j+1,k)=Qllu1(i,j,k)
                  Qcri1(i,j+1,k)=Qcri1(i,j,k)
                  Aer1(i,j+1,k)=Aer1(i,j,k)
                  Fcalo(i,j+1,k)=Fcalo(i,j,k)
1532           continue
               j=0
               U1(i,j,k)=U1(i,j-1,k)
               V1(i,j,k)=V1(i,j-1,k)
               W1(i,j,k)=W1(i,j-1,k)
               Pres1(i,j,k)=Pres1(i,j-1,k)
               U2(i,j,k)=U2(i,j-1,k)
               V2(i,j,k)=V2(i,j-1,k)
               W2(i,j,k)=W2(i,j-1,k)
               Pres2(i,j,k)=Pres2(i,j-1,k)
               Titaa1(i,j,k)=Titaa1(i,j-1,k)
               Qvap1(i,j,k)=Qvap1(i,j-1,k)
               Qgot1(i,j,k)=Qgot1(i,j-1,k)
               Qllu1(i,j,k)=Qllu1(i,j-1,k)
               Qcri1(i,j,k)=Qcri1(i,j-1,k)
               Aer1(i,j,k)=Aer1(i,j-1,k)
               Fcalo(i,j,k)=0.
1531        continue
            j=1
            i=0
            U1(i,j,k)=(U1(i,j+1,k)+U1(i+1,j,k))/2.
            V1(i,j,k)=(V1(i,j+1,k)+V1(i+1,j,k))/2.
            W1(i,j,k)=(W1(i,j+1,k)+W1(i+1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j+1,k)+Pres1(i+1,j,k))/2.
            U2(i,j,k)=(U2(i,j+1,k)+U2(i+1,j,k))/2.
            V2(i,j,k)=(V2(i,j+1,k)+V2(i+1,j,k))/2.
            W2(i,j,k)=(W2(i,j+1,k)+W2(i+1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j+1,k)+Pres2(i+1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j+1,k)+Titaa1(i+1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j+1,k)+Qvap1(i+1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j+1,k)+Qgot1(i+1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j+1,k)+Qllu1(i+1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j+1,k)+Qcri1(i+1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j+1,k)+Aer1(i+1,j,k))/2.
            Fcalo(i,j,k)=0.
            i=nx1+1
            U1(i,j,k)=(U1(i,j+1,k)+U1(i-1,j,k))/2.
            V1(i,j,k)=(V1(i,j+1,k)+V1(i-1,j,k))/2.
            W1(i,j,k)=(W1(i,j+1,k)+W1(i-1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j+1,k)+Pres1(i-1,j,k))/2.
            U2(i,j,k)=(U2(i,j+1,k)+U2(i-1,j,k))/2.
            V2(i,j,k)=(V2(i,j+1,k)+V2(i-1,j,k))/2.
            W2(i,j,k)=(W2(i,j+1,k)+W2(i-1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j+1,k)+Pres2(i-1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j+1,k)+Titaa1(i-1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j+1,k)+Qvap1(i-1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j+1,k)+Qgot1(i-1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j+1,k)+Qllu1(i-1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j+1,k)+Qcri1(i-1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j+1,k)+Aer1(i-1,j,k))/2.
            Fcalo(i,j,k)=0.
1530     continue

      endif

      posxx=posx(tte)*dx1+Xnub(tte)
      posyy=posy(tte)*dx1+Ynub(tte)

   end subroutine corrinu2_init
end module corrinu2
