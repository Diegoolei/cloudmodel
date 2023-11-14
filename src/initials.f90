module velpre01
   integer :: t
   real :: dvx, dvy, dvz, diver
   real :: dprex, dprey, dprez
   real :: vel0, vel1, vel2, vel3
   real :: presi, presix, presiy, presiz
   real :: facx, facy, facz
   integer :: i, j, k
   real :: prom1, prom, kkk, presprom, nnn
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
   integer tt,t1,t2,n,m,l,i,j,k,lll,s,iT,tte
   integer lvapneg,llluneg,lcrineg,laerneg,lnieneg,lgraneg,yy
   real T,P
   real Dv,Lvl,Lvs,Lsl,Vis,Qvap,Qliq,densi,nu
   real Lsl00
   real Eaucn,Eaccn,Eacng
   real*8 qgotaux,qvapaux,qlluaux,qcriaux,qnieaux,qgraaux,aeraux
   real*8 auxx,auxy,auxz
   real*8 Taux,Qvapneg,aerneg
   real Naer,dqgot,dqcri,daer,daer2
   real Fcal
   real elvs,esvs,e1,rl,rs,dden0z
   real aux,aux1,aux2,aux3,aux4
   real*8 ener,ener1,ener2,ener3,ener4,ener5,qv,qg,daitot
   real*8 vapt1,vapt2,vapt3,vapt4
   real*8 gott1,gott2,gott3,gott4
   real*8 aert1,aert2,aert3,aert4
   real*8 totnuc,totmic
   real cks,turbu,lapla
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

module inicio11
   USE dimen
   implicit none

   real :: equis, ygrie, zeta
   real :: G1, centx, centy, centz, cenaerx, cenaery, cenaerz
   real :: sigmat, sigmaa, radiomed, temper, aerper
   real :: elv1, rel1, tem1

   integer :: i, j, k, n

   real :: a0, a1, a2, a3, a4, a5, a6
   real :: b0, b1, b2, b3, b4, b5, b6
   real :: aux, gam, Tk, Qvaptot, aertot
   real :: vb, vc, vh, zeta1, zeta2

contains

   subroutine inicio11_init()
      ! Initialize the variables with their respective values
      centx = (nx1 + 1.) * dx1 / 2.
      centy = (nx1 + 1.) * dx1 / 2.
      centz = 0.
      cenaerx = (nx1 + 1.) * dx1 / 2. + 4000.
      cenaery = (nx1 + 1.) * dx1 / 2. + 1000.
      cenaerz = 0.
      sigmat = 2 * 1000.**2.
      sigmaa = 200.**2.
      radiomed = 2000.
      temper = .7
      aerper = 10000.

      a0 = 6.10780
      a1 = 4.43652e-1
      a2 = 1.42895e-2
      a3 = 2.65065e-4
      a4 = 3.03124e-6
      a5 = 2.03408e-8
      a6 = 6.13682e-11

      b0 = 6.10918
      b1 = 5.03470e-1
      b2 = 1.88601e-2
      b3 = 4.17622e-4
      b4 = 5.82472e-6
      b5 = 4.83880e-8
      b6 = 1.83883e-10
   end subroutine inicio11_init

end module inicio11

module mode20
   implicit none

   integer :: tt, t1, t2, n, m, l, i, j, k, lll, s, iT, tte
   integer :: lvapneg, llluneg, lcrineg, laerneg, lnieneg, lgraneg, yy

   real :: T, P
   real :: Dv, Lvl, Lvs, Lsl, Vis, Qvap, Qliq, densi, nu
   real :: Lsl00
   real :: Eaucn, Eaccn, Eacng

   real(8) :: qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, aeraux
   real(8) :: auxx, auxy, auxz
   real(8) :: Taux, Qvapneg, aerneg
   real :: Naer, dqgot, dqcri, daer, daer2
   real :: Fcal
   real :: elvs, esvs, e1, rl, rs, dden0z
   real :: aux, aux1, aux2, aux3, aux4
   real(8) :: ener, ener1, ener2, ener3, ener4, ener5, qv, qg, daitot

   real(8) :: vapt1, vapt2, vapt3, vapt4
   real(8) :: gott1, gott2, gott3, gott4
   real(8) :: aert1, aert2, aert3, aert4
   real(8) :: totnuc, totmic

   real(8) :: impx, impy, Qagua, Qaguat
   real(8) :: Xnub(5000), Ynub(5000)
   real :: posxx, posyy, zmed
   integer :: posx(-3:5000), posy(-3:5000), spos

   character(len=50) :: bre
   character(len=12) :: nombre
   character(len=2) :: tie

   integer :: laux1, laux2, maux1, maux2, naux2
   integer :: umax, umin, vmax, vmin, wmax, wmin, titamax, titamin
   integer :: qvapmax, qvapmin, qgotmax, qllumax, qcrimax, qniemax, qgramax
   integer :: aermax
   integer :: lumax, mumax, numax, lumin, mumin, numin
   integer :: lvmax, mvmax, nvmax, lvmin, mvmin, nvmin
   integer :: lwmax, mwmax, nwmax, lwmin, mwmin, nwmin
   integer :: ltitamax, mtitamax, ntitamax, ltitamin, mtitamin
   integer :: ntitamin, lqvapmax, mqvapmax, nqvapmax, lqvapmin
   integer :: mqvapmin, nqvapmin, lqgotmax, mqgotmax, nqgotmax
   integer :: lqllumax, mqllumax, nqllumax, laermax, maermax, naermax
   integer :: lqcrimax, mqcrimax, nqcrimax
   integer :: lqniemax, mqniemax, nqniemax
   integer :: lqgramax, mqgramax, nqgramax
   integer :: qgottot, qllutot, qcritot, qnietot, qgratot

   real :: cks, turbu, lapla

contains

   subroutine mode20_init()
      USE cant01
      USE dimen
      USE const
      bre = '01020304050607080910111213141516171819202122232425'
      tie = '31'

      ini = 0                                !inicio por vez primera= 0
      t1 = 0                                 !paso a inicio (si ini=0->t1=0)
      ltt = 250.
      ltg = 10.                              !tiempo de grabacion
      lte = 30.                              !tiempo de grabacion estadistica
      ltb = 10.                              !tiempo de grabacion de backup

      ctur = 0.5

      pro1 = 1. - 2e-2 * (dt1 / 5.)
      pro2 = (1. - pro1) / 6.
      pro3 = 1. - 2e-2 * (dt1 / 5.)
      pro4 = (1. - pro1) / 4.

      lt1 = nint(ltt / dt1)
      lt2 = nint(dt1 / dt2)
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
!$$
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


!            write(*,*) qcrimax,Qcri1(i,j,k),i,j,k
!            pause

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


!        write(*,*) qcrimax,Qcri1(lqcrimax+i,mqcrimax+j,nqcrimax+k),
!     &             lqcrimax+i,mqcrimax+j,nqcrimax+k


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

!$$
      write(30,710) umax,umin,vmax,vmin,wmax,wmin,titamax,titamin&
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


      write(33,715) qgottot,qllutot,qcritot,qnietot,qgratot

710   format(16i5,48i4)
715   format(5i9)
   end subroutine estad03_init
end module estad03
