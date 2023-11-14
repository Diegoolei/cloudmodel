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
!> Condiciones iniciales
   real :: centx,centy,centz,cenaerx,cenaery,cenaerz
   real :: sigmat,sigmaa,radiomed,temper,aerper
   real :: a0,a1,a2,a3,a4,a5,a6
   real :: b0,b1,b2,b3,b4,b5,b6
   real :: G1,equis, ygrie, zeta
   real :: elv1,rel1,tem1
   integer :: i,j,k,n
   real :: aux,gam,Tk,Qvaptot,aertot
   real :: vb,vc,vh,zeta1,zeta2

   public :: inicio11_init

contains
   subroutine inicio11_init()
      USE dimen
      implicit none
      centx=(nx1+1.)*dx1/2.           !Coord x de la perturbacion inicial
      centy=(nx1+1.)*dx1/2.           !Coord y de la perturbacion inicial
      centz=0.                        !Coord z de la perturbacion inicial
      cenaerx=(nx1+1.)*dx1/2.+4000.   !Coord x de la perturbacion de aerosoles
      cenaery=(nx1+1.)*dx1/2.+1000.   !Coord y de la perturbacion de aerosoles
      cenaerz=0.                      !Coord z de la perturbacion de aerosoles
      sigmat=2*1000.**2.              !Decaimiento en z de la perturbacion en T
      sigmaa=200.**2.                 !Decaimiento en z de la perturbacion en A
      radiomed=2000.                  !Ancho de la perturbacion
      temper=.7                       !Perturbacion maxima de temperatura
      aerper=10000.                   !Perturbacion maxima de aerosoles

      a0=6.10780
      a1=4.43652e-1
      a2=1.42895e-2
      a3=2.65065e-4
      a4=3.03124e-6
      a5=2.03408e-8
      a6=6.13682e-11

      b0=6.10918
      b1=5.03470e-1
      b2=1.88601e-2
      b3=4.17622e-4
      b4=5.82472e-6
      b5=4.83880e-8
      b6=1.83883e-10
   end subroutine inicio11_init
end module inicio11
