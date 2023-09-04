!     include para las condiciones iniciales

      real equis,ygrie,zeta
      real G1,centx,centy,centz,cenaerx,cenaery,cenaerz
      real sigmat,sigmaa,radiomed,temper,aerper

      real TT,elv1,rel1,tem1

      integer i,j,k,n

      real a0,a1,a2,a3,a4,a5,a6
      real b0,b1,b2,b3,b4,b5,b6
      real aux,gam,Tk,Qvaptot,aertot
      real vb,vc,vh,zeta1,zeta2

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


