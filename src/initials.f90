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

module model_initial_conditions
contains
   !     Condiciones iniciales para las variables dinamicas
   !     Corresponde a una nube con hielo
   !     Calcula la presion inicial en forma iterativa, primero supone aire
   !      seco y luego integra considerando la densidad total, incluyendo el
   !      vapor.
   !     Solo hay una perturbacion en temperatura para iniciar la conveccion
   !     Incluye viento de corte, tipo frente

   subroutine initial_conditions()
      USE cant01
      USE dimen
      USE permic
      USE perdim
      USE const
      USE estbas
      USE config

      implicit none

      real equis,ygrie,zeta,G1,centx,centy,centz,sigmat&
         ,sigmaa,radiomed,temper,aerper,elv1,rel1,tem1,a0,a1,a2,a3,a4,a5,a6&
         ,b0,b1,b2,b3,b4,b5,b6,aux,gam,Tk,Qvaptot,aertot,vb,vc,vh,zeta1
      integer i,j,k,n,unit

      centx=(nx1+1.)*dx1/2.           !Coord x de la perturbacion inicial
      centy=(nx1+1.)*dx1/2.           !Coord y de la perturbacion inicial
      centz=0.                        !Coord z de la perturbacion inicial
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


      call PP(G,Rd,dx1,nz,Presi0,P00)

      !**   viento de corte
      vc=1.5/3000**2.
      vh=5000./pi

      do concurrent(k=0:nz1)
         zeta=k*dx1
         if (zeta  <=  500.) then
            UU(k)=0.
            VV(k)=0.

         elseif (zeta  <=  2000.) then
            zeta1=zeta-500.
            aux=4.*(zeta1/1500.)**2.
            UU(k)=aux
            VV(k)=0.

         elseif (zeta  <=  9000.) then
            zeta1=zeta-2000.
            vb=zeta1/7000.
            UU(k)=4.-10.*vb**2.
            VV(k)=3.*vb**.5
            !
         else
            zeta1=zeta-9000.
            UU(k)=-6.+4.*(zeta1/9000.)**2.
            VV(k)=3.-5.*(zeta1/9000.)**.5
            !
         endif

         UU(k)=UU(k)*0.7
         VV(k)=VV(k)*0.7
      end do


      !**   calculo de 'constantes' que dependen de T

      open(newunit=unit,file=output_directory//"ccc", access='append')
      do concurrent (k=210:313)
         Tk=k-T0
         Tvis(k)=4.9e-8*Tk+Vis0
         if (k < 273.15) Tvis(k)=Tvis(k)-1.2e-10*Tk**2.

         !calores latentes de evaporacion, fusion y sublimacion
         gam=.167+3.67e-4*k
         Tlvl(k)=Lvl0*(T0/k)**gam
         Tlsl(k)=(Lsl0+0.485*Tk-2.5e-3*Tk**2.)*4180.
         Tlvs(k)=Tlvl(k)+Tlsl(k)

         !tension de vapor de saturacion liquido y solido
         aux=a3+Tk*(a4+Tk*(a5+Tk*a6))
         aux=a0+Tk*(a1+Tk*(a2+Tk*aux))
         Telvs(k)=aux*100.
         aux=b3+Tk*(b4+Tk*(b5+Tk*b6))
         aux=b0+Tk*(b1+Tk*(b2+Tk*aux))
         Tesvs(k)=aux*100.
         if (k < 220) then
            aux=Tlvl(220)/Rv*(1./220.-1./k)
            Telvs(k)=Telvs(220)*exp(aux)
            aux=Tlvs(220)/Rv*(1./220.-1./k)
            Tesvs(k)=Tesvs(220)*exp(aux)
         endif

         !cambio por las expresiones de Straka
         Eautcn(k)=10.**(.035*(Tk)-.7)
         Eacrcn(k)=exp(.09*Tk)
         write(unit,*) k,Tvis(k),Tlvl(k),Tlsl(k),Tlvs(k),Telvs(k),Tesvs(k),&
            Eautcn(k),Eacrcn(k)
      end do
      close(unit)

      !**   condiciones de tiempo bueno
      ! TT_f not pure function, do concurrent not allowed
      do k=-1,nz1+2
         do concurrent(i=-1:nx1+2, j=-1:nx1+1)
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
            Qcri1(i,j,k)=0.
            Qcri2(i,j,k)=0.
            Qnie1(i,j,k)=0.
            Qnie2(i,j,k)=0.
            Qgra1(i,j,k)=0.
            Qgra2(i,j,k)=0.
            aer1(i,j,k)=0.
            aer2(i,j,k)=0.
         end do

         ! cantidades base
         zeta=k*dx1
         Temp0(k)=TT_f(zeta)
         Den0(k)=Presi0(k)/Rd/Temp0(k)
         Tita0(k)=Temp0(k)*(P00/Presi0(k))**Kapa
         Pres00(k)=Temp0(k)/Tita0(k)
         aer0(k)=10000.*exp(-zeta/2500.)
      end do

      Temp0(-1)=Temp0(0)
      Den0(-1)=Den0(0)
      Tita0(-1)=Tita0(0)
      Pres00(-1)=Pres00(0)
      aer0(-1)=-aer0(0)

      do concurrent(k=0:nz1)
         do concurrent(i=1:nx1, j=1:nx1)
            !perturbaciones iniciales en la temperatura y en los aerosoles
            zeta=k*dx1
            equis=i*dx1
            ygrie=j*dx1
            G1=exp(-((centx-equis)**2.+(centy-ygrie)**2.)*.5&
               /radiomed**2.)
            Titaa1(i,j,k)=temper*exp(-(zeta-centz)**2./sigmat)*G1
            if (Titaa1(i,j,k) < 1e-5) Titaa1(i,j,k)=0.
            aer1(i,j,k)=aerper*exp(-zeta**2./sigmaa)*G1
         end do

         !vapor base
         Tem1=Temp0(k)
         if (zeta <= 500) then
            rel1=.55+.05*zeta/500.
         else if (zeta <=  1500.) then
            rel1=.6
         else if (zeta <= 4000) then
            rel1=.6-(zeta-1500)/2500.*.25
         else if (zeta <= 7000) then
            rel1=.35-(zeta-4000.)/3000.*.25
         else if(zeta > 7000) then
            rel1=.1-(zeta-7000)/3000.*.02
         endif
         n=int(Tem1)
         aux=Tem1-n
         elv1=Telvs(n)*(1-aux)+Telvs(n+1)*aux
         Qvap0(k)=rel1*elv1/Rv/Tem1

         !recalculo de la densidad
         Den0(k)=Den0(k)+Qvap0(k)
      end do

      !**   Velocidad terminal para gota de lluvia, cte que depende de P
      do concurrent(k=1:nz1+1)
         Av(2*k-1)=Av0*((P00/Presi0(k-1))**.286+(P00/Presi0(k))**.286)/2. !puntos intermedios
         Av(2*k)=Av0*(P00/Presi0(k))**.286
      end do

      !**   Velocidad terminal para la nieve, cte que depende de P
      do concurrent(k=1:nz1+1)
         Vtnie(2*k-1)=Vtnie0*((P00/Presi0(k-1))**.3+(P00/Presi0(k))**.3)/2. !puntos intermedios&
         Vtnie(2*k)=Vtnie0*(P00/Presi0(k))**.3
      end do

      !**   Velocidad terminal para el granizo, cte que depende de z
      do concurrent(k=0:nz1+1)
         aux=2.754*rhogra**.605
         Vtgra0(2*k)=aux/Tvis(Temp0(k))**.21/Den0(k)**.395
      end do

      do concurrent(k=1:nz1+1)
         Vtgra0(2*k-1)=(Vtgra0(2*k-2)+Vtgra0(2*k))/2.  ! punto intermedio
      end do

      !**************************************************************
      !    Recalculo de la Presion y de Tita

      !    Recalculo de la Presion a partir de la densidad

      call PP2(G,dx1,Den0,Presi0,P00)

      open(newunit=unit,file=output_directory//"inic03.sa", access='append')
      do concurrent(k=0:nz1)
         Tita0(k)=Temp0(k)*(P00/Presi0(k))**Kapa
         Pres00(k)=Temp0(k)/Tita0(k)
         cc2(k)=Cp*Rd*Tita0(k)*Pres00(k)/Cv

         write(unit,210) k,Temp0(k),Tita0(k),Presi0(k),Pres00(k),&
            Den0(k),aer0(k),Qvap0(k),UU(k),VV(k)
      end do
      close(unit)


      Tita0(-1)=Tita0(0)
      Pres00(-1)=Pres00(0)
      Den0(-1)=Den0(0)
      Qvap0(-1)=0

      do concurrent(i=1:nx1, j=1:nx1)
         Pres1(i,j,0)=Pres1(i,j,1)
         Pres1(i,j,-1)=Pres1(i,j,1)
         Pres2(i,j,0)=Pres1(i,j,1)
         Pres2(i,j,-1)=Pres1(i,j,1)
         Titaa1(i,j,0)=Titaa1(i,j,1)
         Titaa1(i,j,-1)=Titaa1(i,j,1)
         Qvap1(i,j,0)=Qvap1(i,j,1)
         Qvap1(i,j,-1)=Qvap1(i,j,1)
      end do

      !     calculo del Qvaprel
      Qvaptot=0.
      do concurrent(k=1:nz1)
         Qvaptot=Qvaptot+Qvap0(k)
      end do
      do concurrent(k=1:nz1)
         Qvaprel(k)=Qvap0(k)/Qvaptot
      end do

      !     calculo del aerrel
      aertot=0.
      do concurrent(k=1:nz1)
         aertot=aertot+aer0(k)
      end do
      do concurrent(k=1:nz1)
         aerrel(k)=aer0(k)/aertot
      end do

210   format(I3,9E12.4)
      return

   end subroutine initial_conditions

   !*********************************************************

   function TT_f (zeta)
      real :: a, xx, TT_f
      a = 298.15
      if (zeta <= 2000) then
         TT_f = a - 9.e-3 * zeta
      elseif (zeta <= 5500) then
         xx = zeta - 2000.
         TT_f = a - 18. - xx * (9.e-3 - 2e-3 * xx / 3500. / 2.)
      elseif (zeta <= 9000) then
         xx = zeta - 5500.
         TT_f = a - 46. - 7e-3 * xx
      elseif (zeta <= 11000) then
         xx = zeta - 9000
         TT_f = a - 70.5 - 7e-3 * xx + 1.75e-6 * xx**2.
      elseif (zeta <= 12000) then
         TT_f = a - 77.5
      else
         xx = zeta - 12000
         TT_f = a - 77.5 + 50. * (xx / 9000.)**2.
      endif
   end function TT_f
   !*****************************************************
   subroutine PP(G,Rd,dx,nz1,Pres,Pres0)
      integer k,nz1,nx4
      parameter (nx4=500)
      real Pres(-3:nz1+3)
      real integ(-2:nx4+2)
      real Pres0
      real G,Rd,dx,dx4
      real zetaa,zetam,zetad
      dx4=dx/4.

      integ(0)=0
      ! TT_f not pure function, do concurrent not allowed
      do k=1,nx4
         zetaa=(2*k-2)*dx4
         zetam=(2*k-1)*dx4
         zetad=(2*k)*dx4
         ya=1/TT_f(zetaa)
         ym=1/TT_f(zetam)
         yd=1/TT_f(zetad)
         integ(k)=integ(k-1)+ya+4*ym+yd
      end do

      do concurrent(k=1:nz1+2)
         Pres(k)=Pres0*exp(-G/Rd*(integ(2*k)*dx4/3))
      end do
      Pres(0)=Pres0
      Pres(-1)=Pres0

      return
   end subroutine PP

   !*****************************************************
   subroutine PP2(G,dx,Den0,Pres00,Pres0)
      USE dimen
      integer k
      real Pres00(-3:nz1+3)
      real Den0(-3:nz1+3)
      real Den00(-3:3*nz1+3)
      real integ(-3:3*nz1+3)
      real Pres0
      real G,dx

      do concurrent(k=0:nz1-1)
         Den00(2*k)=Den0(k)
         Den00(2*k+1)=(Den0(k)+Den0(k+1))/2.
      end do
      Den00(2*nz1)=Den0(nz1)
      Den00(2*nz1+1)=2.*Den0(nz1)-Den00(2*nz1-1)

      integ(0)=0
      do concurrent(k=1:nz1)
         ya=Den00(2*k-1)
         ym=Den00(2*k)
         yd=Den00(2*k+1)
         integ(k)=integ(k-1)+ya+4*ym+yd
      end do
      do concurrent(k=1:nz1)
         Pres00(k)=Pres0-G*integ(k)*dx/6.
      end do
      Pres00(0)=Pres0
      Pres00(-1)=Pres0
      return
   end subroutine PP2
end module model_initial_conditions

module model_var
   real :: T, P, Dv, Lvl, Lvs, Lsl, Vis, Qvap, Qliq, densi, nu, Lsl00, Eaucn,&
      Eaccn, Eacng, Naer, dqgot, dqcri, daer, daer2, Fcal, elvs, esvs, e1, rl,&
      rs, dden0z, aux, aux1, aux2, aux3, aux4, posxx, posyy,cks, turbu, lapla

   real(8) :: qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, aeraux,&
      auxx, auxy, auxz, Taux, Qvapneg, aerneg, ener, ener1, ener2, ener3,&
      ener4, ener5, qv, qg, daitot, vapt1, vapt2, vapt3, vapt4, gott1, gott2,&
      gott3, gott4, aert1, aert2, aert3, aert4, totnuc, totmic

   real(8) :: Xnub(5000), Ynub(5000)
   integer :: posx(-3:5000), posy(-3:5000)

   character(len=3) :: file_number

   integer :: tt, t1, t2, n, m, l, i, j, k, lll, s, iT, tte, lvapneg,&
      llluneg, lcrineg, laerneg, lnieneg, lgraneg, yy

end module model_var

module model_initialization
   implicit none

   real :: zmed
   real(8) :: impx, impy, Qagua, Qaguat
   integer :: spos, laux1, laux2, maux1, maux2, naux2, umax, umin,&
      titamin, qvapmax, qvapmin, qgotmax, qllumax, qcrimax, qniemax,&
      qgramax, aermax, lumax, mumax, numax, lumin, mumin, numin, lvmax,&
      mvmax, nvmax, lvmin, mvmin, nvmin, lwmax, mwmax, nwmax, lwmin,&
      mwmin, nwmin, ltitamax, mtitamax, ntitamax, ltitamin, mtitamin,&
      ntitamin, lqvapmax, mqvapmax, nqvapmax, lqvapmin, mqvapmin, nqvapmin,&
      lqgotmax, mqgotmax, nqgotmax, lqllumax, mqllumax, nqllumax, laermax,&
      maermax, naermax, lqcrimax, mqcrimax, nqcrimax, lqniemax, mqniemax,&
      nqniemax, lqgramax, mqgramax, nqgramax, qgottot, qllutot, qcritot,&
      vmax, vmin, wmax, wmin, titamax, qnietot, qgratot

contains

   subroutine initialize_model()
      USE model_var
      USE cant01
      USE dimen
      USE const
      USE config
      USE estbas, only: Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV,&
         Qvaprel, aerrel
      USE perdim, only: W2, U2, V2, Fcalo, Titaa2, Titaa1, Pres2, Pres1, U1, V1,&
         W1
      USE permic, only: aer1, Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, Qvap2, aer2,&
         Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, Av, Vtgra0, Vtnie
      USE model_initial_conditions, only: initial_conditions
      implicit none
      integer :: unit_number

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
      lt2 = nint(dt1 / dt2)                  ! Proporcion Fisica/Microfisica
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

      if (ini == 0) then
         !Si ini=0 el calculo empieza por primera###
         call initial_conditions()
      else
         !si ini=1 el calculo recomienza desde algun paso
         open(newunit=unit_number,file=output_directory//"inis.da",status='unknown',form='unformatted')
         read(unit_number,*) Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV
         close(unit_number)

         open(newunit=unit_number,file=output_directory//"velos.da",status='unknown',form='unformatted')
         rewind unit_number
         read(unit_number) U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,&
            Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
            Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,&
            aer1,aer2,Fcalo
         close(unit_number)

         open(newunit=unit_number,file=output_directory//"varconz.da",status='unknown',form='unformatted')
         rewind unit_number
         read(unit_number)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn
         close(unit_number)

      endif

      !definicion de calores y presiones de vapor a 0 K
      Lsl00=Lsl0*4180.


   end subroutine initialize_model

   subroutine statistics_init()
      USE model_var
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

      do concurrent(k=1:nz1, i=1:nx1, j=1:nx1)
         if (umax < U1(i,j,k)*100) then
            umax=U1(i,j,k)*100
            lumax=i
            mumax=j
            numax=k
         endif

         if (umin > U1(i,j,k)*100) then
            umin=U1(i,j,k)*100
            lumin=i
            mumin=j
            numin=k
         endif

         if (vmax < V1(i,j,k)*100) then
            vmax=V1(i,j,k)*100
            lvmax=i
            mvmax=j
            nvmax=k
         endif

         if (vmin > V1(i,j,k)*100) then
            vmin=V1(i,j,k)*100
            lvmin=i
            mvmin=j
            nvmin=k
         endif

         if (wmax < W1(i,j,k)*100) then
            wmax=W1(i,j,k)*100
            lwmax=i
            mwmax=j
            nwmax=k
         endif

         if (wmin > W1(i,j,k)*100) then
            wmin=W1(i,j,k)*100
            lwmin=i
            mwmin=j
            nwmin=k
         endif

         if (titamax < Titaa1(i,j,k)*1000) then
            titamax=Titaa1(i,j,k)*1000
            ltitamax=i
            mtitamax=j
            ntitamax=k
         endif

         if (titamin > Titaa1(i,j,k)*1000) then
            titamin=Titaa1(i,j,k)*1000
            ltitamin=i
            mtitamin=j
            ntitamin=k
         endif

         if (qvapmax < Qvap1(i,j,k)*1e6) then
            qvapmax=Qvap1(i,j,k)*1e6
            lqvapmax=i
            mqvapmax=j
            nqvapmax=k
         endif

         if (qvapmin > Qvap1(i,j,k)*1e6) then
            qvapmin=Qvap1(i,j,k)*1e6
            lqvapmin=i
            mqvapmin=j
            nqvapmin=k
         endif

         if (qgotmax < Qgot1(i,j,k)*1e6) then
            qgotmax=Qgot1(i,j,k)*1e6
            lqgotmax=i
            mqgotmax=j
            nqgotmax=k
         endif
         qgottot=qgottot+Qgot1(i,j,k)*1e6

         if (qllumax < Qllu1(i,j,k)*1e6) then
            qllumax=Qllu1(i,j,k)*1e6
            lqllumax=i
            mqllumax=j
            nqllumax=k
         endif
         qllutot=qllutot+Qllu1(i,j,k)*1e6

         if (qcrimax < Qcri1(i,j,k)*1e6) then
            qcrimax=Qcri1(i,j,k)*1e6
            lqcrimax=i
            mqcrimax=j
            nqcrimax=k
         endif
         qcritot=qcritot+Qcri1(i,j,k)*1e6

         if (qniemax < Qnie1(i,j,k)*1e6) then
            qniemax=Qnie1(i,j,k)*1e6
            lqniemax=i
            mqniemax=j
            nqniemax=k
         endif
         qnietot=qnietot+Qnie1(i,j,k)*1e6

         if (qgramax < Qgra1(i,j,k)*1e6) then
            qgramax=Qgra1(i,j,k)*1e6
            lqgramax=i
            mqgramax=j
            nqgramax=k
         endif
         qgratot=qgratot+Qgra1(i,j,k)*1e6

         if (aermax < aer1(i,j,k)/1000) then
            aermax=aer1(i,j,k)/1000
            laermax=i
            maermax=j
            naermax=k
         endif
      end do

      qgotmax=0.
      qcrimax=0.
      qniemax=0.

      do concurrent(i=-1:1, j=-1:1, k=-1:1)
         qgotmax=qgotmax+1e5*Qgot1(lqgotmax+i,mqgotmax+j,nqgotmax+k)
         qcrimax=qcrimax+1e5*Qcri1(lqcrimax+i,mqcrimax+j,nqcrimax+k)
         qniemax=qniemax+1e5*Qnie1(lqniemax+i,mqniemax+j,nqniemax+k)
      end do

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
      write(unit_number,710) umax,umin,vmax,vmin,wmax,wmin,titamax,titamin,&
         qvapmax,qvapmin,qgotmax,qllumax,qcrimax,qniemax,qgramax,aermax,&
         lumax,mumax,numax,lumin,mumin,numin,lvmax,mvmax,nvmax,lvmin,mvmin,&
         nvmin,lwmax,mwmax,nwmax,lwmin,mwmin,nwmin,ltitamax,mtitamax,ntitamax,&
         ltitamin,mtitamin,ntitamin,lqvapmax,mqvapmax,nqvapmax,lqvapmin,mqvapmin,&
         nqvapmin,lqgotmax,mqgotmax,nqgotmax,lqllumax,mqllumax,nqllumax,lqcrimax,&
         mqcrimax,nqcrimax,lqniemax,mqniemax,nqniemax,lqgramax,mqgramax,nqgramax,&
         laermax,maermax,naermax
      close(unit_number)

      open(newunit=unit_number,file=output_directory//"est", ACCESS="append")
      write(unit_number,715) qgottot,qllutot,qcritot,qnietot,qgratot
      close(unit_number)
710   format(16i5,48i4)
715   format(5i9)
   end subroutine statistics_init


   subroutine cloud_position_init()
      USE model_var
      USE lmngot
      USE lmncri
      USE permic
      USE estbas
      USE cant01
!     desplazamientos horizontales a partir de la velocidad media de la nube
!     calculo la altura media de la nube, la velocidad media de la nube
!     es tomada como la velocidad del aire sin perturbar a esa altura
!                  (1/3/2000)

      if (ngot(2) >= 1 .or. ncri(2) > 1) then

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

         do concurrent(k=1:naux2, i=laux1:laux2, j=maux1:maux2)
            Qagua=Qgot1(i,j,k)+Qcri1(i,j,k)+Qllu1(i,j,k)+&
               Qnie1(i,j,k)+Qgra1(i,j,k)
            zmed=zmed+k*Qagua
            Qaguat=Qaguat+Qagua
         end do

         if (Qaguat > 1e-3) then
            zmed=zmed/Qaguat
            Xnub(tte)=Xnub(tte)+UU(nint(zmed))*lte
            Ynub(tte)=Ynub(tte)+VV(nint(zmed))*lte
         endif

      endif

   end subroutine cloud_position_init

   subroutine cloud_movement_init()
      !desplazamiento de la nube
      USE model_var
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

      if (spos  ==  1) then
         posxx=Xnub(tte)
         posyy=Ynub(tte)
      else
         posxx=0.
         posyy=0.
      endif

      posx(tte)=posx(tte-1)
      posy(tte)=posy(tte-1)

!*    corrimiento en x

      if (posxx > dx1) then
         posx(tte)=posx(tte)+1
         Xnub(tte)=Xnub(tte)-dx1

         do concurrent(k=0:nz1+1)
            do concurrent(j=0:nx1+1)
               do concurrent(i=1:nx1+1)
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
               end do
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
            end do
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
         end do
      endif

      if (posxx < -dx1) then
         posx(tte)=posx(tte)-1
         Xnub(tte)=Xnub(tte)+dx1

         do concurrent(k=0:nz1+1)
            do concurrent(j=0:nx1+1)
               do concurrent(i=nx1:0) ! TODO Test this loop
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
               end do
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
            end do
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
         end do
      endif

!*    corrimiento en y

      if (posyy > dx1) then
         posy(tte)=posy(tte)+1
         Ynub(tte)=Ynub(tte)-dx1


         do concurrent(k=0:nz1+1)
            do concurrent(i=0:nx1+1)
               do concurrent(j=1:nx1+1)
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
               end do
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
            end do
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
         end do
      endif

      if (posyy < -dx1) then
         posy(tte)=posy(tte)-1
         Xnub(tte)=Xnub(tte)+dx1

         do concurrent(k=0:nz1+1)
            do concurrent(i=0:nx1+1)
               do concurrent(j=nx1:0) ! TODO Test this loop
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
               end do
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
            end do
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
         end do

      endif

      posxx=posx(tte)*dx1+Xnub(tte)
      posyy=posy(tte)*dx1+Ynub(tte)

   end subroutine cloud_movement_init
end module model_initialization

module aeroana
!> No se usa en MODEL31
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
