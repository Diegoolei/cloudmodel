module model_aux
contains
   subroutine vapor_advection()
      !######################## Adveccion de vapores #######################
      USE advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advaer1
      USE permic, only: aer1, Qvap1, Qllu1, Qnie1, Qgra1
      USE perdim, only: W2
      USE dimen, only: nx1
      implicit none
      integer :: i, j
      do concurrent (i=0:nx1+1, j=0:nx1+1)
         advvap1(i,j)=W2(i,j,1)*(Qvap1(i,j,1)+Qvap1(i,j,0))/4.
         advgot1(i,j)=0.
         advllu1(i,j)=W2(i,j,1)*Qllu1(i,j,1)
         if (W2(i,j,1) > 0) advllu1(i,j)=0.
         advaer1(i,j)=W2(i,j,1)*(aer1(i,j,1)+aer1(i,j,0))/4.
         if(W2(i,j,1) < 0) advaer1(i,j)=advaer1(i,j)*1.5
         advcri1(i,j)=0.
         advnie1(i,j)=W2(i,j,1)*Qnie1(i,j,1)
         if (W2(i,j,1) > 0) advnie1(i,j)=0.
         advgra1(i,j)=W2(i,j,1)*Qgra1(i,j,1)
         if (W2(i,j,1) > 0) advgra1(i,j)=0.
      end do
   end subroutine vapor_advection

   subroutine dinamics()
      !########### calculo de la dinamica y de la termodinamica ############
      USE dimen, only: dt1, dx1, nx1, nz1
      USE model_var, only: dden0z, ener1, s, llluneg, lcrineg, lnieneg, lgraneg,&
         lvapneg, laerneg, Qvapneg, aerneg
      USE estbas, only: Den0, Qvap0, aer0
      USE perdim, only: W2, U2, V2, Fcalo
      USE advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advvap2, advgot2, advllu2, advcri2, advnie2, advgra2, advaer1,&
         advaer2
      USE permic, only: Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, Qvap2, aer2
      USE lmngot, only: lgot, mgot, ngot
      USE lmnllu, only: lllu, mllu, nllu
      USE lmncri, only: lcri, mcri, ncri
      USE lmnnie, only: lnie, mnie, nnie
      USE lmngra, only: lgra, mgra, ngra
      implicit none
      integer :: i, j, k, l, m, n
      s=0
      Qvapneg=0.
      lvapneg=0
      aerneg=0.
      laerneg=0
      llluneg=0
      lcrineg=0.
      lnieneg=0.
      lgraneg=0.
      ener1=0.
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
      do k=1,nz1-1
         n=k
         dden0z=(Den0(k+1)-Den0(k-1))/Den0(k)
         call turbu1(n)
         do i=1,nx1
            l=i
            do j=1,nx1
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
               call daeros(l,m,n)
               advaer1(i,j)=advaer2(i,j)

               !limites de la nube
               if(Qgot2(i,j,k) /= 0) then
                  if (i < lgot(1)) lgot(1)=i
                  if (i > lgot(2)) lgot(2)=i
                  if (j < mgot(1)) mgot(1)=j
                  if (j > mgot(2)) mgot(2)=j
                  if (k < ngot(1)) ngot(1)=k
                  if (k > ngot(2)) ngot(2)=k
                  s=s+1
               endif

               !limites de la lluvia
               if(Qllu2(i,j,k) /= 0) then
                  if (i < lllu(1)) lllu(1)=i
                  if (i > lllu(2)) lllu(2)=i
                  if (j < mllu(1)) mllu(1)=j
                  if (j > mllu(2)) mllu(2)=j
                  if (k < nllu(1)) nllu(1)=k
                  if (k > nllu(2)) nllu(2)=k
                  llluneg=1
               endif

               !limites de los cristales
               if(Qcri2(i,j,k) /= 0) then
                  if (i < lcri(1)) lcri(1)=i
                  if (i > lcri(2)) lcri(2)=i
                  if (j < mcri(1)) mcri(1)=j
                  if (j > mcri(2)) mcri(2)=j
                  if (k < ncri(1)) ncri(1)=k
                  if (k > ncri(2)) ncri(2)=k
                  lcrineg=1
               endif

               !limites de la nieve
               if(Qnie2(i,j,k) /= 0) then
                  if (i < lnie(1)) lnie(1)=i
                  if (i > lnie(2)) lnie(2)=i
                  if (j < mnie(1)) mnie(1)=j
                  if (j > mnie(2)) mnie(2)=j
                  if (k < nnie(1)) nnie(1)=k
                  if (k > nnie(2)) nnie(2)=k
                  lnieneg=1
               endif

               !limites del granizo
               if(Qgra2(i,j,k) /= 0) then
                  if (i < lgra(1)) lgra(1)=i
                  if (i > lgra(2)) lgra(2)=i
                  if (j < mgra(1)) mgra(1)=j
                  if (j > mgra(2)) mgra(2)=j
                  if (k < ngra(1)) ngra(1)=k
                  if (k > ngra(2)) ngra(2)=k
                  lgraneg=1
               endif

               if(Qvap0(k)+Qvap2(i,j,k) < 0) then
                  Qvapneg=Qvapneg+Qvap0(k)+Qvap2(i,j,k)
                  lvapneg=1
               endif

               if(aer0(k)+aer2(i,j,k) < 0) then
                  aerneg=aerneg+aer0(k)+aer2(i,j,k)
                  laerneg=1
               endif
            end do
         end do
      end do

   end subroutine dinamics

   subroutine negative_correction
      USE model_var, only: s, llluneg, lcrineg, lnieneg, lgraneg, lvapneg, laerneg,&
         Qvapneg, aerneg
      implicit none
      if(s >= 1) call corgot
      if (llluneg == 1) call corllu
      if (lcrineg == 1) call corcri
      if (lnieneg == 1) call cornie
      if (lgraneg == 1) call corgra
      if (lvapneg == 1) call corvap(Qvapneg)
      if (laerneg == 1) call coraer(aerneg)
   end subroutine negative_correction

   subroutine water_calculation
      !primer calculo de agua (sin laterales)
      USE dimen, only: nx1, nz1
      USE model_var, only: vapt1, gott1, aert1
      USE permic, only: Qvap2, Qgot2, aer2
      USE estbas, only: Qvap0
      implicit none
      integer :: i, j, k
      vapt1=0.
      gott1=0.
      aert1=0.
      do concurrent (i=1:nx1, j=1:nx1, k=1:nz1-1)
         vapt1=vapt1+Qvap2(i,j,k)
         gott1=gott1+Qgot2(i,j,k)
         aert1=aert1+aer2(i,j,k)

         if(Qvap2(i,j,k)+Qvap0(k) < 0) then
            stop
         endif
      end do
   end subroutine water_calculation

   subroutine microphisics_substring
      USE dimen, only: nx1, nz1, dt1, dx1
      USE model_var, only: aux, P, T, Qvap, Naer, densi, Dv, iT, aux2, Vis,&
         esvs, elvs, Lvl, Lsl, Lvs, Eaccn, Eaucn, Eacng, nu, lll, tt, Qliq, e1,&
         rl, rs, yy, daer, dqgot, dqcri, Taux, totnuc, vapt2, gott2, aert2,&
         qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, t2, Lsl00, Fcal,&
         daer2, totmic, vapt3, gott3, aert3, ener2, ener3, ener4, ener5, qv, qg,&
         daitot
      USE perdim, only: Titaa2, Pres2, Tempa1, Fcalo
      USE estbas, only: Qvap0, aer0, Tita0, Temp0, Pres00
      USE cant01, only: ikapa, AA, lt2
      USE const, only: P00, Rd, Dv0, Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs,&
         Eacrcn, Eautcn, T0, Rv, G, Cp
      USE permic, only: Qvap2, aer2, Qgot2, Qllu2, Qcri2, Qnie2, Qgra2
      implicit none
      integer :: k, n, i, l, j, m
      !####################### Sublazo Microfisico #########################
      totnuc=0.
      vapt2=0.
      gott2=0.
      aert2=0.
      totmic=0.
      vapt3=0.
      gott3=0.
      aert3=0.
      ener2=0.
      ener3=0.
      ener4=0.
      ener5=0.
      qv=0.
      qg=0.
      daitot=0.
      do k=1,nz1-1
         n=k
         do i=1,nx1
            l=i
            do j=1,nx1
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
               if (T >= T0) then
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

               if ((rl > 1e-3 .or. rs > 1e-3).and.Naer > 0) then
                  call nuclea(Qvap,Qliq,Naer,T,densi,e1,elvs,esvs,rl,rs,Lvl,Lvs,daer,dqgot,dqcri)
                  Taux=T-Temp0(k)-Tempa1(i,j,k)
                  Titaa2(i,j,k)=T/aux-Tita0(k)
                  if (dqgot > 0) yy=1
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
               if (Qgot2(i,j,k) > 0 .or. dqgot > 0 .or.Qllu2(i,j,k) > 0 .or. Qcri2(i,j,k) > 0 .or.Qnie2(i,j,k) > 0) then
                  qgotaux=Qgot2(i,j,k)
                  if (Qgot2(i,j,k) == 0) qgotaux=0d0
                  qvapaux=Qvap2(i,j,k)+Qvap0(k)
                  qlluaux=Qllu2(i,j,k)
                  if (Qllu2(i,j,k) == 0) qlluaux=0d0
                  qcriaux=Qcri2(i,j,k)
                  if (Qcri2(i,j,k) == 0) then
                     qcriaux=0d0
                  endif
                  qnieaux=Qnie2(i,j,k)
                  if (Qnie2(i,j,k) == 0) qnieaux=0d0
                  qgraaux=Qgra2(i,j,k)
                  if (Qgra2(i,j,k) == 0) qgraaux=0d0
                  Naer=aer2(i,j,k)+aer0(k)
                  T=Tempa1(i,j,k)+Temp0(k)
                  do t2=1,lt2
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
                     call microfis(elvs,esvs,Lvl,Lvs,Lsl,T,Dv,Eaccn,Eaucn,Eacng,&
                        Lsl00,Fcal,n,qvapaux,qgotaux,qlluaux,qcriaux,qnieaux,&
                        qgraaux,Naer,daer2,nu,yy)
                     Fcalo(l,m,n)=Fcalo(l,m,n)+Fcal/dt1/densi
                     Naer=Naer+daer2
                     totmic=totmic+daer2
                  end do

                  Qgot2(i,j,k)=qgotaux
                  Qllu2(i,j,k)=qlluaux
                  Qcri2(i,j,k)=qcriaux
                  Qnie2(i,j,k)=qnieaux
                  Qgra2(i,j,k)=qgraaux
                  Qvap2(i,j,k)=qvapaux-Qvap0(k)
                  aer2(i,j,k)=Naer-aer0(k)
                  Tempa1(i,j,k)=T-Temp0(k)

               endif

               if (Tita0(k) < abs(Titaa2(i,j,k))+200.or.Temp0(k) < abs(Tempa1(i,j,k))+200) then
                  stop
               endif

               if(aer2(i,j,k)+aer0(k) <= 0) then

                  if (aer2(i,j,k)+aer0(k) < -aer0(k)*.05) then
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
            end do
         end do
      end do
   end subroutine microphisics_substring

   subroutine floor_and_ceiling_contour
      USE dimen, only: nx1, nz1, dt1, dx1
      USE perdim, only: Titaa2, Titaa1, W2, U2, V2
      USE cant01, only: dx2
      USE estbas, only: Tita0, Qvap0, aer0
      USE model_var, only: auxx, auxy, auxz, turbu, aeraux, lapla, cks, Qvap,&
         qv, qg
      USE permic, only: Qvap1, Qvap2, Qgot2, Qcri2, Qllu2, Qnie2, Qgra2, aer1,&
         aer2
      implicit none
      integer :: i, j

      Qvap = (qv+qg)/nx1**2.*nz1*.1
      !contornos en el piso y en el techo
      do concurrent (i=1:nx1, j=1:nx1)
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

         if (W2(i,j,0) > 0) then
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
      end do
   end subroutine floor_and_ceiling_contour

   subroutine lateral_contour
      !contornos laterales
      USE dimen, only: nx1, nz1
      USE perdim, only: Titaa2
      USE permic, only: Qvap2, Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, aer2
      implicit none
      integer :: j, k
      do concurrent (k=1:nz1-1, j=1:nx1)
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
      end do
   end subroutine lateral_contour

   subroutine floor_condition_redefinition()
      USE dimen, only: nx1, nz1, dt1, dx1
      USE perdim, only: Titaa2, Titaa1, W2
      USE cant01, only: pro3, pro4, pro1, pro2
      USE permic, only: Qvap1, Qvap2, Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, aer1,&
         aer2, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1
      USE estbas, only: aer0
      USE model_var, only: aeraux
      implicit none
      integer :: i, j, k
      !*** modificada las condiciones en el piso
      !Redefinicion
      do concurrent (i=1:nx1, j=1:nx1)
         k=0
         Titaa1(i,j,k)=pro3*Titaa2(i,j,k)+&
            pro4*((Titaa2(i+1,j,k)+Titaa2(i-1,j,k))+(Titaa2(i,j+1,k)+Titaa2(i,j-1,k)))

         if (abs(Titaa1(i,j,k)) < 1e-10) Titaa1(i,j,k)=0

         Qvap1(i,j,k)=pro3*Qvap2(i,j,k)+&
            pro4*((Qvap2(i+1,j,k)+Qvap2(i-1,j,k))+(Qvap2(i,j+1,k)+Qvap2(i,j-1,k)))


         if (abs(Qvap1(i,j,k)) < 1e-10) Qvap1(i,j,k)=0

         Qgot1(i,j,k)=pro3*Qgot2(i,j,k)+&
            pro4*((Qgot2(i+1,j,k)+Qgot2(i-1,j,k))+(Qgot2(i,j+1,k)+Qgot2(i,j-1,k)))

         if (Qgot1(i,j,k) < 1e-10) Qgot1(i,j,k)=0

         Qllu1(i,j,k)=Qllu2(i,j,k)

         if (Qllu1(i,j,k) < 1e-10) Qllu1(i,j,k)=0

         Qcri1(i,j,k)=pro3*Qcri2(i,j,k)+&
            pro4*((Qcri2(i+1,j,k)+Qcri2(i-1,j,k))+(Qcri2(i,j+1,k)+Qcri2(i,j-1,k)))

         if (Qcri1(i,j,k) < 1e-10) Qcri1(i,j,k)=0

         Qnie1(i,j,k)=pro3*Qnie2(i,j,k)+pro4*((Qnie2(i+1,j,k)&
            +Qnie2(i-1,j,k))+(Qnie2(i,j+1,k)+Qnie2(i,j-1,k)))

         if (Qnie1(i,j,k) < 1e-10) Qnie1(i,j,k)=0

         Qgra1(i,j,k)=Qgra2(i,j,k)

         if (Qgra1(i,j,k) < 1e-10) Qgra1(i,j,k)=0

         aer1(i,j,k)=pro3*aer2(i,j,k)+pro4*((aer2(i+1,j,k)+aer2(i-1,j,k))+(aer2(i,j+1,k)+aer2(i,j-1,k)))

         !correccion cambiando la absorcion de aerosoles
         if ((Qllu1(i,j,1)+Qgra1(i,j,1)) > 1e-6 .and.W2(i,j,1) < 0) then
            aeraux=-W2(i,j,1)*.5*dt1/(dx1/2)
            aer1(i,j,k)=aer1(i,j,k)-(aer1(i,j,k)+aer0(k))*aeraux
         endif

         if (abs(aer1(i,j,k)) < 1e-10) aer1(i,j,k)=0


         do concurrent(k=1:nz1-1)
            Titaa1(i,j,k)=pro1*Titaa2(i,j,k)+&
               pro2*(&
               (Titaa2(i+1,j,k)+Titaa2(i-1,j,k))+&
               (Titaa2(i,j+1,k)+Titaa2(i,j-1,k))+&
               Titaa2(i,j,k+1)+Titaa2(i,j,k-1))

            if (abs(Titaa1(i,j,k)) < 1e-10) Titaa1(i,j,k)=0

            Qvap1(i,j,k)=pro1*Qvap2(i,j,k)+&
               pro2*((Qvap2(i+1,j,k)+Qvap2(i-1,j,k))+&
               (Qvap2(i,j+1,k)+Qvap2(i,j-1,k))+&
               Qvap2(i,j,k+1)+Qvap2(i,j,k-1))

            if (abs(Qvap1(i,j,k)) < 1e-10) Qvap1(i,j,k)=0

            Qgot1(i,j,k)=pro1*Qgot2(i,j,k)+&
               pro2*((Qgot2(i+1,j,k)+Qgot2(i-1,j,k))+&
               (Qgot2(i,j+1,k)+Qgot2(i,j-1,k))+&
               Qgot2(i,j,k+1)+Qgot2(i,j,k-1))

            if (Qgot1(i,j,k) < 1e-10) Qgot1(i,j,k)=0

            Qllu1(i,j,k)=pro1*Qllu2(i,j,k)+&
               pro2*((Qllu2(i+1,j,k)+Qllu2(i-1,j,k))+&
               (Qllu2(i,j+1,k)+Qllu2(i,j-1,k))+Qllu2(i,j,k+1)+Qllu2(i,j,k-1))

            if (Qllu1(i,j,k) < 1e-10) Qllu1(i,j,k)=0

            Qcri1(i,j,k)=pro1*Qcri2(i,j,k)+&
               pro2*((Qcri2(i+1,j,k)+Qcri2(i-1,j,k))+&
               (Qcri2(i,j+1,k)+Qcri2(i,j-1,k))+Qcri2(i,j,k+1)+Qcri2(i,j,k-1))

            if (Qcri1(i,j,k) < 1e-10) Qcri1(i,j,k)=0

            Qnie1(i,j,k)=pro1*Qnie2(i,j,k)+&
               pro2*((Qnie2(i+1,j,k)+Qnie2(i-1,j,k))+&
               (Qnie2(i,j+1,k)+Qnie2(i,j-1,k))+&
               Qnie2(i,j,k+1)+Qnie2(i,j,k-1))

            if (Qnie1(i,j,k) < 1e-10) Qnie1(i,j,k)=0

            Qgra1(i,j,k)=pro1*Qgra2(i,j,k)+&
               pro2*((Qgra2(i+1,j,k)+Qgra2(i-1,j,k))+&
               (Qgra2(i,j+1,k)+Qgra2(i,j-1,k))+&
               Qgra2(i,j,k+1)+Qgra2(i,j,k-1))

            if (Qgra1(i,j,k) < 1e-10) Qgra1(i,j,k)=0

            aer1(i,j,k)=pro1*aer2(i,j,k)+&
               pro2*((aer2(i+1,j,k)+aer2(i-1,j,k))+&
               (aer2(i,j+1,k)+aer2(i,j-1,k))+aer2(i,j,k+1)+aer2(i,j,k-1))


            if (abs(aer1(i,j,k)) < 1e-10) aer1(i,j,k)=0
         end do
      end do
   end subroutine floor_condition_redefinition

   subroutine floor_and_ceiling_contour_redefinition()
      USE dimen, only: nx1, nz1
      USE perdim, only: Titaa1
      USE permic, only: Qvap1, aer1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1
      USE estbas, only: Qvap0, aer0
      implicit none
      integer :: i, j
      !contornos en el piso y en el techo
      do concurrent(i=1:nx1, j=1:nx1)
         Titaa1(i,j,0)=Titaa1(i,j,0)
         if (Titaa1(i,j,0) > 0.5) Titaa1(i,j,0)=.5
         if (-Titaa1(i,j,0) > 0.5) Titaa1(i,j,0)=-.5

         Titaa1(i,j,nz1)=Titaa1(i,j,nz1-1)

         !corregido para el vapor
         if (Qvap1(i,j,0) > Qvap0(0)*.5) then
            Qvap1(i,j,0)=.8*Qvap0(0)
         endif
         if (-Qvap1(i,j,0) > Qvap0(0)*.5) then
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
         if (-aer1(i,j,0) > 0.8*aer0(0)) then
            aer1(i,j,0)=-.8*aer0(0)
         endif
         aer1(i,j,nz1)=aer1(i,j,nz1-1)
      end do
   end subroutine floor_and_ceiling_contour_redefinition

   subroutine lateral_contour_redefinition()
      !contornos laterales
      USE dimen, only: nx1, nz1
      USE perdim, only: Titaa1
      USE permic, only: Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1
      implicit none
      integer :: j, k
      do concurrent(k=1:nz1-1, j=1:nx1)
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
      end do
   end subroutine lateral_contour_redefinition

   subroutine vapour_negative_correction()
      !correccion de negativos para el vapor
      USE dimen, only: nx1, nz1
      USE permic, only: Qvap1
      USE estbas, only: Qvap0
      implicit none
      integer :: i, j, k
      do concurrent(i=0:nx1+1, j=0:nx1+1, k=0:nz1)
         if(Qvap1(i,j,k)+Qvap0(k) < 0) then
            Qvap1(i,j,k)=-Qvap0(k)
         endif
      end do
   end subroutine vapour_negative_correction

   subroutine save_backup()
      USE model_var, only: tt, tte, posx, posy, Xnub, Ynub, posxx, posyy,&
         file_number, t1
      USE cant01, only: lte, ltg, ltb
      USE dimen, only: dt1
      USE config, only: output_directory
      USE perdim, only: W2, U2, V2, Fcalo, Titaa2, Titaa1, Pres2, Pres1, U1, V1,&
         W1
      USE io, only: str_gen
      USE permic, only: aer1, Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, Qvap2, aer2,&
         Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, Av, Vtgra0, Vtnie
      USE const, only: Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs, Eacrcn, Eautcn
      USE estbas, only: Den0, Qvap0, aer0, Tita0, Temp0, Pres00, aerrel, cc2,&
         Qvaprel, UU, VV
      USE model_initialization, only: cloud_position_init, cloud_movement_init, statistics_init
      implicit none
      integer :: unit_number
      if (tt/nint(lte/dt1)*nint(lte/dt1) == tt) then
         call statistics_init()
         tte=tte+1
         call cloud_position_init()
         call cloud_movement_init()

         open(newunit=unit_number,file=output_directory//"posnub"//'.sa', ACCESS="append")
         write(unit_number,*) tte,posx(tte),posy(tte),Xnub(tte),Ynub(tte),posxx,posyy
         close(unit_number)
      endif

      if (tt/nint(ltg/dt1)*nint(ltg/dt1) == tt) then
         file_number = str_gen(t1)
         !call graba231(k, W2, Titaa1, Qvap1, Qllu1, Qgra1, aer1, Qvap0, aer0)
         call graba320(U1, V1, W1, Titaa1, Pres1, Qvap1, Qgot1, Qllu1, Qcri1, Qnie1, Qgra1, aer1,file_number)
         t1=t1+1
      endif

      if (tt/nint(ltb/dt1)*nint(ltb/dt1) == tt) then
         call graba120(Den0,Temp0,Tita0,Pres00,Qvap0,cc2,aer0,UU,VV,&
            U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,&
            Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,aer1,aer2,Fcalo,&
            Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,Qvaprel,aerrel,Eautcn,Eacrcn)
      endif
   end subroutine save_backup
end module model_aux

