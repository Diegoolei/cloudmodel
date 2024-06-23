module model_aux
contains
   subroutine vapor_advection()
      !######################## Adveccion de vapores #######################
      use advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advaer1
      use microphysics_perturbation, only: spray_amt, vapor_amt, rain_amt, snow_amt, hail_amt
      use dinamic_var_perturbation, only: w_perturbed
      use dimensions, only: nx1
      implicit none
      integer :: i, j
      do concurrent (i=0:nx1+1, j=0:nx1+1)
         advvap1(i,j)=w_perturbed(i,j,1)*(vapor_amt(i,j,1)+vapor_amt(i,j,0))/4.
         advgot1(i,j)=0.
         advllu1(i,j)=w_perturbed(i,j,1)*rain_amt(i,j,1)
         if (w_perturbed(i,j,1) > 0) advllu1(i,j)=0.
         advaer1(i,j)=w_perturbed(i,j,1)*(spray_amt(i,j,1)+spray_amt(i,j,0))/4.
         if(w_perturbed(i,j,1) < 0) advaer1(i,j)=advaer1(i,j)*1.5
         advcri1(i,j)=0.
         advnie1(i,j)=w_perturbed(i,j,1)*snow_amt(i,j,1)
         if (w_perturbed(i,j,1) > 0) advnie1(i,j)=0.
         advgra1(i,j)=w_perturbed(i,j,1)*hail_amt(i,j,1)
         if (w_perturbed(i,j,1) > 0) advgra1(i,j)=0.
      end do
   end subroutine vapor_advection

   subroutine dinamics()
      !########### calculo de la dinamica y de la termodinamica ############
      use dimensions, only: dt1, dx1, nx1, nz1
      use model_var, only: dden0z, ener1, s, llluneg, lcrineg, lnieneg, lgraneg,&
         lvapneg, laerneg, Qvapneg, aerneg
      use estbas, only: Den0, Qvap0, aer0
      use dinamic_var_perturbation, only: w_perturbed, u_perturbed, v_perturbed, heat_force
      use advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advvap2, advgot2, advllu2, advcri2, advnie2, advgra2, advaer1,&
         advaer2
      use microphysics_perturbation, only: perturbed_drop_amt, perturbed_rain_amt,&
         perturbed_crystal_amt, perturbed_snow_amt, perturbed_hail_amt, perturbed_vapor_amt,&
         perturbed_spray_amt
      use lmngot, only: lgot, mgot, ngot
      use lmnllu, only: lllu, mllu, nllu
      use lmncri, only: lcri, mcri, ncri
      use lmnnie, only: lnie, mnie, nnie
      use lmngra, only: lgra, mgra, ngra
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
               ener1=.5*Den0(k)*(u_perturbed(i,j,k)**2.+v_perturbed(i,j,k)**2.+w_perturbed(i,j,k)**2.)+ener1

               !calculo de la temperatura potencial
               call tempot(l,m,n,dden0z,heat_force(i,j,k))
               heat_force(i,j,k)=0.

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
               if(perturbed_drop_amt(i,j,k) /= 0) then
                  if (i < lgot(1)) lgot(1)=i
                  if (i > lgot(2)) lgot(2)=i
                  if (j < mgot(1)) mgot(1)=j
                  if (j > mgot(2)) mgot(2)=j
                  if (k < ngot(1)) ngot(1)=k
                  if (k > ngot(2)) ngot(2)=k
                  s=s+1
               endif

               !limites de la lluvia
               if(perturbed_rain_amt(i,j,k) /= 0) then
                  if (i < lllu(1)) lllu(1)=i
                  if (i > lllu(2)) lllu(2)=i
                  if (j < mllu(1)) mllu(1)=j
                  if (j > mllu(2)) mllu(2)=j
                  if (k < nllu(1)) nllu(1)=k
                  if (k > nllu(2)) nllu(2)=k
                  llluneg=1
               endif

               !limites de los cristales
               if(perturbed_crystal_amt(i,j,k) /= 0) then
                  if (i < lcri(1)) lcri(1)=i
                  if (i > lcri(2)) lcri(2)=i
                  if (j < mcri(1)) mcri(1)=j
                  if (j > mcri(2)) mcri(2)=j
                  if (k < ncri(1)) ncri(1)=k
                  if (k > ncri(2)) ncri(2)=k
                  lcrineg=1
               endif

               !limites de la nieve
               if(perturbed_snow_amt(i,j,k) /= 0) then
                  if (i < lnie(1)) lnie(1)=i
                  if (i > lnie(2)) lnie(2)=i
                  if (j < mnie(1)) mnie(1)=j
                  if (j > mnie(2)) mnie(2)=j
                  if (k < nnie(1)) nnie(1)=k
                  if (k > nnie(2)) nnie(2)=k
                  lnieneg=1
               endif

               !limites del granizo
               if(perturbed_hail_amt(i,j,k) /= 0) then
                  if (i < lgra(1)) lgra(1)=i
                  if (i > lgra(2)) lgra(2)=i
                  if (j < mgra(1)) mgra(1)=j
                  if (j > mgra(2)) mgra(2)=j
                  if (k < ngra(1)) ngra(1)=k
                  if (k > ngra(2)) ngra(2)=k
                  lgraneg=1
               endif

               if(Qvap0(k)+perturbed_vapor_amt(i,j,k) < 0) then
                  Qvapneg=Qvapneg+Qvap0(k)+perturbed_vapor_amt(i,j,k)
                  lvapneg=1
               endif

               if(aer0(k)+perturbed_spray_amt(i,j,k) < 0) then
                  aerneg=aerneg+aer0(k)+perturbed_spray_amt(i,j,k)
                  laerneg=1
               endif
            end do
         end do
      end do

   end subroutine dinamics

   subroutine negative_correction
      use model_var, only: s, llluneg, lcrineg, lnieneg, lgraneg, lvapneg, laerneg,&
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
      use dimensions, only: nx1, nz1
      use model_var, only: vapt1, gott1, aert1
      use microphysics_perturbation, only: perturbed_vapor_amt, perturbed_drop_amt, perturbed_spray_amt
      use estbas, only: Qvap0
      implicit none
      integer :: i, j, k
      vapt1=0.
      gott1=0.
      aert1=0.
      do concurrent (i=1:nx1, j=1:nx1, k=1:nz1-1)
         vapt1=vapt1+perturbed_vapor_amt(i,j,k)
         gott1=gott1+perturbed_drop_amt(i,j,k)
         aert1=aert1+perturbed_spray_amt(i,j,k)

         if(perturbed_vapor_amt(i,j,k)+Qvap0(k) < 0) then
            stop
         endif
      end do
   end subroutine water_calculation

   subroutine microphisics_substring
      use dimensions, only: nx1, nz1, dt1, dx1
      use model_var, only: aux, P, T, Qvap, Naer, densi, Dv, iT, aux2, Vis,&
         esvs, elvs, Lvl, Lsl, Lvs, Eaccn, Eaucn, Eacng, nu, lll, current_time,&
         Qliq, e1, rl, rs, yy, daer, dqgot, dqcri, Taux, totnuc, vapt2, gott2,&
         aert2, qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, t2, Lsl00,&
         Fcal, daer2, totmic, vapt3, gott3, aert3, ener2, ener3, ener4, ener5,&
         qv, qg, daitot
      use dinamic_var_perturbation, only: thermal_property_2, pressure_perturbed,&
         ambient_temperature, heat_force
      use estbas, only: Qvap0, aer0, Tita0, Temp0, Pres00
      use cant01, only: ikapa, AA, lt2
      use constants, only: P00, Rd, Dv0, Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs,&
         Eacrcn, Eautcn, T0, Rv, G, Cp
      use microphysics_perturbation, only: perturbed_vapor_amt, perturbed_spray_amt,&
         perturbed_drop_amt, perturbed_rain_amt, perturbed_crystal_amt,&
         perturbed_snow_amt, perturbed_hail_amt
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
               aux=Pres00(k)+pressure_perturbed(i,j,k)
               P=aux**ikapa*P00
               T=(Tita0(k)+thermal_property_2(i,j,k))*aux
               ambient_temperature(i,j,k)=T-Temp0(k)
               Qvap=Qvap0(k)+perturbed_vapor_amt(i,j,k)
               Naer=aer0(k)+perturbed_spray_amt(i,j,k)
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
               lll = current_time

               Qliq=perturbed_drop_amt(i,j,k)
               e1=Qvap*Rv*T
               rl=(e1-elvs)/elvs
               rs=(e1-esvs)/esvs
               yy=0

               if ((rl > 1e-3 .or. rs > 1e-3).and.Naer > 0) then
                  call nuclea(Qvap,Qliq,Naer,T,densi,e1,elvs,esvs,rl,rs,Lvl,Lvs,daer,dqgot,dqcri)
                  Taux=T-Temp0(k)-ambient_temperature(i,j,k)
                  thermal_property_2(i,j,k)=T/aux-Tita0(k)
                  if (dqgot > 0) yy=1
               else
                  Taux=0.
                  dqgot=0.
                  dqcri=0.
                  daer=0.
               endif

               totnuc=totnuc+daer

               !segundo calculo de agua (sin laterales)
               vapt2=vapt2+perturbed_vapor_amt(i,j,k)
               gott2=gott2+perturbed_drop_amt(i,j,k)
               aert2=aert2+perturbed_spray_amt(i,j,k)
               if (perturbed_drop_amt(i,j,k) > 0 &
                  .or. dqgot > 0 &
                  .or. perturbed_rain_amt(i,j,k) > 0 &
                  .or. perturbed_crystal_amt(i,j,k) > 0 &
                  .or. perturbed_snow_amt(i,j,k) > 0) then

                  qgotaux=perturbed_drop_amt(i,j,k)
                  if (perturbed_drop_amt(i,j,k) == 0) qgotaux=0d0
                  qvapaux=perturbed_vapor_amt(i,j,k)+Qvap0(k)
                  qlluaux=perturbed_rain_amt(i,j,k)
                  if (perturbed_rain_amt(i,j,k) == 0) qlluaux=0d0
                  qcriaux=perturbed_crystal_amt(i,j,k)
                  if (perturbed_crystal_amt(i,j,k) == 0) then
                     qcriaux=0d0
                  endif
                  qnieaux=perturbed_snow_amt(i,j,k)
                  if (perturbed_snow_amt(i,j,k) == 0) qnieaux=0d0
                  qgraaux=perturbed_hail_amt(i,j,k)
                  if (perturbed_hail_amt(i,j,k) == 0) qgraaux=0d0
                  Naer=perturbed_spray_amt(i,j,k)+aer0(k)
                  T=ambient_temperature(i,j,k)+Temp0(k)
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
                     call microfis(elvs, esvs, Lvl, Lvs, Lsl, T, Dv, Eaccn,&
                        Eaucn, Eacng, Lsl00, Fcal, n, qvapaux, qgotaux, qlluaux,&
                        qcriaux,qnieaux, qgraaux, Naer, daer2, nu, yy)
                     heat_force(l,m,n)=heat_force(l,m,n)+Fcal/dt1/densi
                     Naer=Naer+daer2
                     totmic=totmic+daer2
                  end do

                  perturbed_drop_amt(i,j,k)=qgotaux
                  perturbed_rain_amt(i,j,k)=qlluaux
                  perturbed_crystal_amt(i,j,k)=qcriaux
                  perturbed_snow_amt(i,j,k)=qnieaux
                  perturbed_hail_amt(i,j,k)=qgraaux
                  perturbed_vapor_amt(i,j,k)=qvapaux-Qvap0(k)
                  perturbed_spray_amt(i,j,k)=Naer-aer0(k)
                  ambient_temperature(i,j,k)=T-Temp0(k)

               endif

               if (Tita0(k) < abs(thermal_property_2(i,j,k))+200.or.Temp0(k) < abs(ambient_temperature(i,j,k))+200) then
                  stop
               endif

               if(perturbed_spray_amt(i,j,k)+aer0(k) <= 0) then

                  if (perturbed_spray_amt(i,j,k)+aer0(k) < -aer0(k)*.05) then
                     stop
                  endif

                  perturbed_spray_amt(i,j,k)=-aer0(k)
               endif

               !tercer calculo de agua (sin laterales)
               vapt3=vapt3+perturbed_vapor_amt(i,j,k)
               gott3=gott3+perturbed_drop_amt(i,j,k)
               aert3=aert3+perturbed_spray_amt(i,j,k)

               !calculo de la energia
               ener2=densi*G*k*dx1+ener2
               ener3=densi*(Cp-Rd)*T+ener3
               ener4=P+ener4
               ener5=(perturbed_vapor_amt(i,j,k)+perturbed_drop_amt(i,j,k))*G*k*dx1+ener5
               qv=perturbed_vapor_amt(i,j,k)+qv
               qg=perturbed_drop_amt(i,j,k)+qg
               daitot=densi+daitot
            end do
         end do
      end do
   end subroutine microphisics_substring

   subroutine floor_and_ceiling_contour
      use dimensions, only: nx1, nz1, dt1, dx1
      use dinamic_var_perturbation, only: thermal_property_2, thermal_property_1,&
         w_perturbed, u_perturbed, v_perturbed
      use cant01, only: dx2
      use estbas, only: Tita0, Qvap0, aer0
      use model_var, only: auxx, auxy, auxz, turbu, aeraux, lapla, cks, Qvap,&
         qv, qg
      use microphysics_perturbation, only: vapor_amt, perturbed_vapor_amt,&
         perturbed_drop_amt, perturbed_crystal_amt, perturbed_rain_amt,&
         perturbed_snow_amt, perturbed_hail_amt, spray_amt, perturbed_spray_amt
      implicit none
      integer :: i, j

      Qvap = (qv+qg)/nx1**2.*nz1*.1
      !contornos en el piso y en el techo
      do concurrent (i=1:nx1, j=1:nx1)
         thermal_property_2(i,j,0) = thermal_property_1(i,j,0) -&
            w_perturbed(i,j,1)*(Tita0(0)+Tita0(1))*dt1/dx2
         thermal_property_2(i,j,nz1)=thermal_property_2(i,j,nz1-1)

         !suponemos que las velocidades horizontales a nivel de piso son
         !iguales a 1/4 de la correspondiente en el nivel 1
         !#TODO: Check this
         auxx = ((u_perturbed(i+1,j,1) + u_perturbed(i,j,1))*(vapor_amt(i+1,j,0) + vapor_amt(i,j,0))&
            -(u_perturbed(i-1,j,1)+u_perturbed(i,j,1))*(vapor_amt(i-1,j,0)+vapor_amt(i,j,0)))&
            /4.*.25

         auxy=((v_perturbed(i,j+1,1)+v_perturbed(i,j,1))*(vapor_amt(i,j+1,0)+vapor_amt(i,j,0))&
            -(v_perturbed(i,j-1,1)+v_perturbed(i,j,1))*(vapor_amt(i,j-1,0)+vapor_amt(i,j,0)))&
            /4.*.25

         auxz=w_perturbed(i,j,1)*((vapor_amt(i,j,1)+vapor_amt(i,j,0))+Qvap0(1)+Qvap0(0))/2.*.5

         perturbed_vapor_amt(i,j,nz1)=perturbed_vapor_amt(i,j,nz1-1)

         perturbed_drop_amt(i,j,0)=perturbed_drop_amt(i,j,1)
         perturbed_drop_amt(i,j,nz1)=perturbed_drop_amt(i,j,nz1-1)

         perturbed_crystal_amt(i,j,0)=perturbed_crystal_amt(i,j,1)
         perturbed_crystal_amt(i,j,nz1)=perturbed_crystal_amt(i,j,nz1-1)

         !Para que no se acumulen el piso
         perturbed_rain_amt(i,j,0)=perturbed_rain_amt(i,j,1)/2.
         perturbed_rain_amt(i,j,nz1)=perturbed_rain_amt(i,j,nz1-1)

         perturbed_snow_amt(i,j,0)=perturbed_snow_amt(i,j,1)
         perturbed_snow_amt(i,j,nz1)=perturbed_snow_amt(i,j,nz1-1)

         !Para que no se acumulen el piso
         perturbed_hail_amt(i,j,0)=perturbed_hail_amt(i,j,1)/2.
         perturbed_hail_amt(i,j,nz1)=perturbed_hail_amt(i,j,nz1-1)

         !suponemos que las velocidades horizontales a nivel de piso son
         !iguales a 1/4 de la correspondiente en el nivel 1

         auxx=((u_perturbed(i+1,j,1)+u_perturbed(i,j,1))*(spray_amt(i+1,j,0)+spray_amt(i,j,0))&
            -(u_perturbed(i-1,j,1)+u_perturbed(i,j,1))*(spray_amt(i-1,j,0)+spray_amt(i,j,0)))&
            /4.*.25
         auxy=((v_perturbed(i,j+1,1)+v_perturbed(i,j,1))*(spray_amt(i,j+1,0)+spray_amt(i,j,0))&
            -(v_perturbed(i,j-1,1)+v_perturbed(i,j,1))*(spray_amt(i,j-1,0)+spray_amt(i,j,0)))&
            /4.*.25
         auxz=w_perturbed(i,j,1)*((spray_amt(i,j,1)+spray_amt(i,j,0))+aer0(1)+aer0(0))/2.*.5

         if (w_perturbed(i,j,0) > 0) then
            aeraux=-((auxx+auxy)+2.*auxz)*dt1/dx1
         else
            !se refleja un 25 % de los aerosoles que caen
            aeraux=-((auxx+auxy)+.25*2.*auxz)*dt1/dx1
         endif

         !agregamos un termino de turbulencia para los aerosoles
         !a nivel de piso
         turbu = cks/dx1*.25 * (abs(u_perturbed(i,j,1))&
            + abs(v_perturbed(i,j,1))&
            + 2.*abs(w_perturbed(i,j,1)))

         lapla = ((spray_amt(i+1,j,0)+spray_amt(i-1,j,0))&
            + (spray_amt(i,j+1,0)+spray_amt(i,j-1,0)&
            + spray_amt(i,j,1)))&
            - 5.*spray_amt(i,j,0)

         lapla = lapla + (aer0(1) - aer0(0))

         perturbed_spray_amt(i,j,0)=aeraux+spray_amt(i,j,0)+turbu*lapla

         perturbed_spray_amt(i,j,nz1)=perturbed_spray_amt(i,j,nz1-1)
      end do
   end subroutine floor_and_ceiling_contour

   subroutine lateral_contour
      !contornos laterales
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: thermal_property_2
      use microphysics_perturbation, only: perturbed_vapor_amt, perturbed_drop_amt,&
         perturbed_rain_amt, perturbed_crystal_amt, perturbed_snow_amt,&
         perturbed_hail_amt, perturbed_spray_amt
      implicit none
      integer :: j, k
      do concurrent (k=1:nz1-1, j=1:nx1)
         thermal_property_2(0,j,k)=thermal_property_2(1,j,k)
         thermal_property_2(nx1+1,j,k)=thermal_property_2(nx1,j,k)
         thermal_property_2(j,0,k)=thermal_property_2(j,1,k)
         thermal_property_2(j,nx1+1,k)=thermal_property_2(j,nx1,k)
         perturbed_vapor_amt(0,j,k)=0.
         perturbed_vapor_amt(nx1+1,j,k)=0.
         perturbed_vapor_amt(j,0,k)=0.
         perturbed_vapor_amt(j,nx1+1,k)=0.
         perturbed_drop_amt(0,j,k)=perturbed_drop_amt(1,j,k)
         perturbed_drop_amt(nx1+1,j,k)=perturbed_drop_amt(nx1,j,k)
         perturbed_drop_amt(j,0,k)=perturbed_drop_amt(j,1,k)
         perturbed_drop_amt(j,nx1+1,k)=perturbed_drop_amt(j,nx1,k)
         perturbed_rain_amt(0,j,k)=0.
         perturbed_rain_amt(nx1+1,j,k)=0.
         perturbed_rain_amt(j,0,k)=0.
         perturbed_rain_amt(j,nx1+1,k)=0.
         perturbed_crystal_amt(0,j,k)=perturbed_crystal_amt(1,j,k)
         perturbed_crystal_amt(nx1+1,j,k)=perturbed_crystal_amt(nx1,j,k)
         perturbed_crystal_amt(j,0,k)=perturbed_crystal_amt(j,1,k)
         perturbed_crystal_amt(j,nx1+1,k)=perturbed_crystal_amt(j,nx1,k)
         perturbed_snow_amt(0,j,k)=perturbed_snow_amt(1,j,k)
         perturbed_snow_amt(nx1+1,j,k)=perturbed_snow_amt(nx1,j,k)
         perturbed_snow_amt(j,0,k)=perturbed_snow_amt(j,1,k)
         perturbed_snow_amt(j,nx1+1,k)=perturbed_snow_amt(j,nx1,k)
         perturbed_hail_amt(0,j,k)=0.
         perturbed_hail_amt(nx1+1,j,k)=0.
         perturbed_hail_amt(j,0,k)=0.
         perturbed_hail_amt(j,nx1+1,k)=0.
         perturbed_spray_amt(0,j,k)=perturbed_spray_amt(1,j,k)
         perturbed_spray_amt(nx1+1,j,k)=perturbed_spray_amt(nx1,j,k)
         perturbed_spray_amt(j,0,k)=perturbed_spray_amt(j,1,k)
         perturbed_spray_amt(j,nx1+1,k)=perturbed_spray_amt(j,nx1,k)
      end do
   end subroutine lateral_contour

   subroutine floor_condition_redefinition()
      use dimensions, only: nx1, nz1, dt1, dx1
      use dinamic_var_perturbation, only: thermal_property_2, thermal_property_1,&
         w_perturbed
      use cant01, only: pro3, pro4, pro1, pro2
      use microphysics_perturbation, only: vapor_amt, perturbed_vapor_amt,&
         perturbed_drop_amt, perturbed_rain_amt, perturbed_crystal_amt,&
         perturbed_snow_amt, perturbed_hail_amt, spray_amt, perturbed_spray_amt,&
         drop_amt, rain_amt, crystal_amt, snow_amt, hail_amt
      use estbas, only: aer0
      use model_var, only: aeraux
      implicit none
      integer :: i, j, k
      !*** modificada las condiciones en el piso
      !Redefinicion
      do concurrent (i=1:nx1, j=1:nx1)
         k=0
         thermal_property_1(i,j,k) = pro3*thermal_property_2(i,j,k)+&
            pro4*(&
            (thermal_property_2(i+1,j,k) + thermal_property_2(i-1,j,k)) +&
            (thermal_property_2(i,j+1,k) + thermal_property_2(i,j-1,k)))

         if (abs(thermal_property_1(i,j,k)) < 1e-10) thermal_property_1(i,j,k)=0

         vapor_amt(i,j,k) = pro3*perturbed_vapor_amt(i,j,k)&
            + pro4*(&
            (perturbed_vapor_amt(i+1,j,k) + perturbed_vapor_amt(i-1,j,k)) +&
            (perturbed_vapor_amt(i,j+1,k) + perturbed_vapor_amt(i,j-1,k)))


         if (abs(vapor_amt(i,j,k)) < 1e-10) vapor_amt(i,j,k)=0

         drop_amt(i,j,k) = pro3*perturbed_drop_amt(i,j,k) +&
            pro4*(&
            (perturbed_drop_amt(i+1,j,k) + perturbed_drop_amt(i-1,j,k)) +&
            (perturbed_drop_amt(i,j+1,k) + perturbed_drop_amt(i,j-1,k)))

         if (drop_amt(i,j,k) < 1e-10) drop_amt(i,j,k)=0

         rain_amt(i,j,k)=perturbed_rain_amt(i,j,k)

         if (rain_amt(i,j,k) < 1e-10) rain_amt(i,j,k)=0

         crystal_amt(i,j,k) = pro3*perturbed_crystal_amt(i,j,k) +&
            pro4*(&
            (perturbed_crystal_amt(i+1,j,k) + perturbed_crystal_amt(i-1,j,k)) +&
            (perturbed_crystal_amt(i,j+1,k) + perturbed_crystal_amt(i,j-1,k)))

         if (crystal_amt(i,j,k) < 1e-10) crystal_amt(i,j,k)=0

         snow_amt(i,j,k) = pro3*perturbed_snow_amt(i,j,k)+&
            pro4*(&
            (perturbed_snow_amt(i+1,j,k) + perturbed_snow_amt(i-1,j,k)) +&
            (perturbed_snow_amt(i,j+1,k) + perturbed_snow_amt(i,j-1,k)))

         if (snow_amt(i,j,k) < 1e-10) snow_amt(i,j,k)=0

         hail_amt(i,j,k)=perturbed_hail_amt(i,j,k)

         if (hail_amt(i,j,k) < 1e-10) hail_amt(i,j,k)=0

         spray_amt(i,j,k) = pro3*perturbed_spray_amt(i,j,k) +&
            pro4*(&
            (perturbed_spray_amt(i+1,j,k) + perturbed_spray_amt(i-1,j,k))+&
            (perturbed_spray_amt(i,j+1,k) + perturbed_spray_amt(i,j-1,k)))

         !correccion cambiando la absorcion de aerosoles
         if ((rain_amt(i,j,1)+hail_amt(i,j,1)) > 1e-6 .and.w_perturbed(i,j,1) < 0) then
            aeraux=-w_perturbed(i,j,1)*.5*dt1/(dx1/2)
            spray_amt(i,j,k)=spray_amt(i,j,k)-(spray_amt(i,j,k)+aer0(k))*aeraux
         endif

         if (abs(spray_amt(i,j,k)) < 1e-10) spray_amt(i,j,k)=0


         do concurrent(k=1:nz1-1)
            thermal_property_1(i,j,k) = pro1*thermal_property_2(i,j,k) +&
               pro2*(&
               (thermal_property_2(i+1,j,k) + thermal_property_2(i-1,j,k))+&
               (thermal_property_2(i,j+1,k) + thermal_property_2(i,j-1,k))+&
               thermal_property_2(i,j,k+1) +&
               thermal_property_2(i,j,k-1))

            if (abs(thermal_property_1(i,j,k)) < 1e-10) thermal_property_1(i,j,k)=0

            vapor_amt(i,j,k) = pro1*perturbed_vapor_amt(i,j,k) +&
               pro2*(&
               (perturbed_vapor_amt(i+1,j,k) + perturbed_vapor_amt(i-1,j,k)) +&
               (perturbed_vapor_amt(i,j+1,k) + perturbed_vapor_amt(i,j-1,k)) +&
               perturbed_vapor_amt(i,j,k+1) +&
               perturbed_vapor_amt(i,j,k-1))

            if (abs(vapor_amt(i,j,k)) < 1e-10) vapor_amt(i,j,k)=0

            drop_amt(i,j,k) = pro1*perturbed_drop_amt(i,j,k) +&
               pro2*(&
               (perturbed_drop_amt(i+1,j,k) + perturbed_drop_amt(i-1,j,k)) +&
               (perturbed_drop_amt(i,j+1,k) + perturbed_drop_amt(i,j-1,k)) +&
               perturbed_drop_amt(i,j,k+1) +&
               perturbed_drop_amt(i,j,k-1))

            if (drop_amt(i,j,k) < 1e-10) drop_amt(i,j,k)=0

            rain_amt(i,j,k) = pro1*perturbed_rain_amt(i,j,k) +&
               pro2*(&
               (perturbed_rain_amt(i+1,j,k) + perturbed_rain_amt(i-1,j,k)) +&
               (perturbed_rain_amt(i,j+1,k) + perturbed_rain_amt(i,j-1,k)) +&
               perturbed_rain_amt(i,j,k+1) +&
               perturbed_rain_amt(i,j,k-1))

            if (rain_amt(i,j,k) < 1e-10) rain_amt(i,j,k)=0

            crystal_amt(i,j,k) = pro1*perturbed_crystal_amt(i,j,k) +&
               pro2*(&
               (perturbed_crystal_amt(i+1,j,k) + perturbed_crystal_amt(i-1,j,k)) +&
               (perturbed_crystal_amt(i,j+1,k) + perturbed_crystal_amt(i,j-1,k)) +&
               perturbed_crystal_amt(i,j,k+1) +&
               perturbed_crystal_amt(i,j,k-1))

            if (crystal_amt(i,j,k) < 1e-10) crystal_amt(i,j,k)=0

            snow_amt(i,j,k) = pro1*perturbed_snow_amt(i,j,k) +&
               pro2*(&
               (perturbed_snow_amt(i+1,j,k) + perturbed_snow_amt(i-1,j,k)) +&
               (perturbed_snow_amt(i,j+1,k) + perturbed_snow_amt(i,j-1,k)) +&
               perturbed_snow_amt(i,j,k+1) +&
               perturbed_snow_amt(i,j,k-1))

            if (snow_amt(i,j,k) < 1e-10) snow_amt(i,j,k)=0

            hail_amt(i,j,k) = pro1*perturbed_hail_amt(i,j,k) +&
               pro2*(&
               (perturbed_hail_amt(i+1,j,k) + perturbed_hail_amt(i-1,j,k)) +&
               (perturbed_hail_amt(i,j+1,k) + perturbed_hail_amt(i,j-1,k)) +&
               perturbed_hail_amt(i,j,k+1) +&
               perturbed_hail_amt(i,j,k-1))

            if (hail_amt(i,j,k) < 1e-10) hail_amt(i,j,k)=0

            spray_amt(i,j,k) = pro1*perturbed_spray_amt(i,j,k) +&
               pro2*(&
               (perturbed_spray_amt(i+1,j,k) + perturbed_spray_amt(i-1,j,k)) +&
               (perturbed_spray_amt(i,j+1,k) + perturbed_spray_amt(i,j-1,k)) +&
               perturbed_spray_amt(i,j,k+1) +&
               perturbed_spray_amt(i,j,k-1))


            if (abs(spray_amt(i,j,k)) < 1e-10) spray_amt(i,j,k)=0
         end do
      end do
   end subroutine floor_condition_redefinition

   subroutine floor_and_ceiling_contour_redefinition()
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: thermal_property_1
      use microphysics_perturbation, only: vapor_amt, spray_amt, drop_amt,&
         rain_amt, crystal_amt, snow_amt, hail_amt
      use estbas, only: Qvap0, aer0
      implicit none
      integer :: i, j
      !contornos en el piso y en el techo
      do concurrent(i=1:nx1, j=1:nx1)
         thermal_property_1(i,j,0)=thermal_property_1(i,j,0)
         if (thermal_property_1(i,j,0) > 0.5) thermal_property_1(i,j,0)=.5
         if (-thermal_property_1(i,j,0) > 0.5) thermal_property_1(i,j,0)=-.5

         thermal_property_1(i,j,nz1)=thermal_property_1(i,j,nz1-1)

         !corregido para el vapor
         if (vapor_amt(i,j,0) > Qvap0(0)*.5) then
            vapor_amt(i,j,0)=.8*Qvap0(0)
         endif
         if (-vapor_amt(i,j,0) > Qvap0(0)*.5) then
            vapor_amt(i,j,0)=-.8*Qvap0(0)
         endif

         vapor_amt(i,j,nz1)=vapor_amt(i,j,nz1-1)
         drop_amt(i,j,0)=0.
         drop_amt(i,j,nz1)=drop_amt(i,j,nz1-1)
         rain_amt(i,j,0)=rain_amt(i,j,0)
         rain_amt(i,j,nz1)=rain_amt(i,j,nz1-1)
         crystal_amt(i,j,0)=0.
         crystal_amt(i,j,nz1)=crystal_amt(i,j,nz1-1)

         snow_amt(i,j,0)=0.
         snow_amt(i,j,nz1)=snow_amt(i,j,nz1-1)

         hail_amt(i,j,0)=hail_amt(i,j,0)
         hail_amt(i,j,nz1)=hail_amt(i,j,nz1-1)

         !corregido para los aerosoles
         if (-spray_amt(i,j,0) > 0.8*aer0(0)) then
            spray_amt(i,j,0)=-.8*aer0(0)
         endif
         spray_amt(i,j,nz1)=spray_amt(i,j,nz1-1)
      end do
   end subroutine floor_and_ceiling_contour_redefinition

   subroutine lateral_contour_redefinition()
      !contornos laterales
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: thermal_property_1
      use microphysics_perturbation, only: vapor_amt, drop_amt, rain_amt,&
         crystal_amt, snow_amt, hail_amt, spray_amt
      implicit none
      integer :: j, k
      do concurrent(k=1:nz1-1, j=1:nx1)
         thermal_property_1(0,j,k)=thermal_property_1(1,j,k)
         thermal_property_1(nx1+1,j,k)=thermal_property_1(nx1,j,k)
         thermal_property_1(j,0,k)=thermal_property_1(j,1,k)
         thermal_property_1(j,nx1+1,k)=thermal_property_1(j,nx1,k)
         vapor_amt(0,j,k)=0.
         vapor_amt(nx1+1,j,k)=0.
         vapor_amt(j,0,k)=0.
         vapor_amt(j,nx1+1,k)=0.
         drop_amt(0,j,k)=drop_amt(1,j,k)
         drop_amt(nx1+1,j,k)=drop_amt(nx1,j,k)
         drop_amt(j,0,k)=drop_amt(j,1,k)
         drop_amt(j,nx1+1,k)=drop_amt(j,nx1,k)
         rain_amt(0,j,k)=0.
         rain_amt(nx1+1,j,k)=0.
         rain_amt(j,0,k)=0.
         rain_amt(j,nx1+1,k)=0.
         crystal_amt(0,j,k)=crystal_amt(1,j,k)
         crystal_amt(nx1+1,j,k)=crystal_amt(nx1,j,k)
         crystal_amt(j,0,k)=crystal_amt(j,1,k)
         crystal_amt(j,nx1+1,k)=crystal_amt(j,nx1,k)
         snow_amt(0,j,k)=snow_amt(1,j,k)
         snow_amt(nx1+1,j,k)=snow_amt(nx1,j,k)
         snow_amt(j,0,k)=snow_amt(j,1,k)
         snow_amt(j,nx1+1,k)=snow_amt(j,nx1,k)
         hail_amt(0,j,k)=0.
         hail_amt(nx1+1,j,k)=0.
         hail_amt(j,0,k)=0.
         hail_amt(j,nx1+1,k)=0.

         spray_amt(0,j,k)=spray_amt(1,j,k)
         spray_amt(nx1+1,j,k)=spray_amt(nx1,j,k)
         spray_amt(j,0,k)=spray_amt(j,1,k)
         spray_amt(j,nx1+1,k)=spray_amt(j,nx1,k)
      end do
   end subroutine lateral_contour_redefinition

   subroutine vapour_negative_correction()
      !correccion de negativos para el vapor
      use dimensions, only: nx1, nz1
      use microphysics_perturbation, only: vapor_amt
      use estbas, only: Qvap0
      implicit none
      integer :: i, j, k
      do concurrent(i=0:nx1+1, j=0:nx1+1, k=0:nz1)
         if(vapor_amt(i,j,k)+Qvap0(k) < 0) then
            vapor_amt(i,j,k)=-Qvap0(k)
         endif
      end do
   end subroutine vapour_negative_correction

   subroutine save_backup()
      use model_var, only: current_time, tte, posx, posy, Xnub, Ynub, posxx, posyy,&
         file_number, t1
      use cant01, only: lte, ltg, ltb
      use dimensions, only: dt1
      use config, only: output_directory
      use dinamic_var_perturbation, only: w_perturbed, u_perturbed, v_perturbed,&
         heat_force, thermal_property_2, thermal_property_1, pressure_perturbed,&
         pressure_original, u_original, v_original, w_original
      use io, only: str_gen
      use microphysics_perturbation, only: spray_amt, perturbed_drop_amt,&
         perturbed_rain_amt, perturbed_crystal_amt, perturbed_snow_amt,&
         perturbed_hail_amt, perturbed_vapor_amt, perturbed_spray_amt,&
         vapor_amt, drop_amt, rain_amt, crystal_amt, snow_amt, hail_amt,&
         Av, Vtgra0, Vtnie
      use constants, only: Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs, Eacrcn, Eautcn
      use estbas, only: Den0, Qvap0, aer0, Tita0, Temp0, Pres00, aerrel, cc2,&
         Qvaprel, UU, VV
      use model_initialization, only: cloud_position_init, cloud_movement_init,&
         statistics_init
      implicit none
      integer :: unit_number
      if (current_time/nint(lte/dt1)*nint(lte/dt1) == current_time) then
         call statistics_init()
         tte=tte+1
         call cloud_position_init()
         call cloud_movement_init()

         open(newunit=unit_number,file=output_directory//"posnub"//'.sa',&
            ACCESS="append")
         write(unit_number,*) tte,posx(tte),posy(tte),Xnub(tte),Ynub(tte),posxx,posyy
         close(unit_number)
      endif

      if (current_time/nint(ltg/dt1)*nint(ltg/dt1) == current_time) then
         file_number = str_gen(t1)
         !call graba231(k, w_perturbed, thermal_property_1, vapor_amt, rain_amt, hail_amt, spray_amt, Qvap0, aer0)
         call graba320(u_original, v_original, w_original, thermal_property_1,&
            pressure_original, vapor_amt, drop_amt, rain_amt, crystal_amt,&
            snow_amt, hail_amt, spray_amt, file_number)
         t1=t1+1
      endif

      if (current_time/nint(ltb/dt1)*nint(ltb/dt1) == current_time) then
         call graba120(Den0, Temp0, Tita0, Pres00, Qvap0, cc2, aer0, UU, VV,&
            u_original, u_perturbed, v_original, v_perturbed, w_original,&
            w_perturbed, thermal_property_1, thermal_property_2, pressure_original,&
            pressure_perturbed, vapor_amt, perturbed_vapor_amt, drop_amt,&
            perturbed_drop_amt, rain_amt, perturbed_rain_amt, crystal_amt,&
            perturbed_crystal_amt, snow_amt, perturbed_snow_amt, hail_amt,&
            perturbed_hail_amt, spray_amt, perturbed_spray_amt, heat_force,&
            Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, Vtgra0, Qvaprel,&
            aerrel, Eautcn, Eacrcn)
      endif
   end subroutine save_backup
end module model_aux

