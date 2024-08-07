module model_aux
contains
   subroutine vapor_advection()
      !!Adveccion de vapores
      use advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1, &
                        advcri1, advaer1
      use microphysics_perturbation, only: aerosol_base, vapor_base, rain_base, snow_base, hail_base
      use dinamic_var_perturbation, only: w_perturbed_new
      use dimensions, only: nx1
      implicit none
      integer :: i, j
      do concurrent(i=0:nx1 + 1, j=0:nx1 + 1)
         advvap1(i, j) = w_perturbed_new(i, j, 1)*(vapor_base(i, j, 1) + vapor_base(i, j, 0))/4.
         advgot1(i, j) = 0.
         advllu1(i, j) = w_perturbed_new(i, j, 1)*rain_base(i, j, 1)
         if (w_perturbed_new(i, j, 1) > 0) advllu1(i, j) = 0.
         advaer1(i, j) = w_perturbed_new(i, j, 1)*(aerosol_base(i, j, 1) + aerosol_base(i, j, 0))/4.
         if (w_perturbed_new(i, j, 1) < 0) advaer1(i, j) = advaer1(i, j)*1.5
         advcri1(i, j) = 0.
         advnie1(i, j) = w_perturbed_new(i, j, 1)*snow_base(i, j, 1)
         if (w_perturbed_new(i, j, 1) > 0) advnie1(i, j) = 0.
         advgra1(i, j) = w_perturbed_new(i, j, 1)*hail_base(i, j, 1)
         if (w_perturbed_new(i, j, 1) > 0) advgra1(i, j) = 0.
      end do
   end subroutine vapor_advection

   subroutine dinamics()
      !! calculo de la dinamica y de la termodinamica
      use dimensions, only: dt1, dx1, nx1, nz1
      use model_var, only: dden0z, ener1, s, llluneg, lcrineg, lnieneg, lgraneg, &
                           lvapneg, laerneg, Qvapneg, aerneg
      use initial_z_state, only: air_density_z_initial, vapor_z_initial, aerosol_z_initial
      use dinamic_var_perturbation, only: w_perturbed_new, u_perturbed_new, v_perturbed_new, heat_force
      use advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1, &
                        advcri1, advvap2, advgot2, advllu2, advcri2, advnie2, advgra2, advaer1, &
                        advaer2
      use microphysics_perturbation, only: drop_new, rain_new, &
                                           crystal_new, snow_new, hail_new, vapor_new, &
                                           aerosol_new
      use lmngot, only: lgot, mgot, ngot
      use lmnllu, only: lllu, mllu, nllu
      use lmncri, only: lcri, mcri, ncri
      use lmnnie, only: lnie, mnie, nnie
      use lmngra, only: lgra, mgra, ngra
      use extra_subrut, only: turbu1, turbu2, inhomogeneous_velocities, &
                              tempot, dvapor, dgotit, dlluvi, dcrist, dnieve, dgrani, daeros
      implicit none
      integer :: i, j, k, l, m, n
      s = 0
      Qvapneg = 0.
      lvapneg = 0
      aerneg = 0.
      laerneg = 0
      llluneg = 0
      lcrineg = 0.
      lnieneg = 0.
      lgraneg = 0.
      ener1 = 0.
      lgot(1) = nx1
      lgot(2) = 0
      mgot(1) = nx1
      mgot(2) = 0
      ngot(1) = nz1
      ngot(2) = 0
      lllu(1) = nx1
      lllu(2) = 0
      mllu(1) = nx1
      mllu(2) = 0
      nllu(1) = nz1
      nllu(2) = 0
      lcri(1) = nx1
      lcri(2) = 0
      mcri(1) = nx1
      mcri(2) = 0
      ncri(1) = nz1
      ncri(2) = 0
      lnie(1) = nx1
      lnie(2) = 0
      mnie(1) = nx1
      mnie(2) = 0
      nnie(1) = nz1
      nnie(2) = 0
      lgra(1) = nx1
      lgra(2) = 0
      mgra(1) = nx1
      mgra(2) = 0
      ngra(1) = nz1
      ngra(2) = 0
      do k = 1, nz1 - 1
         n = k
         dden0z = (air_density_z_initial(k + 1) - air_density_z_initial(k - 1))/air_density_z_initial(k)
         call turbu1(n)
         do i = 1, nx1
            l = i
            do j = 1, nx1
               m = j
               !calculo del coeficiente de turbulencia y derivadas
               call turbu2(l, m)

               !calculo de las inhomogeneidades para las velocidades
               call inhomogeneous_velocities(l, m, n, dden0z)

               !calculo de la energia cinetica
               ener1 = .5*air_density_z_initial(k)*( &
                       u_perturbed_new(i, j, k)**2.+ &
                       v_perturbed_new(i, j, k)**2.+ &
                       w_perturbed_new(i, j, k)**2.) + ener1

               !calculo de la temperatura potencial
               call tempot(l, m, n, dden0z, heat_force(i, j, k))
               heat_force(i, j, k) = 0.

               !dinamica del vapor y de las gotitas
               call dvapor(l, m, n)
               advvap1(i, j) = advvap2(i, j)
               call dgotit(l, m, n)
               advgot1(i, j) = advgot2(i, j)
               call dlluvi(l, m, n)
               advllu1(i, j) = advllu2(i, j)
               call dcrist(l, m, n)
               advcri1(i, j) = advcri2(i, j)
               call dnieve(l, m, n)
               advnie1(i, j) = advnie2(i, j)
               call dgrani(l, m, n)
               advgra1(i, j) = advgra2(i, j)
               call daeros(l, m, n)
               advaer1(i, j) = advaer2(i, j)

               !limites de la nube
               if (drop_new(i, j, k) /= 0) then
                  if (i < lgot(1)) lgot(1) = i
                  if (i > lgot(2)) lgot(2) = i
                  if (j < mgot(1)) mgot(1) = j
                  if (j > mgot(2)) mgot(2) = j
                  if (k < ngot(1)) ngot(1) = k
                  if (k > ngot(2)) ngot(2) = k
                  s = s + 1
               end if

               !limites de la lluvia
               if (rain_new(i, j, k) /= 0) then
                  if (i < lllu(1)) lllu(1) = i
                  if (i > lllu(2)) lllu(2) = i
                  if (j < mllu(1)) mllu(1) = j
                  if (j > mllu(2)) mllu(2) = j
                  if (k < nllu(1)) nllu(1) = k
                  if (k > nllu(2)) nllu(2) = k
                  llluneg = 1
               end if

               !limites de los cristales
               if (crystal_new(i, j, k) /= 0) then
                  if (i < lcri(1)) lcri(1) = i
                  if (i > lcri(2)) lcri(2) = i
                  if (j < mcri(1)) mcri(1) = j
                  if (j > mcri(2)) mcri(2) = j
                  if (k < ncri(1)) ncri(1) = k
                  if (k > ncri(2)) ncri(2) = k
                  lcrineg = 1
               end if

               !limites de la nieve
               if (snow_new(i, j, k) /= 0) then
                  if (i < lnie(1)) lnie(1) = i
                  if (i > lnie(2)) lnie(2) = i
                  if (j < mnie(1)) mnie(1) = j
                  if (j > mnie(2)) mnie(2) = j
                  if (k < nnie(1)) nnie(1) = k
                  if (k > nnie(2)) nnie(2) = k
                  lnieneg = 1
               end if

               !limites del granizo
               if (hail_new(i, j, k) /= 0) then
                  if (i < lgra(1)) lgra(1) = i
                  if (i > lgra(2)) lgra(2) = i
                  if (j < mgra(1)) mgra(1) = j
                  if (j > mgra(2)) mgra(2) = j
                  if (k < ngra(1)) ngra(1) = k
                  if (k > ngra(2)) ngra(2) = k
                  lgraneg = 1
               end if

               if (vapor_z_initial(k) + vapor_new(i, j, k) < 0) then
                  Qvapneg = Qvapneg + vapor_z_initial(k) + vapor_new(i, j, k)
                  lvapneg = 1
               end if

               if (aerosol_z_initial(k) + aerosol_new(i, j, k) < 0) then
                  aerneg = aerneg + aerosol_z_initial(k) + aerosol_new(i, j, k)
                  laerneg = 1
               end if
            end do
         end do
      end do

   end subroutine dinamics

   subroutine negative_correction
      use model_var, only: s, llluneg, lcrineg, lnieneg, lgraneg, lvapneg, laerneg, &
                           Qvapneg, aerneg
      use extra_subrut, only: corgot, corllu, corcri, cornie, corgra, corvap, coraer
      implicit none
      if (s >= 1) call corgot
      if (llluneg == 1) call corllu
      if (lcrineg == 1) call corcri
      if (lnieneg == 1) call cornie
      if (lgraneg == 1) call corgra
      if (lvapneg == 1) call corvap(Qvapneg)
      if (laerneg == 1) call coraer(aerneg)
   end subroutine negative_correction

   subroutine water_calculation
      !! primer calculo de agua (sin laterales)
      use dimensions, only: nx1, nz1
      use model_var, only: vapt1, gott1, aert1
      use microphysics_perturbation, only: vapor_new, drop_new, aerosol_new
      use initial_z_state, only: vapor_z_initial
      implicit none
      integer :: i, j, k
      vapt1 = 0.
      gott1 = 0.
      aert1 = 0.
      do concurrent(i=1:nx1, j=1:nx1, k=1:nz1 - 1)
         vapt1 = vapt1 + vapor_new(i, j, k)
         gott1 = gott1 + drop_new(i, j, k)
         aert1 = aert1 + aerosol_new(i, j, k)
      end do
   end subroutine water_calculation

   subroutine microphisics_substring
      !! Sublazo Microfisico
      use dimensions, only: nx1, nz1, dt1, dx1
      use model_var, only: aux, P, T, Qvap, Naer, densi, Dv, iT, aux2, Vis, &
                           esvs, elvs, Lvl, Lsl, Lvs, Eaccn, Eaucn, Eacng, nu, lll, current_time, &
                           Qliq, e1, rl, rs, yy, daer, dqgot, dqcri, Taux, totnuc, vapt2, gott2, &
                           aert2, qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, t2, Lsl00, &
                           Fcal, daer2, totmic, vapt3, gott3, aert3, ener2, ener3, ener4, ener5, &
                           qv, qg, daitot
      use dinamic_var_perturbation, only: theta_new, pressure_new, &
                                          temperature, heat_force
      use initial_z_state, only: vapor_z_initial, aerosol_z_initial, theta_z_initial, temperature_z_initial, Pres00
      use cant01, only: ikapa, AA, lt2
      use constants, only: P00, Rd, Dv0, Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs, &
                           Eacrcn, Eautcn, T0, Rv, G, Cp
      use microphysics_perturbation, only: vapor_new, aerosol_new, &
                                           drop_new, rain_new, crystal_new, &
                                           snow_new, hail_new
      use extra_subrut, only: nuclea
      use microphysics, only: microfis
      implicit none
      integer :: k, n, i, l, j, m
      totnuc = 0.
      vapt2 = 0.
      gott2 = 0.
      aert2 = 0.
      totmic = 0.
      vapt3 = 0.
      gott3 = 0.
      aert3 = 0.
      ener2 = 0.
      ener3 = 0.
      ener4 = 0.
      ener5 = 0.
      qv = 0.
      qg = 0.
      daitot = 0.
      do k = 1, nz1 - 1
         n = k
         do i = 1, nx1
            l = i
            do j = 1, nx1
               m = j
               !calculo de T,P,Densi,Dv,Vis
               aux = Pres00(k) + pressure_new(i, j, k)
               P = aux**ikapa*P00
               T = (theta_z_initial(k) + theta_new(i, j, k))*aux
               temperature(i, j, k) = T - temperature_z_initial(k)
               Qvap = vapor_z_initial(k) + vapor_new(i, j, k)
               Naer = aerosol_z_initial(k) + aerosol_new(i, j, k)
               densi = P/T/Rd - AA*Qvap
               Dv = Dv0*(T/273.15)**1.94*(P00/P)

               !calculo de Vis, Lvl, Lsl, Lvs, elvs y  esvs
               iT = int(T)
               aux2 = T - iT
               Vis = Tvis(iT)*(1 - aux2) + Tvis(iT + 1)*aux2
               elvs = Telvs(iT)*(1 - aux2) + Telvs(iT + 1)*aux2
               esvs = Tesvs(iT)*(1 - aux2) + Tesvs(iT + 1)*aux2
               Lvl = Tlvl(iT)*(1 - aux2) + Tlvl(iT + 1)*aux2
               Lsl = Tlsl(iT)*(1 - aux2) + Tlsl(iT + 1)*aux2
               Lvs = Tlvs(iT)*(1 - aux2) + Tlvs(iT + 1)*aux2
               Eaccn = Eacrcn(iT)*(1 - aux2) + Eacrcn(iT + 1)*aux2
               Eaucn = Eautcn(iT)*(1 - aux2) + Eautcn(iT + 1)*aux2
               if (T >= T0) then
                  Eacng = 1.
               else
                  Eacng = exp(.08*(T - T0))
               end if

               nu = Vis/densi

               !nucleacion (de ser necesario tiene otro paso de tiempo)
               lll = current_time

               Qliq = drop_new(i, j, k)
               e1 = Qvap*Rv*T
               rl = (e1 - elvs)/elvs
               rs = (e1 - esvs)/esvs
               yy = 0

               if ((rl > 1e-3 .or. rs > 1e-3) .and. Naer > 0) then
                  call nuclea(Qvap, Qliq, Naer, T, densi, e1, elvs, esvs, rl, rs, Lvl, Lvs, daer, dqgot, dqcri)
                  Taux = T - temperature_z_initial(k) - temperature(i, j, k)
                  theta_new(i, j, k) = T/aux - theta_z_initial(k)
                  if (dqgot > 0) yy = 1
               else
                  Taux = 0.
                  dqgot = 0.
                  dqcri = 0.
                  daer = 0.
               end if

               totnuc = totnuc + daer

               !segundo calculo de agua (sin laterales)
               vapt2 = vapt2 + vapor_new(i, j, k)
               gott2 = gott2 + drop_new(i, j, k)
               aert2 = aert2 + aerosol_new(i, j, k)
               if (drop_new(i, j, k) > 0 &
                   .or. dqgot > 0 &
                   .or. rain_new(i, j, k) > 0 &
                   .or. crystal_new(i, j, k) > 0 &
                   .or. snow_new(i, j, k) > 0) then

                  qgotaux = drop_new(i, j, k)
                  if (drop_new(i, j, k) == 0) qgotaux = 0d0
                  qvapaux = vapor_new(i, j, k) + vapor_z_initial(k)
                  qlluaux = rain_new(i, j, k)
                  if (rain_new(i, j, k) == 0) qlluaux = 0d0
                  qcriaux = crystal_new(i, j, k)
                  if (crystal_new(i, j, k) == 0) then
                     qcriaux = 0d0
                  end if
                  qnieaux = snow_new(i, j, k)
                  if (snow_new(i, j, k) == 0) qnieaux = 0d0
                  qgraaux = hail_new(i, j, k)
                  if (hail_new(i, j, k) == 0) qgraaux = 0d0
                  Naer = aerosol_new(i, j, k) + aerosol_z_initial(k)
                  T = temperature(i, j, k) + temperature_z_initial(k)
                  do t2 = 1, lt2
                     qgotaux = qgotaux + dqgot/float(lt2)
                     qcriaux = qcriaux + dqcri/float(lt2)
                     qvapaux = qvapaux - (dqgot + dqcri)/float(lt2)
                     Naer = Naer + daer/float(lt2)
                     T = T + Taux/float(lt2)

                     !calculo de elvs y esvs
                     iT = int(T)
                     aux2 = T - iT
                     elvs = Telvs(iT)*(1 - aux2) + Telvs(iT + 1)*aux2
                     esvs = Tesvs(iT)*(1 - aux2) + Tesvs(iT + 1)*aux2
                     call microfis(elvs, esvs, Lvl, Lvs, Lsl, T, Dv, Eaccn, &
                                   Eaucn, Eacng, Lsl00, Fcal, n, qvapaux, qgotaux, qlluaux, &
                                   qcriaux, qnieaux, qgraaux, Naer, daer2, nu, yy)
                     heat_force(l, m, n) = heat_force(l, m, n) + Fcal/dt1/densi
                     Naer = Naer + daer2
                     totmic = totmic + daer2
                  end do

                  drop_new(i, j, k) = qgotaux
                  rain_new(i, j, k) = qlluaux
                  crystal_new(i, j, k) = qcriaux
                  snow_new(i, j, k) = qnieaux
                  hail_new(i, j, k) = qgraaux
                  vapor_new(i, j, k) = qvapaux - vapor_z_initial(k)
                  aerosol_new(i, j, k) = Naer - aerosol_z_initial(k)
                  temperature(i, j, k) = T - temperature_z_initial(k)

               end if
               if (aerosol_new(i, j, k) + aerosol_z_initial(k) <= 0) then
                  aerosol_new(i, j, k) = -aerosol_z_initial(k)
               end if

               !tercer calculo de agua (sin laterales)
               vapt3 = vapt3 + vapor_new(i, j, k)
               gott3 = gott3 + drop_new(i, j, k)
               aert3 = aert3 + aerosol_new(i, j, k)

               !calculo de la energia
               ener2 = densi*G*k*dx1 + ener2
               ener3 = densi*(Cp - Rd)*T + ener3
               ener4 = P + ener4
               ener5 = (vapor_new(i, j, k) + drop_new(i, j, k))*G*k*dx1 + ener5
               qv = vapor_new(i, j, k) + qv
               qg = drop_new(i, j, k) + qg
               daitot = densi + daitot
            end do
         end do
      end do
   end subroutine microphisics_substring

   subroutine floor_and_ceiling_contour
      !! contornos en el piso y en el techo
      use dimensions, only: nx1, nz1, dt1, dx1
      use dinamic_var_perturbation, only: theta_new, theta_base, &
                                          w_perturbed_new, u_perturbed_new, v_perturbed_new
      use cant01, only: dx2
      use initial_z_state, only: theta_z_initial, vapor_z_initial, aerosol_z_initial
      use model_var, only: auxx, auxy, auxz, turbu, aeraux, lapla, cks, Qvap, &
                           qv, qg
      use microphysics_perturbation, only: vapor_base, vapor_new, &
                                           drop_new, crystal_new, rain_new, &
                                           snow_new, hail_new, aerosol_base, aerosol_new
      implicit none
      integer :: i, j

      Qvap = (qv + qg)/nx1**2.*nz1*.1
      do concurrent(i=1:nx1, j=1:nx1)
         theta_new(i, j, 0) = theta_base(i, j, 0) - &
                              w_perturbed_new(i, j, 1)*(theta_z_initial(0) + theta_z_initial(1))*dt1/dx2
         theta_new(i, j, nz1) = theta_new(i, j, nz1 - 1)

         !suponemos que las velocidades horizontales a nivel de piso son
         !iguales a 1/4 de la correspondiente en el nivel 1
         auxx = ((u_perturbed_new(i + 1, j, 1) + u_perturbed_new(i, j, 1))*(vapor_base(i + 1, j, 0) + vapor_base(i, j, 0)) &
                 - (u_perturbed_new(i - 1, j, 1) + u_perturbed_new(i, j, 1))*(vapor_base(i - 1, j, 0) + vapor_base(i, j, 0))) &
                /4.*.25

         auxy = ((v_perturbed_new(i, j + 1, 1) + v_perturbed_new(i, j, 1))*(vapor_base(i, j + 1, 0) + vapor_base(i, j, 0)) &
                 - (v_perturbed_new(i, j - 1, 1) + v_perturbed_new(i, j, 1))*(vapor_base(i, j - 1, 0) + vapor_base(i, j, 0))) &
                /4.*.25

         auxz = w_perturbed_new(i, j, 1) &
                *((vapor_base(i, j, 1) + vapor_base(i, j, 0)) + vapor_z_initial(1) + vapor_z_initial(0))/2.*.5

         vapor_new(i, j, nz1) = vapor_new(i, j, nz1 - 1)

         drop_new(i, j, 0) = drop_new(i, j, 1)
         drop_new(i, j, nz1) = drop_new(i, j, nz1 - 1)

         crystal_new(i, j, 0) = crystal_new(i, j, 1)
         crystal_new(i, j, nz1) = crystal_new(i, j, nz1 - 1)

         !Para que no se acumulen el piso
         rain_new(i, j, 0) = rain_new(i, j, 1)/2.
         rain_new(i, j, nz1) = rain_new(i, j, nz1 - 1)

         snow_new(i, j, 0) = snow_new(i, j, 1)
         snow_new(i, j, nz1) = snow_new(i, j, nz1 - 1)

         !Para que no se acumulen el piso
         hail_new(i, j, 0) = hail_new(i, j, 1)/2.
         hail_new(i, j, nz1) = hail_new(i, j, nz1 - 1)

         !suponemos que las velocidades horizontales a nivel de piso son
         !iguales a 1/4 de la correspondiente en el nivel 1

         auxx = ((u_perturbed_new(i + 1, j, 1) + u_perturbed_new(i, j, 1))*(aerosol_base(i + 1, j, 0) + aerosol_base(i, j, 0)) &
                 - (u_perturbed_new(i - 1, j, 1) + u_perturbed_new(i, j, 1))*(aerosol_base(i - 1, j, 0) + aerosol_base(i, j, 0))) &
                /4.*.25
         auxy = ((v_perturbed_new(i, j + 1, 1) + v_perturbed_new(i, j, 1))*(aerosol_base(i, j + 1, 0) + aerosol_base(i, j, 0)) &
                 - (v_perturbed_new(i, j - 1, 1) + v_perturbed_new(i, j, 1))*(aerosol_base(i, j - 1, 0) + aerosol_base(i, j, 0))) &
                /4.*.25
         auxz = w_perturbed_new(i, j, 1)* &
                ((aerosol_base(i, j, 1) + aerosol_base(i, j, 0)) + aerosol_z_initial(1) + aerosol_z_initial(0))/2.*.5

         if (w_perturbed_new(i, j, 0) > 0) then
            aeraux = -((auxx + auxy) + 2.*auxz)*dt1/dx1
         else
            !se refleja un 25 % de los aerosoles que caen
            aeraux = -((auxx + auxy) + .25*2.*auxz)*dt1/dx1
         end if

         !agregamos un termino de turbulencia para los aerosoles
         !a nivel de piso
         turbu = cks/dx1*.25*(abs(u_perturbed_new(i, j, 1)) &
                              + abs(v_perturbed_new(i, j, 1)) &
                              + 2.*abs(w_perturbed_new(i, j, 1)))

         lapla = ((aerosol_base(i + 1, j, 0) + aerosol_base(i - 1, j, 0)) &
                  + (aerosol_base(i, j + 1, 0) + aerosol_base(i, j - 1, 0) &
                     + aerosol_base(i, j, 1))) &
                 - 5.*aerosol_base(i, j, 0)

         lapla = lapla + (aerosol_z_initial(1) - aerosol_z_initial(0))

         aerosol_new(i, j, 0) = aeraux + aerosol_base(i, j, 0) + turbu*lapla

         aerosol_new(i, j, nz1) = aerosol_new(i, j, nz1 - 1)
      end do
   end subroutine floor_and_ceiling_contour

   subroutine lateral_contour
      !! contornos laterales
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: theta_new
      use microphysics_perturbation, only: vapor_new, drop_new, &
                                           rain_new, crystal_new, snow_new, &
                                           hail_new, aerosol_new
      implicit none
      integer :: j, k
      do concurrent(k=1:nz1 - 1, j=1:nx1)
         theta_new(0, j, k) = theta_new(1, j, k)
         theta_new(nx1 + 1, j, k) = theta_new(nx1, j, k)
         theta_new(j, 0, k) = theta_new(j, 1, k)
         theta_new(j, nx1 + 1, k) = theta_new(j, nx1, k)
         vapor_new(0, j, k) = 0.
         vapor_new(nx1 + 1, j, k) = 0.
         vapor_new(j, 0, k) = 0.
         vapor_new(j, nx1 + 1, k) = 0.
         drop_new(0, j, k) = drop_new(1, j, k)
         drop_new(nx1 + 1, j, k) = drop_new(nx1, j, k)
         drop_new(j, 0, k) = drop_new(j, 1, k)
         drop_new(j, nx1 + 1, k) = drop_new(j, nx1, k)
         rain_new(0, j, k) = 0.
         rain_new(nx1 + 1, j, k) = 0.
         rain_new(j, 0, k) = 0.
         rain_new(j, nx1 + 1, k) = 0.
         crystal_new(0, j, k) = crystal_new(1, j, k)
         crystal_new(nx1 + 1, j, k) = crystal_new(nx1, j, k)
         crystal_new(j, 0, k) = crystal_new(j, 1, k)
         crystal_new(j, nx1 + 1, k) = crystal_new(j, nx1, k)
         snow_new(0, j, k) = snow_new(1, j, k)
         snow_new(nx1 + 1, j, k) = snow_new(nx1, j, k)
         snow_new(j, 0, k) = snow_new(j, 1, k)
         snow_new(j, nx1 + 1, k) = snow_new(j, nx1, k)
         hail_new(0, j, k) = 0.
         hail_new(nx1 + 1, j, k) = 0.
         hail_new(j, 0, k) = 0.
         hail_new(j, nx1 + 1, k) = 0.
         aerosol_new(0, j, k) = aerosol_new(1, j, k)
         aerosol_new(nx1 + 1, j, k) = aerosol_new(nx1, j, k)
         aerosol_new(j, 0, k) = aerosol_new(j, 1, k)
         aerosol_new(j, nx1 + 1, k) = aerosol_new(j, nx1, k)
      end do
   end subroutine lateral_contour

   subroutine speed_pressure()
      !! calculo de la velocidad y la presion
      !! Calcula la evolucion del la presion y las velocidades con un paso de tiempo menor lt3
      !! Las cantidades 1 son las presentes en el paso grande y las 2 son las del paso futuro, las 3 son auxiliares
      !! Le resta la perturbacion promedio
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use constants
      use initial_z_state
      use velpre01
      use p3v3
      use sv_inhomogeneous_velocities_and_speed_pressure, only: fu, fv, fw, fp
      implicit none

      call velpre01_init()

      do concurrent(i=0:nx1 + 1, j=0:nx1 + 1, k=0:nz1)
         u_perturbed_new(i, j, k) = u_perturbed_base(i, j, k)
         v_perturbed_new(i, j, k) = v_perturbed_base(i, j, k)
         w_perturbed_new(i, j, k) = w_perturbed_base(i, j, k)
         pressure_new(i, j, k) = pressure_base(i, j, k)
      end do

      do concurrent(t=1:lt3)
         do concurrent(k=1:nz1 - 1)
            presi = -Cp*theta_z_initial(k)*(1.+.61*vapor_z_initial(k)/air_density_z_initial(k))
            vel0 = theta_z_initial(k)*(air_density_z_initial(k) + .61*vapor_z_initial(k))
            vel1 = theta_z_initial(k - 1)*(air_density_z_initial(k - 1) + .61*vapor_z_initial(k - 1))
            vel2 = theta_z_initial(k + 1)*(air_density_z_initial(k + 1) + .61*vapor_z_initial(k + 1))
            vel3 = cc2(k)/presi/vel0
            do concurrent(i=1:nx1, j=1:nx1)
               dprex = pressure_new(i + 1, j, k) - pressure_new(i - 1, j, k)
               dprey = pressure_new(i, j + 1, k) - pressure_new(i, j - 1, k)
               dprez = pressure_new(i, j, k + 1) - pressure_new(i, j, k - 1)

               presix = presi*dprex/dx2
               presiy = presi*dprey/dx2
               presiz = presi*dprez/dx2

               U3(i, j, k) = dt3*(presix + fu(i, j, k)) + u_perturbed_new(i, j, k)
               V3(i, j, k) = dt3*(presiy + fv(i, j, k)) + v_perturbed_new(i, j, k)
               W3(i, j, k) = dt3*(presiz + fw(i, j, k)) + w_perturbed_new(i, j, k)

               dvx = vel0*(u_perturbed_new(i + 1, j, k) - u_perturbed_new(i - 1, j, k))
               dvy = vel0*(v_perturbed_new(i, j + 1, k) - v_perturbed_new(i, j - 1, k))
               if (k == 1) then
                  !      dvz = tiene 80% de (w_perturbed_new(2)-w_perturbed_new(1) y 20% de (w_perturbed_new(1)-w_perturbed_new(0)
                  dvz = (.8*vel2*w_perturbed_new(i, j, k + 1) - .8*vel1*w_perturbed_new(i, j, k))*2.
               else
                  dvz = vel2*w_perturbed_new(i, j, k + 1) - vel1*w_perturbed_new(i, j, k - 1)
               end if

               diver = vel3*((dvx + dvy) + dvz)/dx2

               !      modificado para agrega turbulencia en la P 23/8/97
               Pres3(i, j, k) = dt3*(diver + fp(i, j, k)) + pressure_new(i, j, k)
            end do
         end do

         !*      redefiniciones y contornos
         do concurrent(i=1:nx1, j=1:nx1)
            Pres3(i, j, 0) = Pres3(i, j, 1)
            Pres3(i, j, nz1) = Pres3(i, j, nz1 - 1)
         end do
         do concurrent(i=1:nx1, k=0:nz1)
            Pres3(i, 0, k) = Pres3(i, 1, k)
            Pres3(i, nx1 + 1, k) = Pres3(i, nx1, k)
            Pres3(0, i, k) = Pres3(1, i, k)
            Pres3(nx1 + 1, i, k) = Pres3(nx1, i, k)
         end do

         presprom = 0.
         do concurrent(i=1:nx1, j=1:nx1)
            do k = 1, nz1 - 1
               if (k == 1) then
                  u_perturbed_new(i, j, k) = U3(i, j, k) - kkk* &
                                             (2.*U3(i, j, k) - U3(i, j, k + 1))
                  v_perturbed_new(i, j, k) = V3(i, j, k) - kkk* &
                                             (2.*V3(i, j, k) - V3(i, j, k + 1))
                  w_perturbed_new(i, j, k) = W3(i, j, k) - kkk* &
                                             (2.*W3(i, j, k) - W3(i, j, k + 1))
               else
                  u_perturbed_new(i, j, k) = U3(i, j, k)
                  v_perturbed_new(i, j, k) = V3(i, j, k)
                  w_perturbed_new(i, j, k) = W3(i, j, k)
               end if
               pressure_new(i, j, k) = prom1*Pres3(i, j, k) + prom*( &
                                       ((Pres3(i + 1, j, k) + Pres3(i - 1, j, k)) + &
                                        (Pres3(i, j + 1, k) + Pres3(i, j - 1, k))) + &
                                       Pres3(i, j, k + 1) + Pres3(i, j, k - 1))
               presprom = pressure_new(i, j, k) + presprom
            end do

            u_perturbed_new(i, j, 0) = 0
            v_perturbed_new(i, j, 0) = 0
            w_perturbed_new(i, j, 0) = 0
            pressure_new(i, j, 0) = pressure_new(i, j, 1)
            u_perturbed_new(i, j, nz1) = u_perturbed_new(i, j, nz1 - 1)
            v_perturbed_new(i, j, nz1) = v_perturbed_new(i, j, nz1 - 1)
            w_perturbed_new(i, j, nz1) = w_perturbed_new(i, j, nz1 - 1)
            pressure_new(i, j, nz1) = pressure_new(i, j, nz1 - 1)
         end do
         do concurrent(i=1:nx1, k=0:nz1)
            u_perturbed_new(0, i, k) = u_perturbed_new(1, i, k)
            v_perturbed_new(0, i, k) = v_perturbed_new(1, i, k)
            w_perturbed_new(0, i, k) = w_perturbed_new(1, i, k)
            pressure_new(0, i, k) = pressure_new(1, i, k)
            u_perturbed_new(nx1 + 1, i, k) = u_perturbed_new(nx1, i, k)
            v_perturbed_new(nx1 + 1, i, k) = v_perturbed_new(nx1, i, k)
            w_perturbed_new(nx1 + 1, i, k) = w_perturbed_new(nx1, i, k)
            pressure_new(nx1 + 1, i, k) = pressure_new(nx1, i, k)
            u_perturbed_new(i, 0, k) = u_perturbed_new(i, 1, k)
            v_perturbed_new(i, 0, k) = v_perturbed_new(i, 1, k)
            w_perturbed_new(i, 0, k) = w_perturbed_new(i, 1, k)
            pressure_new(i, 0, k) = pressure_new(i, 1, k)
            u_perturbed_new(i, nx1 + 1, k) = u_perturbed_new(i, nx1, k)
            v_perturbed_new(i, nx1 + 1, k) = v_perturbed_new(i, nx1, k)
            w_perturbed_new(i, nx1 + 1, k) = w_perturbed_new(i, nx1, k)
            pressure_new(i, nx1 + 1, k) = pressure_new(i, nx1, k)
         end do

         presprom = presprom/nnn
         do concurrent(i=0:nx1 + 1, j=0:nx1 + 1, k=0:nz1)
            pressure_new(i, j, k) = pressure_new(i, j, k) - presprom
         end do

         if (t == lt3/2) then
            do concurrent(i=0:nx1 + 1, j=0:nx1 + 1, k=0:nz1)
               u_perturbed_base(i, j, k) = u_perturbed_new(i, j, k)
               v_perturbed_base(i, j, k) = v_perturbed_new(i, j, k)
               w_perturbed_base(i, j, k) = w_perturbed_new(i, j, k)
               pressure_base(i, j, k) = pressure_new(i, j, k)
            end do
         end if

      end do

      !**********************************************************
      !*    suavizado

      call filtro(pressure_base, .15, .15, .1)

      call filtro(pressure_new, .15, .15, .1)

      call filtro(u_perturbed_base, facx, facy, facz)
      call filtro(u_perturbed_new, facx, facy, facz)
      call filtro(v_perturbed_base, facx, facy, facz)
      call filtro(v_perturbed_new, facx, facy, facz)
      call filtro(w_perturbed_base, facx, facy, facz)
      call filtro(w_perturbed_new, facx, facy, facz)

      do concurrent(i=1:nx1, j=1:nx1)
         pressure_base(i, j, 0) = pressure_base(i, j, 1)
         pressure_new(i, j, 0) = pressure_new(i, j, 1)
      end do
      !**********************************************************

      return
   end subroutine speed_pressure

   subroutine filtro(varia1, facx, facy, facz)
      !! filtro para theta_base vapor_base
      !! Esta subrutina filtra componentes de alta frecuencia espacial.
      !! El valor de la variable del punto j se filtra con los valores
      !! extrapolados linalmente de los puntos j-3 y j-1 y similares,
      !! pasando un polinomio de grado 4
      use dimensions
      use filtro01
      implicit none
      character*50 text
      REAL, DIMENSION(-3:NX1 + 3, -3:NX1 + 3, -2:NZ1 + 2), intent(inout) :: varia1
      real, intent(in) :: facx, facy, facz
      fact = 1.-(facx + facy + facz)
      !**********************************************************
      !     Redefiniciones y contornos

      do concurrent(i=0:nx1 + 1, j=0:nx1 + 1, k=0:nz1)
         varia2(i, j, k) = varia1(i, j, k)
      end do

      do concurrent(k=0:nz1, i=0:nx1)
         varia2(i, -1, k) = varia2(i, 1, k)
         varia2(i, -2, k) = varia2(i, 1, k)
         varia2(i, nx1 + 2, k) = varia2(i, nx1, k)
         varia2(i, nx1 + 3, k) = varia2(i, nx1, k)
         varia2(-1, i, k) = varia2(1, i, k)
         varia2(-2, i, k) = varia2(1, i, k)
         varia2(nx1 + 2, i, k) = varia2(nx1, i, k)
         varia2(nx1 + 3, i, k) = varia2(nx1, i, k)
      end do

      do concurrent(i=1:nx1, j=1:nx1)
         varia2(i, j, -1) = varia2(i, j, 0)
         varia2(i, j, -2) = varia2(i, j, 0)
         varia2(i, j, nz1 + 1) = varia2(i, j, nz1)
         varia2(i, j, nz1 + 2) = varia2(i, j, nz1)
      end do

      !**********************************************************
      !     Filtro

      do concurrent(i=1:nx1, j=1:nx1, k=1:nz1 - 1)
         varx = (9.*(varia2(i - 1, j, k) + varia2(i + 1, j, k)) - &
                 (varia2(i - 3, j, k) + varia2(i + 3, j, k)))/16.
         vary = (9.*(varia2(i, j - 1, k) + varia2(i, j + 1, k)) - &
                 (varia2(i, j - 3, k) + varia2(i, j + 3, k)))/16.
         varz = (9.*(varia2(i, j, k - 1) + varia2(i, j, k + 1)) - &
                 (varia2(i, j, k - 3) + varia2(i, j, k + 3)))/16.

         varia1(i, j, k) = ((facx*varx + facy*vary) + facz*varz) + &
                           fact*varia2(i, j, k)
      end do

      do concurrent(k=1:nz1 - 1, i=1:nx1)
         varia1(i, 0, k) = varia1(i, 1, k)
         varia1(i, nx1 + 1, k) = varia1(i, nx1, k)
         varia1(0, i, k) = varia1(1, i, k)
         varia1(nx1 + 1, i, k) = varia1(nx1, i, k)
      end do

      do concurrent(i=1:nx1, j=1:nx1)
         varia1(i, j, nz1) = varia1(i, j, nz1 - 1)
      end do
      !**********************************************************
      return
   end subroutine filtro

   subroutine floor_condition_redefinition()
      !! modificada las condiciones en el piso
      !! Redefinicion
      use dimensions, only: nx1, nz1, dt1, dx1
      use dinamic_var_perturbation, only: theta_new, theta_base, &
                                          w_perturbed_new
      use cant01, only: pro3, pro4, pro1, pro2
      use microphysics_perturbation, only: vapor_base, vapor_new, &
                                           drop_new, rain_new, crystal_new, &
                                           snow_new, hail_new, aerosol_base, aerosol_new, &
                                           drop_base, rain_base, crystal_base, snow_base, hail_base
      use initial_z_state, only: aerosol_z_initial
      use model_var, only: aeraux
      implicit none
      integer :: i, j, k
      do concurrent(i=1:nx1, j=1:nx1)
         k = 0
         theta_base(i, j, k) = pro3*theta_new(i, j, k) + &
                               pro4*( &
                               (theta_new(i + 1, j, k) + theta_new(i - 1, j, k)) + &
                               (theta_new(i, j + 1, k) + theta_new(i, j - 1, k)))

         if (abs(theta_base(i, j, k)) < 1e-10) theta_base(i, j, k) = 0

         vapor_base(i, j, k) = pro3*vapor_new(i, j, k) &
                               + pro4*( &
                               (vapor_new(i + 1, j, k) + vapor_new(i - 1, j, k)) + &
                               (vapor_new(i, j + 1, k) + vapor_new(i, j - 1, k)))

         if (abs(vapor_base(i, j, k)) < 1e-10) vapor_base(i, j, k) = 0

         drop_base(i, j, k) = pro3*drop_new(i, j, k) + &
                              pro4*( &
                              (drop_new(i + 1, j, k) + drop_new(i - 1, j, k)) + &
                              (drop_new(i, j + 1, k) + drop_new(i, j - 1, k)))

         if (drop_base(i, j, k) < 1e-10) drop_base(i, j, k) = 0

         rain_base(i, j, k) = rain_new(i, j, k)

         if (rain_base(i, j, k) < 1e-10) rain_base(i, j, k) = 0

         crystal_base(i, j, k) = pro3*crystal_new(i, j, k) + &
                                 pro4*( &
                                 (crystal_new(i + 1, j, k) + crystal_new(i - 1, j, k)) + &
                                 (crystal_new(i, j + 1, k) + crystal_new(i, j - 1, k)))

         if (crystal_base(i, j, k) < 1e-10) crystal_base(i, j, k) = 0

         snow_base(i, j, k) = pro3*snow_new(i, j, k) + &
                              pro4*( &
                              (snow_new(i + 1, j, k) + snow_new(i - 1, j, k)) + &
                              (snow_new(i, j + 1, k) + snow_new(i, j - 1, k)))

         if (snow_base(i, j, k) < 1e-10) snow_base(i, j, k) = 0

         hail_base(i, j, k) = hail_new(i, j, k)

         if (hail_base(i, j, k) < 1e-10) hail_base(i, j, k) = 0

         aerosol_base(i, j, k) = pro3*aerosol_new(i, j, k) + &
                                 pro4*( &
                                 (aerosol_new(i + 1, j, k) + aerosol_new(i - 1, j, k)) + &
                                 (aerosol_new(i, j + 1, k) + aerosol_new(i, j - 1, k)))

         !correccion cambiando la absorcion de aerosoles
         if ((rain_base(i, j, 0) + hail_base(i, j, 1)) > 1e-6 .and. w_perturbed_new(i, j, 1) < -0.1) then
            aeraux = (-w_perturbed_new(i, j, 1)/1.)*(rain_base(i,j,0)/1e-3)*.2*(dt1/5.)
            aerosol_base(i, j, k) = aerosol_base(i, j, k) - (aerosol_base(i, j, k) + aerosol_z_initial(k))*aeraux
         end if

         if (abs(aerosol_base(i, j, k)) < 1e-10) aerosol_base(i, j, k) = 0

         do concurrent(k=1:nz1 - 1)
            theta_base(i, j, k) = pro1*theta_new(i, j, k) + &
                                  pro2*( &
                                  (theta_new(i + 1, j, k) + theta_new(i - 1, j, k)) + &
                                  (theta_new(i, j + 1, k) + theta_new(i, j - 1, k)) + &
                                  theta_new(i, j, k + 1) + &
                                  theta_new(i, j, k - 1))

            if (abs(theta_base(i, j, k)) < 1e-10) theta_base(i, j, k) = 0

            vapor_base(i, j, k) = pro1*vapor_new(i, j, k) + &
                                  pro2*( &
                                  (vapor_new(i + 1, j, k) + vapor_new(i - 1, j, k)) + &
                                  (vapor_new(i, j + 1, k) + vapor_new(i, j - 1, k)) + &
                                  vapor_new(i, j, k + 1) + &
                                  vapor_new(i, j, k - 1))

            if (abs(vapor_base(i, j, k)) < 1e-10) vapor_base(i, j, k) = 0

            drop_base(i, j, k) = pro1*drop_new(i, j, k) + &
                                 pro2*( &
                                 (drop_new(i + 1, j, k) + drop_new(i - 1, j, k)) + &
                                 (drop_new(i, j + 1, k) + drop_new(i, j - 1, k)) + &
                                 drop_new(i, j, k + 1) + &
                                 drop_new(i, j, k - 1))

            if (drop_base(i, j, k) < 1e-10) drop_base(i, j, k) = 0

            rain_base(i, j, k) = pro1*rain_new(i, j, k) + &
                                 pro2*( &
                                 (rain_new(i + 1, j, k) + rain_new(i - 1, j, k)) + &
                                 (rain_new(i, j + 1, k) + rain_new(i, j - 1, k)) + &
                                 rain_new(i, j, k + 1) + &
                                 rain_new(i, j, k - 1))

            if (rain_base(i, j, k) < 1e-10) rain_base(i, j, k) = 0

            crystal_base(i, j, k) = pro1*crystal_new(i, j, k) + &
                                    pro2*( &
                                    (crystal_new(i + 1, j, k) + crystal_new(i - 1, j, k)) + &
                                    (crystal_new(i, j + 1, k) + crystal_new(i, j - 1, k)) + &
                                    crystal_new(i, j, k + 1) + &
                                    crystal_new(i, j, k - 1))

            if (crystal_base(i, j, k) < 1e-10) crystal_base(i, j, k) = 0

            snow_base(i, j, k) = pro1*snow_new(i, j, k) + &
                                 pro2*( &
                                 (snow_new(i + 1, j, k) + snow_new(i - 1, j, k)) + &
                                 (snow_new(i, j + 1, k) + snow_new(i, j - 1, k)) + &
                                 snow_new(i, j, k + 1) + &
                                 snow_new(i, j, k - 1))

            if (snow_base(i, j, k) < 1e-10) snow_base(i, j, k) = 0

            hail_base(i, j, k) = pro1*hail_new(i, j, k) + &
                                 pro2*( &
                                 (hail_new(i + 1, j, k) + hail_new(i - 1, j, k)) + &
                                 (hail_new(i, j + 1, k) + hail_new(i, j - 1, k)) + &
                                 hail_new(i, j, k + 1) + &
                                 hail_new(i, j, k - 1))

            if (hail_base(i, j, k) < 1e-10) hail_base(i, j, k) = 0

            aerosol_base(i, j, k) = pro1*aerosol_new(i, j, k) + &
                                    pro2*( &
                                    (aerosol_new(i + 1, j, k) + aerosol_new(i - 1, j, k)) + &
                                    (aerosol_new(i, j + 1, k) + aerosol_new(i, j - 1, k)) + &
                                    aerosol_new(i, j, k + 1) + &
                                    aerosol_new(i, j, k - 1))

            if (abs(aerosol_base(i, j, k)) < 1e-10) aerosol_base(i, j, k) = 0
         end do
      end do
   end subroutine floor_condition_redefinition

   subroutine floor_and_ceiling_contour_redefinition()
      !! contornos en el piso y en el techo
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: theta_base
      use microphysics_perturbation, only: vapor_base, aerosol_base, drop_base, &
                                           rain_base, crystal_base, snow_base, hail_base
      use initial_z_state, only: vapor_z_initial, aerosol_z_initial
      implicit none
      integer :: i, j
      do concurrent(i=1:nx1, j=1:nx1)
         theta_base(i, j, 0) = theta_base(i, j, 0)
         if (theta_base(i, j, 0) > 0.5) theta_base(i, j, 0) = .5
         if (-theta_base(i, j, 0) > 0.5) theta_base(i, j, 0) = -.5

         theta_base(i, j, nz1) = theta_base(i, j, nz1 - 1)

         !corregido para el vapor
         if (vapor_base(i, j, 0) > vapor_z_initial(0)*.5) then
            vapor_base(i, j, 0) = .8*vapor_z_initial(0)
         end if
         if (-vapor_base(i, j, 0) > vapor_z_initial(0)*.5) then
            vapor_base(i, j, 0) = -.8*vapor_z_initial(0)
         end if

         vapor_base(i, j, nz1) = vapor_base(i, j, nz1 - 1)
         drop_base(i, j, 0) = 0.
         drop_base(i, j, nz1) = drop_base(i, j, nz1 - 1)
         rain_base(i, j, 0) = rain_base(i, j, 0)
         rain_base(i, j, nz1) = rain_base(i, j, nz1 - 1)
         crystal_base(i, j, 0) = 0.
         crystal_base(i, j, nz1) = crystal_base(i, j, nz1 - 1)

         snow_base(i, j, 0) = 0.
         snow_base(i, j, nz1) = snow_base(i, j, nz1 - 1)

         hail_base(i, j, 0) = hail_base(i, j, 0)
         hail_base(i, j, nz1) = hail_base(i, j, nz1 - 1)

         !corregido para los aerosoles
         if (-aerosol_base(i, j, 0) > 0.8*aerosol_z_initial(0)) then
            aerosol_base(i, j, 0) = -.8*aerosol_z_initial(0)
         end if
         aerosol_base(i, j, nz1) = aerosol_base(i, j, nz1 - 1)
      end do
   end subroutine floor_and_ceiling_contour_redefinition

   subroutine lateral_contour_redefinition()
      !! contornos laterales
      use dimensions, only: nx1, nz1
      use dinamic_var_perturbation, only: theta_base
      use microphysics_perturbation, only: vapor_base, drop_base, rain_base, &
                                           crystal_base, snow_base, hail_base, aerosol_base
      implicit none
      integer :: j, k
      do concurrent(k=1:nz1 - 1, j=1:nx1)
         theta_base(0, j, k) = theta_base(1, j, k)
         theta_base(nx1 + 1, j, k) = theta_base(nx1, j, k)
         theta_base(j, 0, k) = theta_base(j, 1, k)
         theta_base(j, nx1 + 1, k) = theta_base(j, nx1, k)
         vapor_base(0, j, k) = 0.
         vapor_base(nx1 + 1, j, k) = 0.
         vapor_base(j, 0, k) = 0.
         vapor_base(j, nx1 + 1, k) = 0.
         drop_base(0, j, k) = drop_base(1, j, k)
         drop_base(nx1 + 1, j, k) = drop_base(nx1, j, k)
         drop_base(j, 0, k) = drop_base(j, 1, k)
         drop_base(j, nx1 + 1, k) = drop_base(j, nx1, k)
         rain_base(0, j, k) = 0.
         rain_base(nx1 + 1, j, k) = 0.
         rain_base(j, 0, k) = 0.
         rain_base(j, nx1 + 1, k) = 0.
         crystal_base(0, j, k) = crystal_base(1, j, k)
         crystal_base(nx1 + 1, j, k) = crystal_base(nx1, j, k)
         crystal_base(j, 0, k) = crystal_base(j, 1, k)
         crystal_base(j, nx1 + 1, k) = crystal_base(j, nx1, k)
         snow_base(0, j, k) = snow_base(1, j, k)
         snow_base(nx1 + 1, j, k) = snow_base(nx1, j, k)
         snow_base(j, 0, k) = snow_base(j, 1, k)
         snow_base(j, nx1 + 1, k) = snow_base(j, nx1, k)
         hail_base(0, j, k) = 0.
         hail_base(nx1 + 1, j, k) = 0.
         hail_base(j, 0, k) = 0.
         hail_base(j, nx1 + 1, k) = 0.

         aerosol_base(0, j, k) = aerosol_base(1, j, k)
         aerosol_base(nx1 + 1, j, k) = aerosol_base(nx1, j, k)
         aerosol_base(j, 0, k) = aerosol_base(j, 1, k)
         aerosol_base(j, nx1 + 1, k) = aerosol_base(j, nx1, k)
      end do
   end subroutine lateral_contour_redefinition

   subroutine vapour_negative_correction()
      !! correccion de negativos para el vapor
      use dimensions, only: nx1, nz1
      use microphysics_perturbation, only: vapor_base
      use initial_z_state, only: vapor_z_initial
      implicit none
      integer :: i, j, k
      do concurrent(i=0:nx1 + 1, j=0:nx1 + 1, k=0:nz1)
         if (vapor_base(i, j, k) + vapor_z_initial(k) < 0) then
            vapor_base(i, j, k) = -vapor_z_initial(k)
         end if
      end do
   end subroutine vapour_negative_correction

   subroutine save_backup()
      use model_var, only: current_time, tte, posx, posy, Xnub, Ynub, posxx, posyy, &
                           file_number, actual_file
      use cant01, only: lte, ltg, ltb
      use dimensions, only: dt1
      use config, only: output_directory
      use dinamic_var_perturbation, only: w_perturbed_new, u_perturbed_new, &
                                          v_perturbed_new, heat_force, theta_new, theta_base, pressure_new, &
                                          pressure_base, u_perturbed_base, v_perturbed_base, w_perturbed_base
      use io, only: str_gen
      use microphysics_perturbation, only: aerosol_base, drop_new, &
                                           rain_new, crystal_new, snow_new, &
                                           hail_new, vapor_new, aerosol_new, &
                                           vapor_base, drop_base, rain_base, crystal_base, snow_base, hail_base, &
                                           Av, Vtgra0, Vtnie
      use constants, only: Tvis, Telvs, Tesvs, Tlvl, Tlsl, Tlvs, Eacrcn, Eautcn
      use initial_z_state, only: air_density_z_initial, vapor_z_initial, &
                                 aerosol_z_initial, theta_z_initial, temperature_z_initial, &
                                 Pres00, aerosol_z_relative, cc2, &
                                 vapor_z_relative, u_z_initial, v_z_initial
      use model_initialization, only: cloud_position, cloud_movement, &
                                      statistics
      use graba, only: graba320, graba120
      implicit none
      integer :: unit_number
      if (current_time/nint(lte/dt1)*nint(lte/dt1) == current_time) then
         call statistics()
         tte = tte + 1
         call cloud_position()
         call cloud_movement()

         open (newunit=unit_number, file=output_directory//"posnub"//'.sa', &
               ACCESS="append")
         write (unit_number, *) tte, posx(tte), posy(tte), Xnub(tte), Ynub(tte), posxx, posyy
         close (unit_number)
      end if

      if (current_time/nint(ltg/dt1)*nint(ltg/dt1) == current_time) then
         file_number = str_gen(actual_file)
         call graba320(u_perturbed_base, v_perturbed_base, w_perturbed_base, &
                       theta_base, pressure_base, vapor_base, drop_base, &
                       rain_base, crystal_base, snow_base, hail_base, aerosol_base, file_number)
         actual_file = actual_file + 1
      end if

      if (current_time/nint(ltb/dt1)*nint(ltb/dt1) == current_time) then
         call graba120(air_density_z_initial, temperature_z_initial, theta_z_initial, &
                       Pres00, vapor_z_initial, cc2, &
                       aerosol_z_initial, u_z_initial, v_z_initial, &
                       u_perturbed_base, u_perturbed_new, v_perturbed_base, v_perturbed_new, &
                       w_perturbed_base, w_perturbed_new, theta_base, &
                       theta_new, pressure_base, pressure_new, vapor_base, &
                       vapor_new, drop_base, drop_new, rain_base, rain_new, crystal_base, &
                       crystal_new, snow_base, snow_new, hail_base, hail_new, aerosol_base, &
                       aerosol_new, heat_force, Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Av, &
                       Vtnie, Vtgra0, vapor_z_relative, aerosol_z_relative, Eautcn, Eacrcn)
      end if
   end subroutine save_backup
end module model_aux

