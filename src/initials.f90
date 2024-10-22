module velpre01
   integer :: t, i, j, k
   real :: dvx, dvy, dvz, diver, dprex, dprey, dprez, vel0, vel1, vel2, &
      vel3, presi, presix, presiy, presiz, facx, facy, facz, prom1, prom, &
      kkk, presprom, nnn
contains
   subroutine velpre01_init()
      use dimensions
      implicit none
      facx = .05
      facy = .05
      facz = .05

      prom = .3/6.*(dt3/.2)
      prom1 = 1.-prom*6.
      kkk = .01
      nnn = (nx1 + 2)**2.*(nz1 + 1)

   end subroutine velpre01_init
end module velpre01

module model_initial_conditions
contains

   subroutine initial_conditions()
      !!     Condiciones iniciales para las variables dinamicas
      !!     Corresponde a una nube con hielo
      !!     Calcula la presion inicial en forma iterativa, primero supone aire
      !!      seco y luego integra considerando la densidad total, incluyendo el
      !!      vapor.
      !!     Solo hay una perturbacion en temperatura para iniciar la conveccion
      !!     Incluye viento de corte, tipo frente
      use cant01
      use dimensions
      use microphysics_perturbation
      use dinamic_var_perturbation
      use constants
      use initial_z_state
      use config

      implicit none

      real x_aux, y_aux, z_aux
      real :: cenaerx=(nx1+1.)*300/2.+3000.   !Coord x de la perturbacion de aerosoles
      real :: cenaery=(nx1+1.)*300/2.   !Coord y de la perturbacion de aerosoles

      real vapor_total, aerosol_total, gaussian
      real :: initial_x_perturbation  !! Initial disturbance’s x-coordinate
      real :: initial_y_perturbation  !! Initial disturbance’s y-coordinate
      real :: initial_z_perturbation = 0. !! Initial disturbance’s z-coordinate
      real :: sigma_t = 2*1000.**2.   !! z decay of the perturbation in T
      real :: sigma_a = 200.**2. !! z decay of the perturbation in A
      real :: perturbation_width = 2000. !! perturbation_width
      real :: temperature_max_perturbation = .7  !! Maximum temperature perturbation
      real :: aerosol_max_perturbation = 10000. !! Maximum aerosol perturbation

      integer i, j, k, unit
      initial_x_perturbation = (nx1 + 1.)*dx1/2.  !! Initial disturbance’s x-coordinate
      initial_y_perturbation = (nx1 + 1.)*dx1/2.
      
      do concurrent(k=0:nz1)
         z_aux = k*dx1
         if (z_aux <= 500.) then
            u_z_initial(k) = 0.
            v_z_initial(k) = 0.

         elseif (z_aux <= 2000.) then
            u_z_initial(k) = 4.*((z_aux - 500.)/1500.)**2.
            v_z_initial(k) = 0.

         elseif (z_aux <= 9000.) then
            u_z_initial(k) = 4.-10.*((z_aux - 2000.)/7000.)**2.
            v_z_initial(k) = 3.*((z_aux - 2000.)/7000.)**.5
            !
         else
            u_z_initial(k) = 4.*((z_aux - 9000.)/9000.)**2.-6.
            v_z_initial(k) = 3.-5.*((z_aux - 9000.)/9000.)**.5
            !
         end if

         u_z_initial(k) = u_z_initial(k)*0.7
         v_z_initial(k) = v_z_initial(k)*0.7
      end do

      open (newunit=unit, file=output_directory//"initial_z_arrays", access='append')
      do concurrent(k=210:313)
         write (unit, *) k, Tvis(k), Tlvl(k), Tlsl(k), Tlvs(k), Telvs(k), Tesvs(k), &
            Eautcn(k), Eacrcn(k)
      end do
      close (unit)

      do concurrent(k=0:nz1, i=1:nx1, j=1:nx1)
         !perturbaciones iniciales en la temperatura y en los aerosoles
         z_aux = k*dx1
         x_aux = i*dx1
         y_aux = j*dx1

         gaussian = exp(-((initial_x_perturbation - x_aux)**2. &
                           +(initial_y_perturbation - y_aux)**2.)*.5/perturbation_width**2.)

         theta_base(i, j, k) = temperature_max_perturbation*exp(-(z_aux &
                                                                  - initial_z_perturbation)**2./sigma_t)*gaussian

         if (theta_base(i, j, k) < 1e-5) theta_base(i, j, k) = 0.

         gaussian = exp(-((cenaerx - x_aux)**2. &
                           +(cenaery - y_aux)**2.)*.5/perturbation_width**2.)

         aerosol_base(i, j, k) = aerosol_max_perturbation*exp( &
                                 -z_aux**2./sigma_a)*gaussian
      end do

      open (newunit=unit, file=output_directory//"inic03.sa", access='append')
      do concurrent(k=0:nz1)
         write (unit, 210) k, temperature_z_initial(k), theta_z_initial(k), &
            Presi0(k), Pres00(k), air_density_z_initial(k), aerosol_z_initial(k), &
            vapor_z_initial(k), u_z_initial(k), v_z_initial(k)
      end do
      close (unit)

      do concurrent(i=1:nx1, j=1:nx1)
         pressure_base(i, j, 0) = pressure_base(i, j, 1)
         pressure_base(i, j, -1) = pressure_base(i, j, 1)
         pressure_new(i, j, 0) = pressure_base(i, j, 1)
         pressure_new(i, j, -1) = pressure_base(i, j, 1)
         theta_base(i, j, 0) = theta_base(i, j, 1)
         theta_base(i, j, -1) = theta_base(i, j, 1)
         vapor_base(i, j, 0) = vapor_base(i, j, 1)
         vapor_base(i, j, -1) = vapor_base(i, j, 1)
      end do

      ! calculo del vapor_z_relative
      vapor_total = 0.
      do concurrent(k=1:nz1)
         vapor_total = vapor_total + vapor_z_initial(k)
      end do

      do concurrent(k=1:nz1)
         vapor_z_relative(k) = vapor_z_initial(k)/vapor_total
      end do

      !     calculo del aerosol_z_relative
      aerosol_total = 0.
      do concurrent(k=1:nz1)
         aerosol_total = aerosol_total + aerosol_z_initial(k)
      end do
      do concurrent(k=1:nz1)
         aerosol_z_relative(k) = aerosol_z_initial(k)/aerosol_total
      end do

210   format(I3, 9E12.4)
      return

   end subroutine initial_conditions
end module model_initial_conditions

module model_var
   real :: T, P, Dv, Lvl, Lvs, Lsl, Vis, Qvap, Qliq, densi, nu, Lsl00, Eaucn, &
      Eaccn, Eacng, Naer, dqgot, dqcri, daer, daer2, Fcal, elvs, esvs, e1, rl, &
      rs, dden0z, aux, aux1, aux2, aux3, aux4, posxx, posyy, cks, turbu, lapla

   real(8) :: qgotaux, qvapaux, qlluaux, qcriaux, qnieaux, qgraaux, aeraux, &
      auxx, auxy, auxz, Taux, Qvapneg, aerneg, ener, ener1, ener2, ener3, &
      ener4, ener5, qv, qg, daitot, vapt1, vapt2, vapt3, vapt4, gott1, gott2, &
      gott3, gott4, aert1, aert2, aert3, aert4, totnuc, totmic

   real(8) :: Xnub(5000), Ynub(5000)
   integer :: posx(-3:5000), posy(-3:5000)

   character(len=3) :: file_number

   integer :: current_time
   integer :: actual_file  !! actual file number
   integer :: t2, n, m, l, i, j, k, lll, s, iT, tte, lvapneg, &
      llluneg, lcrineg, laerneg, lnieneg, lgraneg, yy

end module model_var

module model_initialization
   implicit none

   real :: zmed
   real(8) :: impx, impy, Qagua, Qaguat
   integer :: spos, laux1, laux2, maux1, maux2, naux2, umax, umin, &
      titamin, qvapmax, qvapmin, qgotmax, qllumax, qcrimax, qniemax, &
      qgramax, aermax, lumax, mumax, numax, lumin, mumin, numin, lvmax, &
      mvmax, nvmax, lvmin, mvmin, nvmin, lwmax, mwmax, nwmax, lwmin, &
      mwmin, nwmin, ltitamax, mtitamax, ntitamax, ltitamin, mtitamin, &
      ntitamin, lqvapmax, mqvapmax, nqvapmax, lqvapmin, mqvapmin, nqvapmin, &
      lqgotmax, mqgotmax, nqgotmax, lqllumax, mqllumax, nqllumax, laermax, &
      maermax, naermax, lqcrimax, mqcrimax, nqcrimax, lqniemax, mqniemax, &
      nqniemax, lqgramax, mqgramax, nqgramax, qgottot, qllutot, qcritot, &
      vmax, vmin, wmax, wmin, titamax, qnietot, qgratot

contains

   subroutine initialize_model()
      use model_var
      use cant01, only: total_time, lt2, lt3, cteturb, dx2, dx8, dx12, &
         AA, ikapa, cteqgot, cteqllu, cteqnie, cteqgra, ltt, ltg, lte, &
         ltb, ctur, pro1, pro2, pro3, pro4, cteqnie
      use dimensions
      use constants
      use config, only: sim_time_minutes, save_lapse_minutes, statistic_time_minutes, &
         output_directory, backup_time_minutes, restore_backup
      use initial_z_state, only: air_density_z_initial, temperature_z_initial, &
         theta_z_initial, Pres00, vapor_z_initial, cc2, &
         aerosol_z_initial, u_z_initial, v_z_initial, &
         vapor_z_relative, aerosol_z_relative
      use dinamic_var_perturbation, only: w_perturbed_new, u_perturbed_new, &
         v_perturbed_new, heat_force, theta_new, theta_base, pressure_new, &
         pressure_base, u_perturbed_base, v_perturbed_base, w_perturbed_base
      use microphysics_perturbation, only: aerosol_base, drop_new, &
         rain_new, crystal_new, snow_new, &
         hail_new, vapor_new, aerosol_new, &
         vapor_base, drop_base, rain_base, crystal_base, snow_base, hail_base, &
         Av, Vtgra0, Vtnie
      use model_initial_conditions, only: initial_conditions
      implicit none
      integer :: unit_number

      ltt = sim_time_minutes*90.         !tiempo total de simulacion
      ltg = save_lapse_minutes*90.       !tiempo de grabacion
      lte = statistic_time_minutes*90.   !tiempo de grabacion estadistica
      ltb = backup_time_minutes*90.      !tiempo de backup

      ctur = 0.5
      pro1 = 1.-2e-2*(dt1/5.)
      pro2 = (1.-pro1)/6.
      pro3 = 1.-2e-2*(dt1/5.)
      pro4 = (1.-pro1)/4.

      total_time = nint(ltt/dt1)
      lt2 = nint(dt1/dt2)
      lt3 = 2*nint(dt1/dt3)
      cteturb = ctur/2.**.5

      cks = cteturb*2.

      dx2 = 2.*dx1
      dx8 = 8.*dx1
      dx12 = 12.*dx1
      AA = 1./Eps - 1.
      ikapa = 1./Kapa

      cteqgot = 160./3**6.*pi*rhow*N0got
      cteqllu = 8.*pi*rhow*N0llu
      cteqnie = 0.6*pi*rhonie*N0nie
      cteqgra = 8.*pi*rhogra*N0gra

      spos = 0
      posxx = 0.
      posyy = 0.
      posx(0) = 0
      posy(0) = 0
      tte = 0

      if (.not. restore_backup) then
         actual_file = 1
         call initial_conditions()
      else
         open (newunit=unit_number, file=output_directory//"inis.da", status= &
            'unknown', form='unformatted')
         read (unit_number, *) air_density_z_initial, temperature_z_initial, &
            theta_z_initial, Pres00, vapor_z_initial, cc2, &
            aerosol_z_initial, u_z_initial, v_z_initial
         close (unit_number)

         open (newunit=unit_number, file=output_directory//"velos.da", status= &
            'unknown', form='unformatted')
         rewind unit_number
         read (unit_number) u_perturbed_base, u_perturbed_new, v_perturbed_base, &
            v_perturbed_new, w_perturbed_base, w_perturbed_new, &
            theta_base, theta_new, pressure_base, pressure_new, vapor_base, &
            vapor_new, drop_base, drop_new, rain_base, rain_new, crystal_base, &
            crystal_new, snow_base, snow_new, hail_base, hail_new, aerosol_base, &
            aerosol_new, heat_force
         close (unit_number)

         open (newunit=unit_number, file=output_directory//"varconz.da", &
            status='unknown', form='unformatted')
         rewind unit_number
         read (unit_number) Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Av, Vtnie, &
            Vtgra0, vapor_z_relative, aerosol_z_relative, Eautcn, Eacrcn
         close (unit_number)
      end if

      !definicion de calores y presiones de vapor a 0 K
      Lsl00 = Lsl0*4180.

   end subroutine initialize_model

   subroutine statistics()
      use model_var
      use microphysics_perturbation
      use dinamic_var_perturbation
      use config
      integer :: unit_number
      umax = 0.
      lumax = 0
      mumax = 0
      numax = 0
      vmax = 0.
      lvmax = 0
      mvmax = 0
      nvmax = 0
      wmax = 0.
      lwmax = 0
      mwmax = 0
      nwmax = 0
      titamax = 0.
      ltitamax = 0
      mtitamax = 0
      ntitamax = 0
      qvapmax = 0.
      lqvapmax = 0
      mqvapmax = 0
      nqvapmax = 0
      qgotmax = 0.
      lqgotmax = 0
      mqgotmax = 0
      nqgotmax = 0
      qllumax = 0.
      lqllumax = 0
      mqllumax = 0
      nqllumax = 0
      qcrimax = 0.
      lqcrimax = 0
      mqcrimax = 0
      nqcrimax = 0
      qniemax = 0.
      lqniemax = 0
      mqniemax = 0
      nqniemax = 0
      qgramax = 0.
      lqgramax = 0
      mqgramax = 0
      nqgramax = 0
      aermax = 0.
      laermax = 0
      maermax = 0
      naermax = 0
      umin = 0.
      lumin = 0
      mumin = 0
      numin = 0
      vmin = 0.
      lvmin = 0
      mvmin = 0
      nvmin = 0
      wmin = 0.
      lwmin = 0
      mwmin = 0
      nwmin = 0
      titamin = 0.
      ltitamin = 0
      mtitamin = 0
      ntitamin = 0
      qvapmin = 0.
      lqvapmin = 0
      mqvapmin = 0
      nqvapmin = 0
      qgottot = 0.
      qllutot = 0.
      qcritot = 0.
      qnietot = 0.
      qgratot = 0.

      do concurrent(k=1:nz1, i=1:nx1, j=1:nx1)
         if (umax < u_perturbed_base(i, j, k)*100) then
            umax = u_perturbed_base(i, j, k)*100
            lumax = i
            mumax = j
            numax = k
         end if

         if (umin > u_perturbed_base(i, j, k)*100) then
            umin = u_perturbed_base(i, j, k)*100
            lumin = i
            mumin = j
            numin = k
         end if

         if (vmax < v_perturbed_base(i, j, k)*100) then
            vmax = v_perturbed_base(i, j, k)*100
            lvmax = i
            mvmax = j
            nvmax = k
         end if

         if (vmin > v_perturbed_base(i, j, k)*100) then
            vmin = v_perturbed_base(i, j, k)*100
            lvmin = i
            mvmin = j
            nvmin = k
         end if

         if (wmax < w_perturbed_base(i, j, k)*100) then
            wmax = w_perturbed_base(i, j, k)*100
            lwmax = i
            mwmax = j
            nwmax = k
         end if

         if (wmin > w_perturbed_base(i, j, k)*100) then
            wmin = w_perturbed_base(i, j, k)*100
            lwmin = i
            mwmin = j
            nwmin = k
         end if

         if (titamax < theta_base(i, j, k)*1000) then
            titamax = theta_base(i, j, k)*1000
            ltitamax = i
            mtitamax = j
            ntitamax = k
         end if

         if (titamin > theta_base(i, j, k)*1000) then
            titamin = theta_base(i, j, k)*1000
            ltitamin = i
            mtitamin = j
            ntitamin = k
         end if

         if (qvapmax < vapor_base(i, j, k)*1e6) then
            qvapmax = vapor_base(i, j, k)*1e6
            lqvapmax = i
            mqvapmax = j
            nqvapmax = k
         end if

         if (qvapmin > vapor_base(i, j, k)*1e6) then
            qvapmin = vapor_base(i, j, k)*1e6
            lqvapmin = i
            mqvapmin = j
            nqvapmin = k
         end if

         if (qgotmax < drop_base(i, j, k)*1e6) then
            qgotmax = drop_base(i, j, k)*1e6
            lqgotmax = i
            mqgotmax = j
            nqgotmax = k
         end if
         qgottot = qgottot + drop_base(i, j, k)*1e6

         if (qllumax < rain_base(i, j, k)*1e6) then
            qllumax = rain_base(i, j, k)*1e6
            lqllumax = i
            mqllumax = j
            nqllumax = k
         end if
         qllutot = qllutot + rain_base(i, j, k)*1e6

         if (qcrimax < crystal_base(i, j, k)*1e6) then
            qcrimax = crystal_base(i, j, k)*1e6
            lqcrimax = i
            mqcrimax = j
            nqcrimax = k
         end if
         qcritot = qcritot + crystal_base(i, j, k)*1e6

         if (qniemax < snow_base(i, j, k)*1e6) then
            qniemax = snow_base(i, j, k)*1e6
            lqniemax = i
            mqniemax = j
            nqniemax = k
         end if
         qnietot = qnietot + snow_base(i, j, k)*1e6

         if (qgramax < hail_base(i, j, k)*1e6) then
            qgramax = hail_base(i, j, k)*1e6
            lqgramax = i
            mqgramax = j
            nqgramax = k
         end if
         qgratot = qgratot + hail_base(i, j, k)*1e6

         if (aermax < aerosol_base(i, j, k)/1000) then
            aermax = aerosol_base(i, j, k)/1000
            laermax = i
            maermax = j
            naermax = k
         end if
      end do

      qgotmax = 0.
      qcrimax = 0.
      qniemax = 0.

      do concurrent(i=-1:1, j=-1:1, k=-1:1)
         qgotmax = qgotmax + 1e5*drop_base(lqgotmax + i, mqgotmax + j, nqgotmax + k)
         qcrimax = qcrimax + 1e5*crystal_base(lqcrimax + i, mqcrimax + j, nqcrimax + k)
         qniemax = qniemax + 1e5*snow_base(lqniemax + i, mqniemax + j, nqniemax + k)
      end do

      qgotmax = qgotmax/27.
      qcrimax = qcrimax/27.
      qniemax = qniemax/27.
      umax = umax/10
      umin = umin/10
      vmax = vmax/10
      vmin = vmin/10
      wmax = wmax/10
      wmin = wmin/10
      titamax = titamax/10
      titamin = titamin/10
      qvapmax = qvapmax/10
      qvapmin = qvapmin/10
      qgottot = qgottot/1000
      qllutot = qllutot/1000
      qcritot = qcritot/1000
      qnietot = qnietot/1000
      qgratot = qgratot/1000

      open (newunit=unit_number, file=output_directory//"esta", ACCESS="append")
      write (unit_number, 710) umax, umin, vmax, vmin, wmax, wmin, titamax, titamin, &
         qvapmax, qvapmin, qgotmax, qllumax, qcrimax, qniemax, qgramax, aermax, &
         lumax, mumax, numax, lumin, mumin, numin, lvmax, mvmax, nvmax, lvmin, mvmin, &
         nvmin, lwmax, mwmax, nwmax, lwmin, mwmin, nwmin, ltitamax, mtitamax, ntitamax, &
         ltitamin, mtitamin, ntitamin, lqvapmax, mqvapmax, nqvapmax, lqvapmin, mqvapmin, &
         nqvapmin, lqgotmax, mqgotmax, nqgotmax, lqllumax, mqllumax, nqllumax, lqcrimax, &
         mqcrimax, nqcrimax, lqniemax, mqniemax, nqniemax, lqgramax, mqgramax, nqgramax, &
         laermax, maermax, naermax
      close (unit_number)

      open (newunit=unit_number, file=output_directory//"est", ACCESS="append")
      write (unit_number, 715) qgottot, qllutot, qcritot, qnietot, qgratot
      close (unit_number)
710   format(16i5, 48i4)
715   format(5i9)
   end subroutine statistics

   subroutine cloud_position()
      !! desplazamientos horizontales a partir de la velocidad media de la nube
      !! calculo la altura media de la nube, la velocidad media de la nube
      !! es tomada como la velocidad del aire sin perturbar a esa altura
      use model_var
      use lmngot
      use lmncri
      use microphysics_perturbation
      use initial_z_state
      use cant01

      if (ngot(2) >= 1 .or. ncri(2) > 1) then

         impx = 0.
         impy = 0.
         zmed = 0.
         Qaguat = 0.
         spos = 1
         Xnub(tte) = Xnub(tte - 1)
         Ynub(tte) = Ynub(tte - 1)

         laux1 = min(lgot(1), lcri(1))
         laux2 = max(lgot(2), lcri(2))
         maux1 = min(mgot(1), mcri(1))
         maux2 = max(mgot(2), mcri(2))
         naux2 = max(ngot(2), ncri(2))

         do concurrent(k=1:naux2, i=laux1:laux2, j=maux1:maux2)
            Qagua = drop_base(i, j, k) + crystal_base(i, j, k) + rain_base(i, j, k) + &
               snow_base(i, j, k) + hail_base(i, j, k)
            zmed = zmed + k*Qagua
            Qaguat = Qaguat + Qagua
         end do

         if (Qaguat > 1e-3) then
            zmed = zmed/Qaguat
            Xnub(tte) = Xnub(tte) + u_z_initial(nint(zmed))*lte
            Ynub(tte) = Ynub(tte) + v_z_initial(nint(zmed))*lte
         end if

      end if

   end subroutine cloud_position

   subroutine cloud_movement()
      !! desplazamiento de la nube
      !! Redefine el valor de todas las variables (deberian
      !! ser las que se graban solamente)
      !! calculo de la posicion media de la nube, es decir de las gotitas
      !! el centro esta inicialmente en nx1/2+.5, nx1/2+.5
      !! posxx y posyy son siempre en modulo menores que dx1
      !! En posx y posy se guarda para cada current_time la posicion en
      !! puntos de red
      use model_var
      use dinamic_var_perturbation
      use microphysics_perturbation

      if (spos == 1) then
         posxx = Xnub(tte)
         posyy = Ynub(tte)
      else
         posxx = 0.
         posyy = 0.
      end if

      posx(tte) = posx(tte - 1)
      posy(tte) = posy(tte - 1)

      !*    corrimiento en x
      if (posxx > dx1) then
         posx(tte) = posx(tte) + 1
         Xnub(tte) = Xnub(tte) - dx1

         do concurrent(k=0:nz1 + 1)
            do concurrent(j=0:nx1 + 1)
               do concurrent(i=1:nx1 + 1)
                  u_perturbed_base(i - 1, j, k) = u_perturbed_base(i, j, k)
                  v_perturbed_base(i - 1, j, k) = v_perturbed_base(i, j, k)
                  w_perturbed_base(i - 1, j, k) = w_perturbed_base(i, j, k)
                  pressure_base(i - 1, j, k) = pressure_base(i, j, k)
                  u_perturbed_new(i - 1, j, k) = u_perturbed_new(i, j, k)
                  v_perturbed_new(i - 1, j, k) = v_perturbed_new(i, j, k)
                  w_perturbed_new(i - 1, j, k) = w_perturbed_new(i, j, k)
                  pressure_new(i - 1, j, k) = pressure_new(i, j, k)
                  theta_base(i - 1, j, k) = theta_base(i, j, k)
                  vapor_base(i - 1, j, k) = vapor_base(i, j, k)
                  drop_base(i - 1, j, k) = drop_base(i, j, k)
                  rain_base(i - 1, j, k) = rain_base(i, j, k)
                  crystal_base(i - 1, j, k) = crystal_base(i, j, k)
                  aerosol_base(i - 1, j, k) = aerosol_base(i, j, k)
                  heat_force(i - 1, j, k) = heat_force(i, j, k)
               end do
               i = nx1 + 1
               u_perturbed_base(i, j, k) = u_perturbed_base(i - 1, j, k)
               v_perturbed_base(i, j, k) = v_perturbed_base(i - 1, j, k)
               w_perturbed_base(i, j, k) = w_perturbed_base(i - 1, j, k)
               pressure_base(i, j, k) = pressure_base(i - 1, j, k)
               u_perturbed_new(i, j, k) = u_perturbed_new(i - 1, j, k)
               v_perturbed_new(i, j, k) = v_perturbed_new(i - 1, j, k)
               w_perturbed_new(i, j, k) = w_perturbed_new(i - 1, j, k)
               pressure_new(i, j, k) = pressure_new(i - 1, j, k)
               theta_base(i, j, k) = theta_base(i - 1, j, k)
               vapor_base(i, j, k) = vapor_base(i - 1, j, k)
               drop_base(i, j, k) = drop_base(i - 1, j, k)
               rain_base(i, j, k) = rain_base(i - 1, j, k)
               crystal_base(i, j, k) = crystal_base(i - 1, j, k)
               aerosol_base(i, j, k) = aerosol_base(i - 1, j, k)
               heat_force(i, j, k) = 0.
            end do
            i = nx1
            j = 0
            u_perturbed_base(i, j, k) = (u_perturbed_base(i - 1, j, k) &
               + u_perturbed_base(i, j + 1, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i - 1, j, k) &
               + v_perturbed_base(i, j + 1, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i - 1, j, k) &
               + w_perturbed_base(i, j + 1, k))/2.
            pressure_base(i, j, k) = (pressure_base(i - 1, j, k) &
               + pressure_base(i, j + 1, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i - 1, j, k) &
               + u_perturbed_new(i, j + 1, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i - 1, j, k) &
               + v_perturbed_new(i, j + 1, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i - 1, j, k) &
               + w_perturbed_new(i, j + 1, k))/2.
            pressure_new(i, j, k) = (pressure_new(i - 1, j, k) &
               + pressure_new(i, j + 1, k))/2.
            theta_base(i, j, k) = (theta_base(i - 1, j, k) &
               + theta_base(i, j + 1, k))/2.
            vapor_base(i, j, k) = (vapor_base(i - 1, j, k) &
               + vapor_base(i, j + 1, k))/2.
            drop_base(i, j, k) = (drop_base(i - 1, j, k) &
               + drop_base(i, j + 1, k))/2.
            rain_base(i, j, k) = (rain_base(i - 1, j, k) &
               + rain_base(i, j + 1, k))/2.
            crystal_base(i, j, k) = (crystal_base(i - 1, j, k) &
               + crystal_base(i, j + 1, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i - 1, j, k) &
               + aerosol_base(i, j + 1, k))/2.
            heat_force(i, j, k) = 0.
            j = nx1 + 1
            u_perturbed_base(i, j, k) = (u_perturbed_base(i - 1, j, k) + u_perturbed_base(i, j - 1, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i - 1, j, k) + v_perturbed_base(i, j - 1, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i - 1, j, k) + w_perturbed_base(i, j - 1, k))/2.
            pressure_base(i, j, k) = (pressure_base(i - 1, j, k) + pressure_base(i, j - 1, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i - 1, j, k) + u_perturbed_new(i, j - 1, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i - 1, j, k) + v_perturbed_new(i, j - 1, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i - 1, j, k) + w_perturbed_new(i, j - 1, k))/2.
            pressure_new(i, j, k) = (pressure_new(i - 1, j, k) + pressure_new(i, j - 1, k))/2.
            theta_base(i, j, k) = (theta_base(i - 1, j, k) + theta_base(i, j - 1, k))/2.
            vapor_base(i, j, k) = (vapor_base(i - 1, j, k) + vapor_base(i, j - 1, k))/2.
            drop_base(i, j, k) = (drop_base(i - 1, j, k) + drop_base(i, j - 1, k))/2.
            rain_base(i, j, k) = (rain_base(i - 1, j, k) + rain_base(i, j - 1, k))/2.
            crystal_base(i, j, k) = (crystal_base(i - 1, j, k) + crystal_base(i, j - 1, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i - 1, j, k) + aerosol_base(i, j - 1, k))/2.
            heat_force(i, j, k) = 0.
         end do
      end if

      if (posxx < -dx1) then
         posx(tte) = posx(tte) - 1
         Xnub(tte) = Xnub(tte) + dx1

         do concurrent(k=0:nz1 + 1)
            do concurrent(j=0:nx1 + 1)
               do concurrent(i=nx1:0) ! TODO Test this loop
                  u_perturbed_base(i + 1, j, k) = u_perturbed_base(i, j, k)
                  v_perturbed_base(i + 1, j, k) = v_perturbed_base(i, j, k)
                  w_perturbed_base(i + 1, j, k) = w_perturbed_base(i, j, k)
                  pressure_base(i + 1, j, k) = pressure_base(i, j, k)
                  u_perturbed_new(i + 1, j, k) = u_perturbed_new(i, j, k)
                  v_perturbed_new(i + 1, j, k) = v_perturbed_new(i, j, k)
                  w_perturbed_new(i + 1, j, k) = w_perturbed_new(i, j, k)
                  pressure_new(i + 1, j, k) = pressure_new(i, j, k)
                  theta_base(i + 1, j, k) = theta_base(i, j, k)
                  vapor_base(i + 1, j, k) = vapor_base(i, j, k)
                  drop_base(i + 1, j, k) = drop_base(i, j, k)
                  rain_base(i + 1, j, k) = rain_base(i, j, k)
                  crystal_base(i + 1, j, k) = crystal_base(i, j, k)
                  aerosol_base(i + 1, j, k) = aerosol_base(i, j, k)
                  heat_force(i + 1, j, k) = heat_force(i, j, k)
               end do
               i = 0
               u_perturbed_base(i, j, k) = u_perturbed_base(i + 1, j, k)
               v_perturbed_base(i, j, k) = v_perturbed_base(i + 1, j, k)
               w_perturbed_base(i, j, k) = w_perturbed_base(i + 1, j, k)
               pressure_base(i, j, k) = pressure_base(i + 1, j, k)
               u_perturbed_new(i, j, k) = u_perturbed_new(i + 1, j, k)
               v_perturbed_new(i, j, k) = v_perturbed_new(i + 1, j, k)
               w_perturbed_new(i, j, k) = w_perturbed_new(i + 1, j, k)
               pressure_new(i, j, k) = pressure_new(i + 1, j, k)
               theta_base(i, j, k) = theta_base(i + 1, j, k)
               vapor_base(i, j, k) = vapor_base(i + 1, j, k)
               drop_base(i, j, k) = drop_base(i + 1, j, k)
               rain_base(i, j, k) = rain_base(i + 1, j, k)
               crystal_base(i, j, k) = crystal_base(i + 1, j, k)
               aerosol_base(i, j, k) = aerosol_base(i + 1, j, k)
               heat_force(i, j, k) = 0.
            end do
            i = 1
            j = 0
            u_perturbed_base(i, j, k) = (u_perturbed_base(i + 1, j, k) &
               + u_perturbed_base(i, j + 1, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i + 1, j, k) &
               + v_perturbed_base(i, j + 1, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i + 1, j, k) &
               + w_perturbed_base(i, j + 1, k))/2.
            pressure_base(i, j, k) = (pressure_base(i + 1, j, k) &
               + pressure_base(i, j + 1, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i + 1, j, k) &
               + u_perturbed_new(i, j + 1, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i + 1, j, k) &
               + v_perturbed_new(i, j + 1, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i + 1, j, k) &
               + w_perturbed_new(i, j + 1, k))/2.
            pressure_new(i, j, k) = (pressure_new(i + 1, j, k) &
               + pressure_new(i, j + 1, k))/2.
            theta_base(i, j, k) = (theta_base(i + 1, j, k) &
               + theta_base(i, j + 1, k))/2.
            vapor_base(i, j, k) = (vapor_base(i + 1, j, k) &
               + vapor_base(i, j + 1, k))/2.
            drop_base(i, j, k) = (drop_base(i + 1, j, k) &
               + drop_base(i, j + 1, k))/2.
            rain_base(i, j, k) = (rain_base(i + 1, j, k) &
               + rain_base(i, j + 1, k))/2.
            crystal_base(i, j, k) = (crystal_base(i + 1, j, k) &
               + crystal_base(i, j + 1, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i + 1, j, k) &
               + aerosol_base(i, j + 1, k))/2.
            heat_force(i, j, k) = 0.
            j = nx1 + 1
            u_perturbed_base(i, j, k) = (u_perturbed_base(i + 1, j, k) &
               + u_perturbed_base(i, j - 1, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i + 1, j, k) &
               + v_perturbed_base(i, j - 1, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i + 1, j, k) &
               + w_perturbed_base(i, j - 1, k))/2.
            pressure_base(i, j, k) = (pressure_base(i + 1, j, k) &
               + pressure_base(i, j - 1, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i + 1, j, k) &
               + u_perturbed_new(i, j - 1, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i + 1, j, k) &
               + v_perturbed_new(i, j - 1, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i + 1, j, k) &
               + w_perturbed_new(i, j - 1, k))/2.
            pressure_new(i, j, k) = (pressure_new(i + 1, j, k) &
               + pressure_new(i, j - 1, k))/2.
            theta_base(i, j, k) = (theta_base(i + 1, j, k) &
               + theta_base(i, j - 1, k))/2.
            vapor_base(i, j, k) = (vapor_base(i + 1, j, k) &
               + vapor_base(i, j - 1, k))/2.
            drop_base(i, j, k) = (drop_base(i + 1, j, k) &
               + drop_base(i, j - 1, k))/2.
            rain_base(i, j, k) = (rain_base(i + 1, j, k) &
               + rain_base(i, j - 1, k))/2.
            crystal_base(i, j, k) = (crystal_base(i + 1, j, k) &
               + crystal_base(i, j - 1, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i + 1, j, k) &
               + aerosol_base(i, j - 1, k))/2.
            heat_force(i, j, k) = 0.
         end do
      end if

      !*    corrimiento en y

      if (posyy > dx1) then
         posy(tte) = posy(tte) + 1
         Ynub(tte) = Ynub(tte) - dx1

         do concurrent(k=0:nz1 + 1)
            do concurrent(i=0:nx1 + 1)
               do concurrent(j=1:nx1 + 1)
                  u_perturbed_base(i, j - 1, k) = u_perturbed_base(i, j, k)
                  v_perturbed_base(i, j - 1, k) = v_perturbed_base(i, j, k)
                  w_perturbed_base(i, j - 1, k) = w_perturbed_base(i, j, k)
                  pressure_base(i, j - 1, k) = pressure_base(i, j, k)
                  u_perturbed_new(i, j - 1, k) = u_perturbed_new(i, j, k)
                  v_perturbed_new(i, j - 1, k) = v_perturbed_new(i, j, k)
                  w_perturbed_new(i, j - 1, k) = w_perturbed_new(i, j, k)
                  pressure_new(i, j - 1, k) = pressure_new(i, j, k)
                  theta_base(i, j - 1, k) = theta_base(i, j, k)
                  vapor_base(i, j - 1, k) = vapor_base(i, j, k)
                  drop_base(i, j - 1, k) = drop_base(i, j, k)
                  rain_base(i, j - 1, k) = rain_base(i, j, k)
                  crystal_base(i, j - 1, k) = crystal_base(i, j, k)
                  aerosol_base(i, j - 1, k) = aerosol_base(i, j, k)
                  heat_force(i, j - 1, k) = heat_force(i, j, k)
               end do
               j = nx1 + 1
               u_perturbed_base(i, j, k) = u_perturbed_base(i, j - 1, k)
               v_perturbed_base(i, j, k) = v_perturbed_base(i, j - 1, k)
               w_perturbed_base(i, j, k) = w_perturbed_base(i, j - 1, k)
               pressure_base(i, j, k) = pressure_base(i, j - 1, k)
               u_perturbed_new(i, j, k) = u_perturbed_new(i, j - 1, k)
               v_perturbed_new(i, j, k) = v_perturbed_new(i, j - 1, k)
               w_perturbed_new(i, j, k) = w_perturbed_new(i, j - 1, k)
               pressure_new(i, j, k) = pressure_new(i, j - 1, k)
               theta_base(i, j, k) = theta_base(i, j - 1, k)
               vapor_base(i, j, k) = vapor_base(i, j - 1, k)
               drop_base(i, j, k) = drop_base(i, j - 1, k)
               rain_base(i, j, k) = rain_base(i, j - 1, k)
               crystal_base(i, j, k) = crystal_base(i, j - 1, k)
               aerosol_base(i, j, k) = aerosol_base(i, j - 1, k)
               heat_force(i, j, k) = 0.
            end do
            j = nx1
            i = 0
            u_perturbed_base(i, j, k) = (u_perturbed_base(i, j - 1, k) &
               + u_perturbed_base(i + 1, j, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i, j - 1, k) &
               + v_perturbed_base(i + 1, j, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i, j - 1, k) &
               + w_perturbed_base(i + 1, j, k))/2.
            pressure_base(i, j, k) = (pressure_base(i, j - 1, k) &
               + pressure_base(i + 1, j, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i, j - 1, k) &
               + u_perturbed_new(i + 1, j, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i, j - 1, k) &
               + v_perturbed_new(i + 1, j, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i, j - 1, k) &
               + w_perturbed_new(i + 1, j, k))/2.
            pressure_new(i, j, k) = (pressure_new(i, j - 1, k) &
               + pressure_new(i + 1, j, k))/2.
            theta_base(i, j, k) = (theta_base(i, j - 1, k) &
               + theta_base(i + 1, j, k))/2.
            vapor_base(i, j, k) = (vapor_base(i, j - 1, k) &
               + vapor_base(i + 1, j, k))/2.
            drop_base(i, j, k) = (drop_base(i, j - 1, k) &
               + drop_base(i + 1, j, k))/2.
            rain_base(i, j, k) = (rain_base(i, j - 1, k) &
               + rain_base(i + 1, j, k))/2.
            crystal_base(i, j, k) = (crystal_base(i, j - 1, k) &
               + crystal_base(i + 1, j, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i, j - 1, k) &
               + aerosol_base(i + 1, j, k))/2.
            heat_force(i, j, k) = 0.
            i = nx1 + 1
            u_perturbed_base(i, j, k) = (u_perturbed_base(i, j - 1, k) &
               + u_perturbed_base(i - 1, j, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i, j - 1, k) &
               + v_perturbed_base(i - 1, j, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i, j - 1, k) &
               + w_perturbed_base(i - 1, j, k))/2.
            pressure_base(i, j, k) = (pressure_base(i, j - 1, k) &
               + pressure_base(i - 1, j, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i, j - 1, k) &
               + u_perturbed_new(i - 1, j, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i, j - 1, k) &
               + v_perturbed_new(i - 1, j, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i, j - 1, k) &
               + w_perturbed_new(i - 1, j, k))/2.
            pressure_new(i, j, k) = (pressure_new(i, j - 1, k) &
               + pressure_new(i - 1, j, k))/2.
            theta_base(i, j, k) = (theta_base(i, j - 1, k) &
               + theta_base(i - 1, j, k))/2.
            vapor_base(i, j, k) = (vapor_base(i, j - 1, k) &
               + vapor_base(i - 1, j, k))/2.
            drop_base(i, j, k) = (drop_base(i, j - 1, k) &
               + drop_base(i - 1, j, k))/2.
            rain_base(i, j, k) = (rain_base(i, j - 1, k) &
               + rain_base(i - 1, j, k))/2.
            crystal_base(i, j, k) = (crystal_base(i, j - 1, k) &
               + crystal_base(i - 1, j, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i, j - 1, k) &
               + aerosol_base(i - 1, j, k))/2.
            heat_force(i, j, k) = 0.
         end do
      end if

      if (posyy < -dx1) then
         posy(tte) = posy(tte) - 1
         Xnub(tte) = Xnub(tte) + dx1

         do concurrent(k=0:nz1 + 1)
            do concurrent(i=0:nx1 + 1)
               do concurrent(j=nx1:0) ! TODO Test this loop
                  u_perturbed_base(i, j + 1, k) = u_perturbed_base(i, j, k)
                  v_perturbed_base(i, j + 1, k) = v_perturbed_base(i, j, k)
                  w_perturbed_base(i, j + 1, k) = w_perturbed_base(i, j, k)
                  pressure_base(i, j + 1, k) = pressure_base(i, j, k)
                  u_perturbed_new(i, j + 1, k) = u_perturbed_new(i, j, k)
                  v_perturbed_new(i, j + 1, k) = v_perturbed_new(i, j, k)
                  w_perturbed_new(i, j + 1, k) = w_perturbed_new(i, j, k)
                  pressure_new(i, j + 1, k) = pressure_new(i, j, k)
                  theta_base(i, j + 1, k) = theta_base(i, j, k)
                  vapor_base(i, j + 1, k) = vapor_base(i, j, k)
                  drop_base(i, j + 1, k) = drop_base(i, j, k)
                  rain_base(i, j + 1, k) = rain_base(i, j, k)
                  crystal_base(i, j + 1, k) = crystal_base(i, j, k)
                  aerosol_base(i, j + 1, k) = aerosol_base(i, j, k)
                  heat_force(i, j + 1, k) = heat_force(i, j, k)
               end do
               j = 0
               u_perturbed_base(i, j, k) = u_perturbed_base(i, j - 1, k)
               v_perturbed_base(i, j, k) = v_perturbed_base(i, j - 1, k)
               w_perturbed_base(i, j, k) = w_perturbed_base(i, j - 1, k)
               pressure_base(i, j, k) = pressure_base(i, j - 1, k)
               u_perturbed_new(i, j, k) = u_perturbed_new(i, j - 1, k)
               v_perturbed_new(i, j, k) = v_perturbed_new(i, j - 1, k)
               w_perturbed_new(i, j, k) = w_perturbed_new(i, j - 1, k)
               pressure_new(i, j, k) = pressure_new(i, j - 1, k)
               theta_base(i, j, k) = theta_base(i, j - 1, k)
               vapor_base(i, j, k) = vapor_base(i, j - 1, k)
               drop_base(i, j, k) = drop_base(i, j - 1, k)
               rain_base(i, j, k) = rain_base(i, j - 1, k)
               crystal_base(i, j, k) = crystal_base(i, j - 1, k)
               aerosol_base(i, j, k) = aerosol_base(i, j - 1, k)
               heat_force(i, j, k) = 0.
            end do
            j = 1
            i = 0
            u_perturbed_base(i, j, k) = (u_perturbed_base(i, j + 1, k) &
               + u_perturbed_base(i + 1, j, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i, j + 1, k) &
               + v_perturbed_base(i + 1, j, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i, j + 1, k) &
               + w_perturbed_base(i + 1, j, k))/2.
            pressure_base(i, j, k) = (pressure_base(i, j + 1, k) &
               + pressure_base(i + 1, j, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i, j + 1, k) &
               + u_perturbed_new(i + 1, j, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i, j + 1, k) &
               + v_perturbed_new(i + 1, j, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i, j + 1, k) &
               + w_perturbed_new(i + 1, j, k))/2.
            pressure_new(i, j, k) = (pressure_new(i, j + 1, k) &
               + pressure_new(i + 1, j, k))/2.
            theta_base(i, j, k) = (theta_base(i, j + 1, k) &
               + theta_base(i + 1, j, k))/2.
            vapor_base(i, j, k) = (vapor_base(i, j + 1, k) &
               + vapor_base(i + 1, j, k))/2.
            drop_base(i, j, k) = (drop_base(i, j + 1, k) &
               + drop_base(i + 1, j, k))/2.
            rain_base(i, j, k) = (rain_base(i, j + 1, k) &
               + rain_base(i + 1, j, k))/2.
            crystal_base(i, j, k) = (crystal_base(i, j + 1, k) &
               + crystal_base(i + 1, j, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i, j + 1, k) &
               + aerosol_base(i + 1, j, k))/2.
            heat_force(i, j, k) = 0.

            i = nx1 + 1
            u_perturbed_base(i, j, k) = (u_perturbed_base(i, j + 1, k) &
               + u_perturbed_base(i - 1, j, k))/2.
            v_perturbed_base(i, j, k) = (v_perturbed_base(i, j + 1, k) &
               + v_perturbed_base(i - 1, j, k))/2.
            w_perturbed_base(i, j, k) = (w_perturbed_base(i, j + 1, k) &
               + w_perturbed_base(i - 1, j, k))/2.
            pressure_base(i, j, k) = (pressure_base(i, j + 1, k) &
               + pressure_base(i - 1, j, k))/2.
            u_perturbed_new(i, j, k) = (u_perturbed_new(i, j + 1, k) &
               + u_perturbed_new(i - 1, j, k))/2.
            v_perturbed_new(i, j, k) = (v_perturbed_new(i, j + 1, k) &
               + v_perturbed_new(i - 1, j, k))/2.
            w_perturbed_new(i, j, k) = (w_perturbed_new(i, j + 1, k) &
               + w_perturbed_new(i - 1, j, k))/2.
            pressure_new(i, j, k) = (pressure_new(i, j + 1, k) &
               + pressure_new(i - 1, j, k))/2.
            theta_base(i, j, k) = (theta_base(i, j + 1, k) &
               + theta_base(i - 1, j, k))/2.
            vapor_base(i, j, k) = (vapor_base(i, j + 1, k) &
               + vapor_base(i - 1, j, k))/2.
            drop_base(i, j, k) = (drop_base(i, j + 1, k) &
               + drop_base(i - 1, j, k))/2.
            rain_base(i, j, k) = (rain_base(i, j + 1, k) &
               + rain_base(i - 1, j, k))/2.
            crystal_base(i, j, k) = (crystal_base(i, j + 1, k) &
               + crystal_base(i - 1, j, k))/2.
            aerosol_base(i, j, k) = (aerosol_base(i, j + 1, k) &
               + aerosol_base(i - 1, j, k))/2.
            heat_force(i, j, k) = 0.
         end do

      end if

      posxx = posx(tte)*dx1 + Xnub(tte)
      posyy = posy(tte)*dx1 + Ynub(tte)

   end subroutine cloud_movement
end module model_initialization

