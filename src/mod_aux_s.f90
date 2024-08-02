module extra_subrut
contains
   subroutine turbu1(kk)
      !! Esta subrutina calcula los Dnm para cada plano Z
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use constants
      use turbvar1
      use turbu1_vars
      implicit none

      integer, intent(in) :: kk

      k = kk + 1
      do concurrent(i=0:nx1 + 1, j=0:nx1 + 1)
         if (kk == 1) then
            do concurrent(n=1:2)
               do concurrent(m=1:n)
                  D(n, m, i, j, 1) = 0.
               end do
            end do
            D(3, 1, i, j, 1) = u_perturbed_new(i, j, 1)
            D(3, 2, i, j, 1) = v_perturbed_new(i, j, 1)
            D(3, 3, i, j, 1) = w_perturbed_new(i, j, 1)*2./3.
            do concurrent(n=1:3)
               do concurrent(m=1:n)
                  D(m, n, i, j, 1) = D(n, m, i, j, 1)
               end do
            end do
            do concurrent(lx=-1:1, ly=-1:1, lz=-1:1)
               ldis = abs(lx) + abs(ly) + abs(lz)
               if (ldis <= 1) then
                  vel(1, lx, ly, lz) = u_perturbed_new(lx + i, ly + j, lz + 1)
                  vel(2, lx, ly, lz) = v_perturbed_new(lx + i, ly + j, lz + 1)
                  vel(3, lx, ly, lz) = w_perturbed_new(lx + i, ly + j, lz + 1)
               end if
            end do
            !     calculo de Dij
            do concurrent(n=1:3)
               dv(n, 1) = vel(n, 1, 0, 0) - vel(n, -1, 0, 0)
               dv(n, 2) = vel(n, 0, 1, 0) - vel(n, 0, -1, 0)
               dv(n, 3) = vel(n, 0, 0, 1) - vel(n, 0, 0, -1)
            end do
            do concurrent(n=1:3)
               do concurrent(m=1:n)
                  D(n, m, i, j, 2) = (dv(n, m) + dv(m, n))
                  D(m, n, i, j, 2) = D(n, m, i, j, 2)
                  if (n == m) D(n, n, i, j, 2) = 2./3.*D(n, n, i, j, 2)
               end do
            end do
         else
            do concurrent(n=1:3, m=1:3, lz=1:2)
               D(n, m, i, j, lz) = D(n, m, i, j, lz + 1)
            end do
         end if
         !*********************************************************

         !     Lectura de las velocidades necesarias
         do concurrent(lx=-1:1, ly=-1:1, lz=-1:1)
            ldis = abs(lx) + abs(ly) + abs(lz)
            if (ldis <= 1) then
               vel(1, lx, ly, lz) = u_perturbed_new(lx + i, ly + j, lz + k)
               vel(2, lx, ly, lz) = v_perturbed_new(lx + i, ly + j, lz + k)
               vel(3, lx, ly, lz) = w_perturbed_new(lx + i, ly + j, lz + k)
            end if
         end do
         !     calculo de Dij
         do concurrent(n=1:3)
            dv(n, 1) = vel(n, 1, 0, 0) - vel(n, -1, 0, 0)
            dv(n, 2) = vel(n, 0, 1, 0) - vel(n, 0, -1, 0)
            dv(n, 3) = vel(n, 0, 0, 1) - vel(n, 0, 0, -1)
         end do
         do concurrent(n=1:3)
            do concurrent(m=1:n)
               D(n, m, i, j, 3) = (dv(n, m) + dv(m, n))
               D(m, n, i, j, 3) = D(n, m, i, j, 3)
               if (n == m) D(n, n, i, j, 3) = 2./3.*D(n, n, i, j, 3)
            end do
         end do
      end do
      return
   end subroutine turbu1

   subroutine turbu2(i, j)
      !!      Esta subrutina calcula las cantidades referida a los terminos de
      !!      turbulencia: K, DK, DDij
      !!      Dij viene de turbu1
      !!     En realidad para tener los valores de las 4 cantidades falta
      !!     multiplicarlas por:
      !!                      KMM : cteturb*dx1/2
      !!                      DK  : cteturb/4
      !!                      Dij : 1/dx2 (simetrico)
      !!                      DDij: 1/dx2**2
      use dimensions
      use turbvar
      use turbvar1
      use turbu2_vars
      implicit none
      integer, intent(in) :: i, j
      real aux

      !     calculo de KM
      do lx = -1, 1
         do ly = -1, 1
            do lz = -1, 1
               ldis = abs(lx) + abs(ly) + abs(lz)
               if (ldis <= 2) then
                  call suma(aux, &
                            D(1, 1, i + lx, j + ly, 2 + lz)**2., &
                            D(2, 2, i + lx, j + ly, 2 + lz)**2., &
                            D(3, 3, i + lx, j + ly, 2 + lz)**2.)
                  sum = aux/2.
                  call suma(aux, &
                            D(1, 2, i + lx, j + ly, 2 + lz)**2., &
                            D(1, 3, i + lx, j + ly, 2 + lz)**2., &
                            D(2, 3, i + lx, j + ly, 2 + lz)**2.)
                  sum = sum + aux

                  KM(lx, ly, lz) = sum**.5
               end if
            end do
         end do
      end do

      !     calculo de las derivadas (sin la distancia abajo)
      KM1 = KM(1, 0, 0) - KM(-1, 0, 0)
      KM2 = KM(0, 1, 0) - KM(0, -1, 0)
      KM3 = KM(0, 0, 1) - KM(0, 0, -1)
      do concurrent(n=1:3)
         D1(n) = D(n, 1, i + 1, j, 2) - D(n, 1, i - 1, j, 2)
         D2(n) = D(n, 2, i, j + 1, 2) - D(n, 2, i, j - 1, 2)
         D3(n) = D(n, 3, i, j, 3) - D(n, 3, i, j, 1)
      end do
      KMM = KM(0, 0, 0)
      do concurrent(n=1:3, m=1:3)
         DD(n, m) = D(n, m, i, j, 2)
      end do

      return
   end subroutine turbu2

   subroutine inhomogeneous_velocities(i, j, k, dden0z)
      !! Esta subrutina calcula los terminos inomogeneos para las velocidades
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use turbvar
      use sv_inhomogeneous_velocities_and_speed_pressure, only: fu, fv, fw, fp
      implicit none

      integer, intent(in) :: i, j, k
      real, intent(in) :: dden0z

      real(8) :: velocity_xx !! Velocity component xx
      real(8) :: velocity_xy !! Velocity component xy
      real(8) :: velocity_xz !! Velocity component xz
      real(8) :: velocity_yx !! Velocity component yx
      real(8) :: velocity_yy !! Velocity component yy
      real(8) :: velocity_yz !! Velocity component yz
      real(8) :: velocity_zx !! Velocity component zx
      real(8) :: velocity_zy !! Velocity component zy
      real(8) :: velocity_zz !! Velocity component zz

      real(8) :: coefficient_a1 !! Coefficient a1
      real(8) :: coefficient_a2 !! Coefficient a2
      real(8) :: coefficient_a3 !! Coefficient a3
      real(8) :: turbulence_x !! Turbulence x
      real(8) :: turbulence_y !! Turbulence y
      real(8) :: turbulence_z !! Turbulence z

      real(8) :: divergence_x !! Divergence x
      real(8) :: divergence_y !! Divergence y
      real(8) :: divergence_z !! Divergence z
      real(8) :: gravitational_acceleration !! Gravitational acceleration
      real(8) :: laplacian_of_laplacian !! Laplacian of laplacian

      velocity_xx = (u_perturbed_new(i - 2, j, k) - u_perturbed_new(i + 2, j, k)) &
                    + 8.*(u_perturbed_new(i + 1, j, k) - u_perturbed_new(i - 1, j, k))
      velocity_xy = (u_perturbed_new(i, j - 2, k) - u_perturbed_new(i, j + 2, k)) &
                    + 8.*(u_perturbed_new(i, j + 1, k) - u_perturbed_new(i, j - 1, k))
      velocity_xz = (u_perturbed_new(i, j, k - 2) - u_perturbed_new(i, j, k + 2)) &
                    + 8.*(u_perturbed_new(i, j, k + 1) - u_perturbed_new(i, j, k - 1))

      velocity_yx = (v_perturbed_new(i - 2, j, k) - v_perturbed_new(i + 2, j, k)) &
                    + 8.*(v_perturbed_new(i + 1, j, k) - v_perturbed_new(i - 1, j, k))
      velocity_yy = (v_perturbed_new(i, j - 2, k) - v_perturbed_new(i, j + 2, k)) &
                    + 8.*(v_perturbed_new(i, j + 1, k) - v_perturbed_new(i, j - 1, k))
      velocity_yz = (v_perturbed_new(i, j, k - 2) - v_perturbed_new(i, j, k + 2)) &
                    + 8.*(v_perturbed_new(i, j, k + 1) - v_perturbed_new(i, j, k - 1))

      velocity_zx = (w_perturbed_new(i - 2, j, k) - w_perturbed_new(i + 2, j, k)) &
                    + 8.*(w_perturbed_new(i + 1, j, k) - w_perturbed_new(i - 1, j, k))
      velocity_zy = (w_perturbed_new(i, j - 2, k) - w_perturbed_new(i, j + 2, k)) &
                    + 8.*(w_perturbed_new(i, j + 1, k) - w_perturbed_new(i, j - 1, k))
      velocity_zz = (w_perturbed_new(i, j, k - 2) - w_perturbed_new(i, j, k + 2)) &
                    + 8.*(w_perturbed_new(i, j, k + 1) - w_perturbed_new(i, j, k - 1))

      divergence_x = (u_perturbed_new(i, j, k)*velocity_xx + v_perturbed_new(i, j, k)*velocity_xy) &
                     + w_perturbed_new(i, j, k)*velocity_xz
      divergence_y = (u_perturbed_new(i, j, k)*velocity_yx + v_perturbed_new(i, j, k)*velocity_yy) &
                     + w_perturbed_new(i, j, k)*velocity_yz
      divergence_z = (u_perturbed_new(i, j, k)*velocity_zx + v_perturbed_new(i, j, k)*velocity_zy) &
                     + w_perturbed_new(i, j, k)*velocity_zz

      coefficient_a1 = (KM1*DD(1, 1) + KM2*DD(1, 2)) + KM3*DD(1, 3)
      coefficient_a2 = KMM*((D1(1) + D2(1)) + D3(1))
      coefficient_a3 = KMM*DD(1, 3)*dden0z
      turbulence_x = cteturb*((coefficient_a1 + coefficient_a2) + coefficient_a3)

      coefficient_a1 = (KM1*DD(2, 1) + KM2*DD(2, 2)) + KM3*DD(2, 3)
      coefficient_a2 = KMM*((D1(2) + D2(2)) + D3(2))
      coefficient_a3 = KMM*DD(2, 3)*dden0z
      turbulence_y = cteturb*((coefficient_a1 + coefficient_a2) + coefficient_a3)

      coefficient_a1 = (KM1*DD(3, 1) + KM2*DD(3, 2)) + KM3*DD(3, 3)
      coefficient_a2 = KMM*((D1(3) + D2(3)) + D3(3))
      coefficient_a3 = KMM*DD(3, 3)*dden0z
      turbulence_z = cteturb*((coefficient_a1 + coefficient_a2) + coefficient_a3)

      !$$
      gravitational_acceleration = G*(theta_base(i, j, k)/theta_z_initial(k) &
                                      + (AA*vapor_base(i, j, k) - drop_base(i, j, k) - rain_base(i, j, k) &
                                         - crystal_base(i, j, k) - snow_base(i, j, k) - hail_base(i, j, k) &
                                         )/air_density_z_initial(k))

      fu(i, j, k) = turbulence_x/dx8 - divergence_x/dx12
      fv(i, j, k) = turbulence_y/dx8 - divergence_y/dx12
      fw(i, j, k) = turbulence_z/dx8 - divergence_z/dx12 + gravitational_acceleration

      !     agregado para la P (23/8/97)
      laplacian_of_laplacian = (pressure_new(i + 1, j, k) + pressure_new(i - 1, j, k)) &
                               + (pressure_new(i, j + 1, k) + pressure_new(i, j - 1, k)) &
                               + pressure_new(i, j, k + 1) + pressure_new(i, j, k - 1) &
                               - 6.*pressure_new(i, j, k)

      fp(i, j, k) = cteturb*KMM/dx2*laplacian_of_laplacian

      return
   end subroutine inhomogeneous_velocities

   subroutine tempot(i, j, k, dden0z, Fcal)
      !! heat_force es el calor liberado por cambio de fase, por unidad de masa de aire
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use constants
      use initial_z_state
      use turbvar
      use tempe01
      implicit none

      real, intent(in) :: dden0z, Fcal
      integer, intent(in) :: i, j, k

      adv(1) = (u_perturbed_new(i, j, k) + u_z_initial(k))*(theta_base(i + 1, j, k) - theta_base(i - 1, j, k))
      adv(2) = (v_perturbed_new(i, j, k) + v_z_initial(k))*(theta_base(i, j + 1, k) - theta_base(i, j - 1, k))
      adv(3) = w_perturbed_new(i, j, k)*(theta_base(i, j, k + 1) - theta_base(i, j, k - 1))

      advec = -((adv(1) + adv(2)) + adv(3))

      verti = -w_perturbed_new(i, j, k)*(theta_z_initial(k + 1) - theta_z_initial(k - 1))

      calor = Fcal*theta_z_initial(k)/(temperature_z_initial(k)*Cp)

      dtita(1) = theta_base(i + 1, j, k) - theta_base(i - 1, j, k)
      dtita(2) = theta_base(i, j + 1, k) - theta_base(i, j - 1, k)
      dtita(3) = theta_base(i, j, k + 1) - theta_base(i, j, k - 1)

      escal = (dtita(1)*KM1 + dtita(2)*KM2) + dtita(3)*KM3

      lapla = (theta_base(i + 1, j, k) + theta_base(i - 1, j, k)) &
              + (theta_base(i, j + 1, k) + theta_base(i, j - 1, k)) &
              + theta_base(i, j, k + 1) + theta_base(i, j, k - 1) &
              - 6*theta_base(i, j, k)
      lapla = lapla + (theta_z_initial(k + 1) + theta_z_initial(k - 1) - 2.*theta_z_initial(k))

      turden = dden0z*(theta_z_initial(k + 1) - theta_z_initial(k - 1))

      turbul = 3.*cteturb/dx8*(escal + KMM*(4.*lapla + turden))

      theta_new(i, j, k) = dt1*((advec + verti)/dx2 + turbul + calor) + &
                           theta_base(i, j, k)
      return
   end subroutine tempot

   subroutine dvapor(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dvapor_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqvap(1) = vapor_base(l + 1, m, n) - vapor_base(l - 1, m, n)
      dqvap(2) = vapor_base(l, m + 1, n) - vapor_base(l, m - 1, n)
      dqvap(3) = vapor_base(l, m, n + 1) - vapor_base(l, m, n - 1)

      adv(1) = (((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))* &
                 (vapor_base(l + 1, m, n) + vapor_base(l, m, n))) - &
                ((u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))* &
                 (vapor_base(l - 1, m, n) + vapor_base(l, m, n))))/4.
      adv(1) = adv(1) + dqvap(1)/2.*u_z_initial(n)

      adv(2) = (((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))* &
                 (vapor_base(l, m + 1, n) + vapor_base(l, m, n))) - &
                ((v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))* &
                 (vapor_base(l, m - 1, n) + vapor_base(l, m, n))))/4.
      adv(2) = adv(2) + dqvap(2)/2.*v_z_initial(n)

      advvap2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (vapor_base(l, m, n) + vapor_base(l, m, n + 1))/4.

      adv(3) = advvap2(l, m) - advvap1(l, m)

      advec = -((adv(1) + adv(2)) + adv(3))

      verti = -((w_perturbed_new(l, m, n + 1) + w_perturbed_new(l, m, n))*(vapor_z_initial(n + 1) + vapor_z_initial(n)) - &
                (w_perturbed_new(l, m, n - 1) + w_perturbed_new(l, m, n))*(vapor_z_initial(n - 1) + vapor_z_initial(n)))/4.

      aux = -(u_perturbed_new(l + 1, m, n) - u_perturbed_new(l - 1, m, n) &
              + v_perturbed_new(l, m + 1, n) - v_perturbed_new(l, m - 1, n))* &
            vapor_z_initial(n)/2.

      escal = dqvap(1)*KM1 + dqvap(2)*KM2 + dqvap(3)*KM3

      lapla = ((vapor_base(l + 1, m, n) + vapor_base(l - 1, m, n)) + &
               (vapor_base(l, m + 1, n) + vapor_base(l, m - 1, n))) + &
              vapor_base(l, m, n + 1) + vapor_base(l, m, n - 1) - &
              6.*vapor_base(l, m, n)

      lapla = lapla + (vapor_z_initial(n + 1) + vapor_z_initial(n - 1) - 2.*vapor_z_initial(n))

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      vapor_new(l, m, n) = dt1*((advec + verti + aux)/dx1 + turbul) + vapor_base(l, m, n)

      aux = dt1/dx1
      aux = aux*(verti + advec) + dt1*turbul + vapor_base(l, m, n)

      return
   end subroutine dvapor

   subroutine dgotit(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dgotit_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqgot(1) = drop_base(l + 1, m, n) - drop_base(l - 1, m, n)
      dqgot(2) = drop_base(l, m + 1, n) - drop_base(l, m - 1, n)
      dqgot(3) = drop_base(l, m, n + 1) - drop_base(l, m, n - 1)

      adv(1) = ((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(drop_base(l + 1, m, n) + drop_base(l, m, n)) &
                - (u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(drop_base(l - 1, m, n) + drop_base(l, m, n)))/4.
      adv(1) = adv(1) + dqgot(1)/2.*u_z_initial(n)

      adv(2) = ((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))*(drop_base(l, m + 1, n) + drop_base(l, m, n)) &
                - (v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))*(drop_base(l, m - 1, n) + drop_base(l, m, n)))/4.
      adv(2) = adv(2) + dqgot(2)/2.*v_z_initial(n)

      advgot2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (drop_base(l, m, n) + drop_base(l, m, n + 1))/4.
      if ((advgot2(l, m) - advgot1(l, m))*dt1/dx1 > drop_base(l, m, n) .and. &
          w_perturbed_new(l, m, n) > 0) then
         advgot2(l, m) = advgot1(l, m) + drop_base(l, m, n)*dx1/dt1
      end if
      adv(3) = advgot2(l, m) - advgot1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))/dx1

      escal = dqgot(1)*KM1 + dqgot(2)*KM2 + dqgot(3)*KM3

      lapla = drop_base(l + 1, m, n) + drop_base(l, m + 1, n) + drop_base(l, m, n + 1) + &
              drop_base(l - 1, m, n) + drop_base(l, m - 1, n) + drop_base(l, m, n - 1) - &
              6.*drop_base(l, m, n)

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      drop_new(l, m, n) = dt1*(advec + turbul) + drop_base(l, m, n)

      return
   end subroutine dgotit

   subroutine dlluvi(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dlluvi_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqllu(1) = rain_base(l + 1, m, n) - rain_base(l - 1, m, n)
      dqllu(2) = rain_base(l, m + 1, n) - rain_base(l, m - 1, n)
      dqllu(3) = rain_base(l, m, n + 1) - rain_base(l, m, n - 1)

      adv(1) = ((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(rain_base(l + 1, m, n) + rain_base(l, m, n)) &
                - (u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(rain_base(l - 1, m, n) + rain_base(l, m, n)))/4.
      adv(1) = adv(1) + dqllu(1)/2.*u_z_initial(n)

      adv(2) = ((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))*(rain_base(l, m + 1, n) + rain_base(l, m, n)) - &
                (v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))*(rain_base(l, m - 1, n) + rain_base(l, m, n)))/4.
      adv(2) = adv(2) + dqllu(2)/2.*v_z_initial(n)

      advllu2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (rain_base(l, m, n) + rain_base(l, m, n + 1))/4.

      adv(3) = advllu2(l, m) - advllu1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))

      escal = dqllu(1)*KM1 + dqllu(2)*KM2 + dqllu(3)*KM3

      lapla = rain_base(l + 1, m, n) + rain_base(l, m + 1, n) + rain_base(l, m, n + 1) + &
              rain_base(l - 1, m, n) + rain_base(l, m - 1, n) + rain_base(l, m, n - 1) - &
              6.*rain_base(l, m, n)

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      !***  termino de sedimentacion

      Qllus = (rain_base(l, m, n + 1) + rain_base(l, m, n))/2.
      Qllui = (rain_base(l, m, n - 1) + rain_base(l, m, n))/2.
      Rms = (Qllus/cteqllu)**.25
      Rmm = (rain_base(l, m, n)/cteqllu)**.25
      Rmi = (Qllui/cteqllu)**.25
      Vtllus = (Av(2*n + 1)*Rms**.8 + Av(2*n)*Rmm**.8)/2.
      Vtllui = (Av(2*n - 1)*Rmi**.8 + Av(2*n)*Rmm**.8)/2.
      if (n == 1) then
         Vtllui = Av(2*n)*Rmm**.8
         Qllui = rain_base(l, m, n)
      end if

      sedim = gam4p8/6.*(Qllus*Vtllus - Qllui*Vtllui)

      rain_new(l, m, n) = dt1*((advec + sedim)/dx1 + turbul) + rain_base(l, m, n)
      return
   end subroutine dlluvi

   subroutine dcrist(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dcrist_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqcri(1) = crystal_base(l + 1, m, n) - crystal_base(l - 1, m, n)
      dqcri(2) = crystal_base(l, m + 1, n) - crystal_base(l, m - 1, n)
      dqcri(3) = crystal_base(l, m, n + 1) - crystal_base(l, m, n - 1)

      adv(1) = ((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(crystal_base(l + 1, m, n) + crystal_base(l, m, n)) &
                - (u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(crystal_base(l - 1, m, n) + crystal_base(l, m, n)))/4.
      adv(1) = adv(1) + dqcri(1)/2.*u_z_initial(n)

      adv(2) = ((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))*(crystal_base(l, m + 1, n) + crystal_base(l, m, n)) - &
                (v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))*(crystal_base(l, m - 1, n) + crystal_base(l, m, n)))/4.
      adv(2) = adv(2) + dqcri(2)/2.*v_z_initial(n)

      advcri2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (crystal_base(l, m, n) + crystal_base(l, m, n + 1))/4.

      adv(3) = advcri2(l, m) - advcri1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))/dx1

      escal = dqcri(1)*KM1 + dqcri(2)*KM2 + dqcri(3)*KM3

      lapla = crystal_base(l + 1, m, n) + crystal_base(l, m + 1, n) + crystal_base(l, m, n + 1) + &
              crystal_base(l - 1, m, n) + crystal_base(l, m - 1, n) + crystal_base(l, m, n - 1) - &
              6.*crystal_base(l, m, n)

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      crystal_new(l, m, n) = dt1*(advec + turbul) + crystal_base(l, m, n)

      return
   end subroutine dcrist

   subroutine dnieve(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dnieve_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqnie(1) = snow_base(l + 1, m, n) - snow_base(l - 1, m, n)
      dqnie(2) = snow_base(l, m + 1, n) - snow_base(l, m - 1, n)
      dqnie(3) = snow_base(l, m, n + 1) - snow_base(l, m, n - 1)

      adv(1) = ((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(snow_base(l + 1, m, n) + snow_base(l, m, n)) &
                - (u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(snow_base(l - 1, m, n) + snow_base(l, m, n)))/4.
      adv(1) = adv(1) + dqnie(1)/2.*u_z_initial(n)

      adv(2) = ((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))*(snow_base(l, m + 1, n) + snow_base(l, m, n)) - &
                (v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))*(snow_base(l, m - 1, n) + snow_base(l, m, n)))/4.
      adv(2) = adv(2) + dqnie(2)/2.*v_z_initial(n)

      advnie2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (snow_base(l, m, n) + snow_base(l, m, n + 1))/4.

      adv(3) = advnie2(l, m) - advnie1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))

      escal = dqnie(1)*KM1 + dqnie(2)*KM2 + dqnie(3)*KM3

      lapla = snow_base(l + 1, m, n) + snow_base(l, m + 1, n) + snow_base(l, m, n + 1) + &
              snow_base(l - 1, m, n) + snow_base(l, m - 1, n) + snow_base(l, m, n - 1) - &
              6.*snow_base(l, m, n)

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      !***  termino de sedimentacion

      Qnies = (snow_base(l, m, n + 1) + snow_base(l, m, n))/2.
      Qniei = (snow_base(l, m, n - 1) + snow_base(l, m, n))/2.

      sedim = Vtnie(2*n + 1)*Qnies - Vtnie(2*n - 1)*Qniei
      snow_new(l, m, n) = dt1*((advec + sedim)/dx1 + turbul) + snow_base(l, m, n)

      return
   end subroutine dnieve

   subroutine dgrani(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use dgrani_vars
      implicit none

      integer, intent(in) :: l, m, n

      dqgra(1) = hail_base(l + 1, m, n) - hail_base(l - 1, m, n)
      dqgra(2) = hail_base(l, m + 1, n) - hail_base(l, m - 1, n)
      dqgra(3) = hail_base(l, m, n + 1) - hail_base(l, m, n - 1)

      adv(1) = ((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(hail_base(l + 1, m, n) + hail_base(l, m, n)) &
                - (u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(hail_base(l - 1, m, n) + hail_base(l, m, n)))/4.
      adv(1) = adv(1) + dqgra(1)/2.*u_z_initial(n)

      adv(2) = ((v_perturbed_new(l, m + 1, n) + v_perturbed_new(l, m, n))*(hail_base(l, m + 1, n) + hail_base(l, m, n)) - &
                (v_perturbed_new(l, m - 1, n) + v_perturbed_new(l, m, n))*(hail_base(l, m - 1, n) + hail_base(l, m, n)))/4.
      adv(2) = adv(2) + dqgra(2)/2.*v_z_initial(n)

      advgra2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (hail_base(l, m, n) + hail_base(l, m, n + 1))/4.

      adv(3) = advgra2(l, m) - advgra1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))

      escal = dqgra(1)*KM1 + dqgra(2)*KM2 + dqgra(3)*KM3

      lapla = hail_base(l + 1, m, n) + hail_base(l, m + 1, n) + hail_base(l, m, n + 1) + &
              hail_base(l - 1, m, n) + hail_base(l, m - 1, n) + hail_base(l, m, n - 1) - &
              6.*hail_base(l, m, n)

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      !***  termino de sedimentacion
      Qgras = (hail_base(l, m, n + 1) + hail_base(l, m, n))/2.
      Qgrai = (hail_base(l, m, n - 1) + hail_base(l, m, n))/2.
      Rms = (Qgras/cteqgra)**.25
      Rmm = (hail_base(l, m, n)/cteqgra)**.25
      Rmi = (Qgrai/cteqgra)**.25
      Vtgras = (Vtgra0(2*n + 1)*Rms**.8 + Vtgra0(2*n)*Rmm**.8)/2.
      Vtgrai = (Vtgra0(2*n - 1)*Rmi**.8 + Vtgra0(2*n)*Rmm**.8)/2.
      if (n == 1) then
         Vtgrai = Vtgra0(2*n)*Rmm**.8
         Qgrai = hail_base(l, m, n)
      end if

      sedim = gam4p8/6.*(Qgras*Vtgras - Qgrai*Vtgrai)
      !***

      hail_new(l, m, n) = dt1*((advec + sedim)/dx1 + turbul) + hail_base(l, m, n)

      return
   end subroutine dgrani

   subroutine daeros(l, m, n)
      use cant01
      use dimensions
      use dinamic_var_perturbation
      use microphysics_perturbation
      use constants
      use initial_z_state
      use advecs
      use turbvar
      use daeros_vars
      implicit none

      integer, intent(in) :: l, m, n

      daer(1) = aerosol_base(l + 1, m, n) - aerosol_base(l - 1, m, n)
      daer(2) = aerosol_base(l, m + 1, n) - aerosol_base(l, m - 1, n)
      daer(3) = aerosol_base(l, m, n + 1) - aerosol_base(l, m, n - 1)

      adv(1) = (((u_perturbed_new(l + 1, m, n) + u_perturbed_new(l, m, n))*(aerosol_base(l + 1, m, n) + aerosol_base(l, m, n))) - &
                ((u_perturbed_new(l - 1, m, n) + u_perturbed_new(l, m, n))*(aerosol_base(l - 1, m, n) + aerosol_base(l, m, n))))/4.
      adv(1) = adv(1) + daer(1)/2.*u_z_initial(n)

      adv(2) = (((v_perturbed_new(l, m + 1, n) &
                  + v_perturbed_new(l, m, n))*(aerosol_base(l, m + 1, n) &
                                               + aerosol_base(l, m, n))) &
                - ((v_perturbed_new(l, m - 1, n) &
                    + v_perturbed_new(l, m, n))*(aerosol_base(l, m - 1, n) &
                                                 + aerosol_base(l, m, n))))/4.
      adv(2) = adv(2) + daer(2)/2.*v_z_initial(n)

      advaer2(l, m) = (w_perturbed_new(l, m, n) + w_perturbed_new(l, m, n + 1))* &
                      (aerosol_base(l, m, n) + aerosol_base(l, m, n + 1))/4.

      adv(3) = advaer2(l, m) - advaer1(l, m)

      advec = -(adv(1) + adv(2) + adv(3))

      verti = -((w_perturbed_new(l, m, n + 1) + w_perturbed_new(l, m, n))*(aerosol_z_initial(n + 1) + aerosol_z_initial(n)) - &
                (w_perturbed_new(l, m, n - 1) + w_perturbed_new(l, m, n))*(aerosol_z_initial(n - 1) + aerosol_z_initial(n)))/4.

      aux = -((u_perturbed_new(l + 1, m, n) - u_perturbed_new(l - 1, m, n)) &
              + (v_perturbed_new(l, m + 1, n) - v_perturbed_new(l, m - 1, n)))* &
            aerosol_z_initial(n)/2.

      escal = daer(1)*KM1 + daer(2)*KM2 + daer(3)*KM3

      lapla = aerosol_base(l + 1, m, n) + aerosol_base(l, m + 1, n) + aerosol_base(l, m, n + 1) + &
              aerosol_base(l - 1, m, n) + aerosol_base(l, m - 1, n) + aerosol_base(l, m, n - 1) - &
              6.*aerosol_base(l, m, n)

      lapla = ((aerosol_base(l + 1, m, n) + aerosol_base(l - 1, m, n)) &
               + (aerosol_base(l, m + 1, n) + aerosol_base(l, m - 1, n))) &
              + aerosol_base(l, m, n - 1) + aerosol_base(l, m, n + 1) &
              - 6.*aerosol_base(l, m, n)

      lapla = lapla + (aerosol_z_initial(n + 1) + aerosol_z_initial(n - 1) - 2.*aerosol_z_initial(n))

      turbul = cteturb*(escal/dx8 + KMM/dx2*lapla)

      aerosol_new(l, m, n) = dt1*((advec + verti + aux)/dx1 + turbul) + aerosol_base(l, m, n)
      return
   end subroutine daeros

   subroutine corgot
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de gotitas
      use dimensions
      use microphysics_perturbation
      use lmngot
      use corgot_vars
      implicit none

      neg1 = 0.
      pos1 = 0.
      do concurrent(n=ngot(1):ngot(2), l=lgot(1):lgot(2), m=mgot(1):mgot(2))
         if (drop_new(l, m, n) < 0.) then
            neg1 = neg1 + drop_new(l, m, n)
            drop_new(l, m, n) = 0
         else
            pos1 = pos1 + drop_new(l, m, n)
         end if
      end do

      if (pos1 <= -neg1) then
         do concurrent(l=lgot(1):lgot(2), m=mgot(1):mgot(2), n=ngot(1):ngot(2))
            drop_new(l, m, n) = 0.
         end do
      else
         aux1 = neg1/pos1
         do concurrent(l=lgot(1):lgot(2), m=mgot(1):mgot(2), n=ngot(1):ngot(2))
            if (drop_new(l, m, n) > 0) then
               drop_new(l, m, n) = drop_new(l, m, n)*(1.+aux1)
            end if
         end do
      end if

      return
   end subroutine corgot

   subroutine corllu
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de gotas
      use dimensions
      use microphysics_perturbation
      use lmnllu
      use corgot_vars
      implicit none

      neg1 = 0.
      pos1 = 0.
      do concurrent(n=nllu(1):nllu(2), l=lllu(1):lllu(2), m=mllu(1):mllu(2))
         if (rain_new(l, m, n) < 0.) then
            neg1 = neg1 + rain_new(l, m, n)
            rain_new(l, m, n) = 0
         else
            pos1 = pos1 + rain_new(l, m, n)
         end if
      end do

      if (pos1 <= -neg1) then
         do concurrent(l=lllu(1):lllu(2), m=mllu(1):mllu(2), n=nllu(1):nllu(2))
            rain_new(l, m, n) = 0.
         end do
      else
         aux1 = neg1/pos1
         do concurrent(l=lllu(1):lllu(2), m=mllu(1):mllu(2), n=nllu(1):nllu(2))
            if (rain_new(l, m, n) > 0) then
               rain_new(l, m, n) = rain_new(l, m, n)*(1.+aux1)
            end if
         end do
      end if
      return
   end subroutine corllu

   subroutine corcri
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de cristales
      use dimensions
      use microphysics_perturbation
      use lmncri
      use corgot_vars
      implicit none
      neg1 = 0.
      pos1 = 0.
      do concurrent(n=ncri(1):ncri(2), l=lcri(1):lcri(2), m=mcri(1):mcri(2))
         if (crystal_new(l, m, n) < 0.) then
            neg1 = neg1 + crystal_new(l, m, n)
            crystal_new(l, m, n) = 0
         else
            pos1 = pos1 + crystal_new(l, m, n)
         end if
      end do

      if (pos1 <= -neg1) then
         do concurrent(l=lcri(1):lcri(2), m=mcri(1):mcri(2), n=ncri(1):ncri(2))
            crystal_new(l, m, n) = 0.
         end do
      else
         aux1 = neg1/pos1
         do concurrent(l=lcri(1):lcri(2), m=mcri(1):mcri(2), n=ncri(1):ncri(2))
            if (crystal_new(l, m, n) > 0) then
               crystal_new(l, m, n) = crystal_new(l, m, n)*(1.+aux1)
            end if
         end do
      end if

      return
   end subroutine corcri

   subroutine cornie
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de nieve
      use dimensions
      use microphysics_perturbation
      use lmnnie
      use corgot_vars
      implicit none
      neg1 = 0.
      pos1 = 0.
      do concurrent(n=nnie(1):nnie(2), l=lnie(1):lnie(2), m=mnie(1):mnie(2))
         if (snow_new(l, m, n) < 0.) then
            neg1 = neg1 + snow_new(l, m, n)
            snow_new(l, m, n) = 0
         else
            pos1 = pos1 + snow_new(l, m, n)
         end if
      end do

      if (pos1 <= -neg1) then
         do concurrent(l=lnie(1):lnie(2), m=mnie(1):mnie(2), n=nnie(1):nnie(2))
            snow_new(l, m, n) = 0.
         end do
      else
         aux1 = neg1/pos1
         do concurrent(l=lnie(1):lnie(2), m=mnie(1):mnie(2), n=nnie(1):nnie(2))
            if (snow_new(l, m, n) > 0) then
               snow_new(l, m, n) = snow_new(l, m, n)*(1.+aux1)
            end if
         end do
      end if
      return
   end subroutine cornie

   subroutine corgra
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de granizos
      use dimensions
      use microphysics_perturbation
      use lmngra
      use corgot_vars
      implicit none
      neg1 = 0.
      pos1 = 0.
      do concurrent(n=ngra(1):ngra(2), l=lgra(1):lgra(2), m=mgra(1):mgra(2))
         if (hail_new(l, m, n) < 0.) then
            neg1 = neg1 + hail_new(l, m, n)
            hail_new(l, m, n) = 0
         else
            pos1 = pos1 + hail_new(l, m, n)
         end if
      end do

      if (pos1 <= -neg1) then
         do concurrent(l=lgra(1):lgra(2), m=mgra(1):mgra(2), n=ngra(1):ngra(2))
            hail_new(l, m, n) = 0.
         end do
      else
         aux1 = neg1/pos1
         do concurrent(l=lgra(1):lgra(2), m=mgra(1):mgra(2), n=ngra(1):ngra(2))
            if (hail_new(l, m, n) > 0) then
               hail_new(l, m, n) = hail_new(l, m, n)*(1.+aux1)
            end if
         end do
      end if
      return
   end subroutine corgra

   subroutine corvap(Qvapneg)
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de vapor
      use dimensions
      use microphysics_perturbation
      use initial_z_state
      use corvap_vars
      implicit none

      real(8), intent(in) :: Qvapneg

      do concurrent(k=1:nz1)
         dq = Qvapneg*vapor_z_relative(k)/nx1**2.
         do concurrent(i=1:nx1, j=1:nx1)
            vapor_new(i, j, k) = vapor_new(i, j, k) + dq
            if (vapor_new(i, j, k) + vapor_z_initial(k) < 0) vapor_new(i, j, k) = -vapor_z_initial(k)
         end do
      end do

      return
   end subroutine corvap

   subroutine coraer(aerneg)
      !! Esta subrutina corrige los lugares en donde la dinamica da
      !! negativa la cantidad de aerosoles
      use dimensions
      use microphysics_perturbation
      use initial_z_state
      use coraer_vars
      implicit none

      real(8), intent(in) :: aerneg

      do concurrent(k=1:nz1)
         dq = aerneg*aerosol_z_relative(k)/nx1**2.
         do concurrent(i=1:nx1, j=1:nx1)
            aerosol_new(i, j, k) = aerosol_new(i, j, k) + dq
            if (aerosol_new(i, j, k) + aerosol_z_initial(k) < 0) aerosol_new(i, j, k) = -aerosol_z_initial(k)
         end do
      end do
      return
   end subroutine coraer

   subroutine suma(sum, a1, a2, a3)
      implicit none
      real a1, a2, a3, sum, aux
      integer j

      do concurrent(j=1:2)
         if (a1 > a2) then
            aux = a1
            a1 = a2
            a2 = aux
         end if
         if (a2 > a3) then
            aux = a2
            a2 = a3
            a3 = aux
         end if
      end do

      sum = (a1 + a2) + a3
      return
   end subroutine suma

   subroutine nuclea(Qvap, Qliq, Naer, TT, rhoa, e1, esl, ess, rl, rs, Lvl, Lvs, Naux, auxl, auxs)
      use cant01
      use dimensions
      use constants
      use nuclea61
      implicit none

      real, intent(in) ::  Naer, rhoa, rs, Lvl, Lvs
      real, intent(inout) :: Qliq, Qvap, TT, e1, esl, ess, rl
      real, intent(inout) :: Naux, auxl, auxs   !     Numero de aesosoles

      B = Lvl/Rv
      auxl = 0.
      auxs = 0.
      Naux = 0.

      Rcri = 5e-5
      if (TT < T0) Rcri = Rcri - 4e-5*(T0 - TT)/40.
      if (Rcri < 1e-5) Rcri = 1e-5

      !     nucleacion sobre cristales
      Tc = T0 - TT
      if (Tc > 0 .and. rs > 0 .and. Naer > 0) then

         mcri = pi*Rcri**3./10.*rhocri

         Naux = Acri*exp(Bcri*Tc)
         if (Naux > .9) Naux = .9
         Naux = Naux*Naer*1e6   ! en m3
         auxs = Naux*mcri

         if (auxs > (Qvap - ess/Rv/TT*.95)) then
            auxs = (Qvap - ess/Rv/TT)*.95
            Naux = auxs/mcri     ! en m3
         end if
         Qvap = Qvap - auxs
         TT1 = TT
         TT2 = TT + (auxs*Lvs)/(Cp*rhoa)
         TT = TT2
         e1 = Qvap*Rv*TT
         esl = esl*exp(B*(TT2 - TT1)/TT2/TT1)
         ess = ess*exp(B*(TT2 - TT1)/TT2/TT1)
         Naux = -Naux/1e6        ! en cm3

      end if

      !     nucleacion sobre gotitas
      if (e1 > esl) then
         s = 0
         hhh = 0
         xxx = 0
         TT1 = TT
         Ti = TT
         ei = e1
         esli = esl
         auxl = 0.
         caux = B*esli/Ti
         !TODO: Check this loop
10       continue
         F0 = Lvl/Rv*(esl/TT1 - ei/Ti) + Cp*rhoa*(TT1 - Ti)
         F0p = Cp*rhoa + B/TT1**2.*caux
         TT2 = TT1 - F0/F0p
         auxl = (ei/Ti - esl/TT2)/Rv
         Qliq1 = Qliq + auxl
         e1 = esl

         if (Qliq1 < 0) then
            Qliq1 = 0.
            auxl = -Qliq
            TT2 = esl/(e1/TT1 - auxl*Rv)
            e1 = (Qvap - auxl)*Rv*TT2
            hhh = 1
         end if

         s = 1
         esl = esl*exp(B*(TT2 - TT1)/TT2/TT1)
         TT1 = TT2
         rl = abs(e1 - esl)/esl
         xxx = xxx + 1
         if (rl > 1e-3 .and. hhh == 0) goto 10

         Qvap = Qvap - auxl
         Qliq = Qliq + auxl
         !*
         !      variacion en los aerosoles
         !      considerando que las nuevas gotitas tienen un radio Rgotmin

         Naux = Naux - auxl/(4./3.*pi*rhow*Rgotmin**3.)/1e6

         TT = TT2

      end if

      return
   end subroutine nuclea
end module extra_subrut
