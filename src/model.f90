
module cloud_model
   !! This module contains the cloud model implementation.
contains
   subroutine model()
      !! Incluye microfisica con vapor, gotitas, lluvia, cristales, nieve
      !! y granizos (por unidad de volumen)
      !! Con viento de corte.
      !! Este modelo simula una nube tridimensional, con diferencias
      !! finitas adelantadas en el tiempo y centradas en el espacio
      !! Este modelo (asi como sus variantes) sirve de test
      !! Las grillas son las mismas para las cantidades dinamicas y
      !! microfisicas, como asi tambien los intervalos de tiempo.
      !! Todas las variables son reales*4
      !! Condiciones de contorno homogeneas para las variables microfisicas
      !! Graba el valor de todas las variables cada ltb segundos.
      !! Graba el valor de las variables para analisis cada ltg segundos.
      !! Condicion de contorno nula para el vapor
      !! Contempla el desplazamiento de la nube.
      !! Mejora la condicion en el piso para los aerosoles cuando hay agua (815)

      use cant01, only: total_time
      use model_var, only: current_time
      use dinamic_var_perturbation, only: theta_base
      use model_initialization, only: initialize_model
      use model_aux, only: vapor_advection, dinamics, negative_correction, water_calculation, &
                           microphisics_substring, floor_and_ceiling_contour, lateral_contour, &
                           floor_condition_redefinition, floor_and_ceiling_contour_redefinition, &
                           lateral_contour_redefinition, vapour_negative_correction, save_backup, &
                           speed_pressure, filtro
      use memory_managment, only: allocate_model, deallocate_model
      use, intrinsic :: iso_fortran_env, only: I4P => int32, R8P => real64
      use, intrinsic :: iso_fortran_env
      use forbear, only: bar_object
      use get_cut, only: get_cutted
      implicit none
      type(bar_object) :: progress_bar
      real(R8P)        :: progress_percent
      call allocate_model()
      call initialize_model()
      call get_cutted()
      !call progress_bar%initialize(filled_char_string='㊂', empty_char_string='●', &
      !                             suffix_string='| ', add_progress_percent=.true., prefix_string='Progress |', &
      !                             scale_bar_color_fg='blue', scale_bar_style='underline_on', spinner_string='(  ●   )')
      !call progress_bar%start
      !do current_time = 1, total_time
      !   call vapor_advection()
      !   call dinamics()
      !   call negative_correction()
      !   call water_calculation()
      !   call microphisics_substring()
      !   call floor_and_ceiling_contour()
      !   call lateral_contour()
      !  call speed_pressure()
      !   call floor_condition_redefinition()
      !   call floor_and_ceiling_contour_redefinition()
      !   call lateral_contour_redefinition()
      !   call filtro(theta_base, .01, .01, .02)
      !   call vapour_negative_correction()
      !   call save_backup()
      !   progress_percent = real(current_time, R8P)/real(total_time, R8P)
      !   call progress_bar%update(current=progress_percent)
      !end do

      !call progress_bar%destroy
      call deallocate_model()
   end subroutine model
end module cloud_model
