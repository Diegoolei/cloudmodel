!     Incluye microfisica con vapor, gotitas, lluvia, cristales, nieve
!      y granizos (por unidad de volumen)
!     Con viento de corte.
!     Este modelo simula una nube tridimensional, con diferencias
!      finitas adelantadas en el tiempo y centradas en el espacio
!     Este modelo (asi como sus variantes) sirve de test
!     Las grillas son las mismas para las cantidades dinamicas y
!      microfisicas, como asi tambien los intervalos de tiempo.
!     Todas las variables son reales*4
!     Condiciones de contorno homogeneas para las variables microfisicas
!     Graba el valor de todas las variables cada ltb segundos.
!     Graba el valor de las variables para analisis cada ltg segundos.
!     Condicion de contorno nula para el vapor
!     Contempla el desplazamiento de la nube.
!     Mejora la condicion en el piso para los aerosoles cuando hay agua (815)

!> @brief This module contains the cloud model implementation.
module cloud_model
contains
   subroutine model()
      USE cant01, only: lt1
      USE model_var, only: tt
      USE perdim, only: Titaa1
      USE model_initialization, only: initialize_model
      USE model_aux, only: vapor_advection, dinamics, negative_correction, water_calculation,&
         microphisics_substring, floor_and_ceiling_contour, lateral_contour,&
         floor_condition_redefinition, floor_and_ceiling_contour_redefinition,&
         lateral_contour_redefinition, vapour_negative_correction, save_backup
      implicit none
      
      call initialize_model()
      do tt=1,lt1
         call vapor_advection()
         call dinamics()
         call negative_correction()
         call water_calculation()
         call microphisics_substring()
         call floor_and_ceiling_contour()
         call lateral_contour()
         call speed_pressure()
         call floor_condition_redefinition()
         call lateral_contour_redefinition()
         call filtro(Titaa1,.01,.01,.02)
         call vapour_negative_correction()
         call save_backup()
         write(*,*) '----Tiempo transcurrido:',tt,'de',lt1,'----'
      end do

   end subroutine model
end module cloud_model
