!     Include para las cantidades no perturbadas
module estbas
      USE dimen
      real, dimension(-3:nz1+3) :: Temp0,Tita0,Pres00,Presi0,UU,VV,cc2,Den0,aer0,Qvap0
      real, dimension(nz1) :: Qvaprel,aerrel
   
      common /estbase/ Temp0,Tita0,Pres00,Presi0,UU,VV,cc2,Den0,aer0,Qvap0,Qvaprel,aerrel
   end module estbas