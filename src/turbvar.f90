!      include para las variables usadas en turbulencia
module turbvar
      real :: KMM,KM1,KM2,KM3
      real, dimension(3,3) :: DD
      real, dimension(3) :: D1,D2,D3
      common /turbcom/ KMM,KM1,KM2,KM3,DD,D1,D2,D3
end module turbvar