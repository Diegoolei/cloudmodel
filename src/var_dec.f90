
module turbvar1
!> Variables usadas en turbulencia
   USE dimen
   real :: D(3, 3, -3:nx1+2, -3:nx1 + 2, 3)
end module turbvar1

module cant01
!> Cantidades auxiliares, pasos de tiempos, distancias, etc.
   implicit none
   real :: ltt, ltg, lte, ltb
   integer :: ini, lt1, lt2, lt3
   real :: ctur, cteturb
   real :: dx2, dx8, dx12
   real :: AA, ikapa
   real :: pro1, pro2, pro3, pro4
   real :: cteqgot, cteqllu, cteqnie, cteqgra
end module cant01

module turbvar
!> Variables usadas en turbulencia
   real :: KMM, KM1, KM2, KM3
   real, dimension(3,3) :: DD
   real, dimension(3) :: D1, D2, D3
end module turbvar

module advecs
!> Terminos de adveccion
   USE dimen
   real, dimension(-2:nx1+2,-2:nx1+2) :: advaer1, advaer2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advgot1, advgot2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advllu1, advllu2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advcri1, advcri2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advnie1, advnie2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advgra1, advgra2
   real, dimension(-2:nx1+2,-2:nx1+2) :: advvap1, advvap2
end module advecs

module estbas
!> Cantidades no perturbadas
   USE dimen
   real, dimension(-3:nz1+3) :: Temp0,Tita0,Pres00,Presi0,UU,VV,cc2,Den0,aer0,Qvap0
   real, dimension(nz1) :: Qvaprel,aerrel
end module estbas

module lmncri
!> Posiciones en las cuales Qcri<0
   integer, dimension(2) :: lcri,mcri,ncri
end module lmncri

module lmngot
   !> Posiciones en las cuales Qgot<0
   integer, dimension(2) :: lgot,mgot,ngot
end module lmngot

module lmngra
   !> Posiciones en las cuales Qgra<0
   integer, dimension(2) :: lgra,mgra,ngra
end module lmngra

