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

module lmnllu
!> Posiciones en las cuales Qllu<0
   integer, dimension(2) :: lllu,mllu,nllu
end module lmnllu

module lmnnie
!> Posiciones en las cuales Qnie<0
   integer, dimension(2) :: lnie,mnie,nnie
end module lmnnie

module coraer_vars
!> Correccion de vapor
   integer i,j,k
   real dq
end module coraer_vars

module corgot_vars
   integer l,m,n
   real aux1,pos1,neg1
end module corgot_vars

module corvap_vars
!> Correccion de vapor
   integer i,j,k
   real dq
end module corvap_vars

module turbu1_vars
!> Turbu
   real dv(3,3)
   integer lx,ly,lz,ldis
   integer i,j,k,n,m
   real vel(3,-5:5,-5:5,-5:5)
end module turbu1_vars

module turbu2_vars
!> Turbu
   real sum
   real KM(-3:3,-3:3,-3:3)
   integer lx,ly,lz
   integer n,m,ldis
end module turbu2_vars

module tempe01
!> Tempe

   real dtita(3)

   real adv(3)
   real advec,verti,escal,lapla,turden,turbul,calor
end module tempe01

module nuclea61
!> Nuclea
   real Qliq1,TT1,TT2
   real B,Tc
   real Ti,ei,esli
   real F0,F0p
   real mcri,Rcri,caux
   integer hhh,s,xxx

!> Parametros de las particulas
   real Rgotmin,Acri,Bcri
   parameter (Rgotmin=5e-6,Acri=1e-11,Bcri=.6)   !A en cm^-3
end module nuclea61


module p3v3
!> Velocidades y las presiones
   USE dimen
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: U3
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: V3
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: W3
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Pres3
end module p3v3
