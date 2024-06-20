module turbvar1
!> Variables usadas en turbulencia
   USE dimen
   real :: D(3, 3, -3:nx1+2, -3:nx1 + 2, 3)
end module turbvar1

module cant01
!> Cantidades auxiliares, pasos de tiempos, distancias, etc.
   implicit none
   real :: ltt, ltg, lte, ltb, ctur, cteturb, dx2, dx8, dx12, AA, ikapa, pro1,&
      pro2, pro3, pro4, cteqgot, cteqllu, cteqnie, cteqgra
   integer :: ini, lt1, lt2, lt3
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
   real, dimension(-2:nx1+2,-2:nx1+2) :: advaer1, advaer2, advgot1, advgot2, advllu1,&
      advllu2, advcri1, advcri2, advnie1, advnie2, advgra1, advgra2, advvap1, advvap2

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
   real vel(3,-5:5,-5:5,-5:5)
   integer lx, ly, lz, ldis, i, j, k, n, m
end module turbu1_vars

module turbu2_vars
!> Turbu
   real sum
   real KM(-3:3,-3:3,-3:3)
   integer lx, ly, lz, n, m, ldis
end module turbu2_vars

module tempe01
!> Tempe
   real dtita(3),adv(3)
   real advec,verti,escal,lapla,turden,turbul,calor
end module tempe01

module nuclea61
!> Nuclea
   real Qliq1,TT1,TT2,B,Tc,Ti,ei,esli,F0,F0p,mcri,Rcri,caux
   integer hhh,s,xxx

!> Parametros de las particulas
   real Rgotmin,Acri,Bcri
   parameter (Rgotmin=5e-6,Acri=1e-11,Bcri=.6)   !A en cm^-3
end module nuclea61


module p3v3
!> Velocidades y las presiones
   USE dimen
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: U3, V3, W3, Pres3
end module p3v3

module daeros_vars
   real daer(3), adv(3)
   real advec,verti,escal,lapla,turbul, aux
end module daeros_vars

module dcrist_vars
   real dqcri(3), adv(3)
   real advec,escal,lapla,turbul
end module dcrist_vars

module dgotit_vars
   real dqgot(3), adv(3)
   real advec,escal,lapla,turbul
end module dgotit_vars

module dgrani_vars
   real dqgra(3), adv(3)
   real advec,escal,lapla,turbul,sedim, Qgras,Qgrai,Rms,Rmm,Rmi,Vtgras,Vtgrai
end module dgrani_vars

module dimlee_vars
   integer nx2,ny2,nz2
   parameter(nx2=32,ny2=32,nz2=32)
end module dimlee_vars

module dlluvi_vars
   real dqllu(3), adv(3)
   real advec,escal,lapla,turbul,sedim, Qllus,Qllui,Rms,Rmm,Rmi,Vtllus,Vtllui
end module dlluvi_vars

module dnieve_vars
   real dqnie(3), adv(3)
   real advec,escal,lapla,turbul,sedim, Qnies,Qniei
end module dnieve_vars

module dvapor_vars
   real dqvap(3), adv(3)
   real advec,verti,escal,lapla,turbul, aux
end module dvapor_vars

module filtro01
!> filtro01
   USE dimen
   real varia2(-5:nx1+5,-5:nx1+5,-4:nz1+4)
   real varx,vary,varz, fact
   integer i,j,k
end module filtro01

module fuvw
   USE dimen
   real, dimension(-1:nx1+2,-1:nx1+2,-1:nz1+2) :: fu, fv, fw, fp
end module fuvw

module inomo_var
!> Inomo
   real(8) dvelxx,dvelxy,dvelxz,dvelyx,dvelyy,dvelyz,dvelzx,dvelzy,dvelzz,&
      a1,a2,a3, turbulx,turbuly,turbulz
   real diverx,divery,diverz,grave,laplap
end module inomo_var

module permic
   USE dimen
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Qvap1,Qvap2,Qgot1,Qgot2,aer1,&
      aer2,Qllu1,Qllu2,Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2
   real, dimension(-3:2*nz1+5) :: Av,Vtnie,Vtgra0
end module permic

module perdim
!> Perturbaciones de las variables dinamicas
   USE dimen
   real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: U1,V1,W1,U2,V2,W2,Pres1,Pres2,&
      Tempa1,Titaa1,Titaa2,Fcalo
end module perdim

!> @file var_dec.f90
!! @brief This file contains the module declaration for the microf05 module.
!! @details The microf05 module is responsible for declaring variables used in the Fortran 77 cloud model.
!!          This module is part of the Fortran 77 Cloud Model project.
module microf05
!     variables auxiliares
   real(8) qauxl,qauxs,aux,qauxl0,Naux
!     variables de procesos microfisicos
   real(8) coevgot,coevllu,coevcri,coevnie,coevgra, liqconv,hieconv,nieconv,acgollu,&
      acgonie,accrnie,acgogra,acllgra,accrgra,acnigra, cfgotcri,cfllunie,cfllugra, &
      fugrallu,libaer, coaergot,coaerllu,coaercri,coaernie,coaergra, ccnigra,colilc,&
      coliln,congagua, mucrgrni,mucrgrgr, invapgot,invapllu,invapcri,invapnie,invapgra,&
      ingotllu,ingotcri,ingotnie,ingotgra, incrinie,incrigra,inllunie,inllugra,inniegra,&
      Intvap,Intgot,Intllu,Intcri,Intnie,Intgra,Intaer, cfln1,cfln2
!     numero de cristales por colision
   real Ncrgrni,Ncrgrgr
   parameter (Ncrgrni=20.,Ncrgrgr=20.)
!     variables para la evaporacion y la condensacion
   real Qvls,Qvss,Qvls0
!     parametros de las particulas
   real(8) Rgot,Ngot,Rllu,Nllu,Rcri,Ncri,Rnie,Nnie,Rgra,Ngra, Vtm,Vtgra,Nre,Nsc,A,fventl,&
      fventn,fventgs,fventgl,Rgotmin, qvapaux1,qgotaux1,qlluaux1,qcriaux1,qnieaux1,qgraaux1
   parameter (Rgotmin=5e-6)
   integer s
!     variables para tgra
   real Tg,agual,hielo,alfagra,fugra,Aalfa,Balfa, A1,A2,A3,A3b,A4,B1,B2b,B3,B4,BB,CC2,CC3,&
      Q1,Q2,Q3,Q4,Qt,dQt,Qvaux,esaux,Taux, Fcalgra
   integer crecigra,i
!     parametros para los cristales
   real Acri,Bcri
   parameter (Acri=1e-2,Bcri=.6)  !A en m^-3
!     difusion de aerosoles
   real Dfaer,Efcaer
   parameter (Dfaer=1e-10,Efcaer=.01)
end module microf05
