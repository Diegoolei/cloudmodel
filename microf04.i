*     include para la microfisica

*     cantidades 
      real els,ess,Lvl,Lvs,Lsl,T,Fcal,Dv,nu,daer2,Naer
      real Lsl00
      real Eaucn,Eaccn,Eacng

*     variables auxiliares
      integer l,m,n,yy
      real*8 qauxl,qauxs,aux,qauxl0
*$$
      real*8 qlluaux,qgotaux,qvapaux,qcriaux,qnieaux,qgraaux,Naux

*     variables de procesos microfisicos
*$$
      real*8 coevgot,coevllu,coevcri,coevnie,coevgra
      real*8 liqconv,hieconv,nieconv,acgollu
      real*8 acgonie,accrnie,acgogra,acllgra,accrgra,acnigra
      real*8 cfgotcri,cfllunie,cfllugra,fugrallu,libaer
      real*8 coaergot,coaerllu,coaercri,coaernie,coaergra
      real*8 ccnigra,colilc,coliln,congagua
      real*8 mucrgrni,mucrgrgr
      real*8 invapgot,invapllu,invapcri,invapnie,invapgra
      real*8 ingotllu,ingotcri,ingotnie,ingotgra
      real*8 incrinie,incrigra,inllunie,inllugra,inniegra
      real*8 Intvap,Intgot,Intllu,Intcri,Intnie,Intgra,Intaer
      real*8 cfln1,cfln2

*     numero de cristales por colision
      real Ncrgrni,Ncrgrgr
      parameter (Ncrgrni=20.,Ncrgrgr=20.)
            
*     variables para la evaporacion y la condensacion
      real Qvls,Qvss,Qvls0
      
*     parametros de las particulas
*$$
      real*8 Rgot,Ngot,Rllu,Nllu,Rcri,Ncri,Rnie,Nnie,Rgra,Ngra
      real*8 Vtm,Vtgra,Nre,Nsc,A,fventl,fventn,fventgs,fventgl,Rgotmin
      parameter (Rgotmin=7e-6)
      

      integer s
*$$
      real*8 qvapaux1,qgotaux1,qlluaux1,qcriaux1,qnieaux1,qgraaux1

*     variables para tgra
      real Tg,agual,hielo,alfagra,fugra,Aalfa,Balfa
      real A1,A2,A3,A3b,A4,B1,B2b,B3,B4,BB,CC2,CC3
      real Q1,Q2,Q3,Q4,Qt,dQt,Qvaux,esaux,Taux
      real Fcalgra
      integer crecigra,i

*     parametros para los cristales
      real Acri,Bcri
      parameter (Acri=1e-2,Bcri=.6)  !A en m^-3

*     difusion de aerosoles
      real Dfaer,Efcaer
      parameter (Dfaer=1e-10,Efcaer=.01)