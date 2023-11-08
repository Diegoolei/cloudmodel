module const
!> Constantes matematicas y fisicas para precipitacion
   real, parameter :: pi = 3.1415926, gam1p8 = .9134, gam2p8 = 1.6765, gam3p8 = 4.6941742, &
      gam4p8 = 17.837862, G = 9.8, Rd = 287.04, Rv = 461.05, Kapa = .2857, Eps = .622646, &
      T0 = 273.15, P00 = 101300, Kair = 2.40e-2, Cwl = 4218., &
      Cwv = 1839., Cwi = 2106., Lvl0 = 2.500e6, Lsl0 = 79.7, &
      Lvs0 = 2.835e6, Cp = 1003., Cv = 716.,&
      elvs0 = 610.78, esvs0 = 610.918,&
      Dv0 = 2.11e-5, Vis0 = 1.718e-5, rhow = 1000., rhocri = 900.&
      , rhonie = 100., rhogra = 500., N0got = 2.9e24, N0llu = 4912189.&
      , N0nie = 1.66e5, N0gra = 310.&
      , Av0 = 1455., Vtnie0 = .5, Efcol = .8, Efcolgn = .7
   real, dimension(210:320) :: Tvis, Tlvl, Tlsl, Tlvs, Telvs, Tesvs, Eautcn, Eacrcn
end module const

module dimen
!> Numero de puntos y los intervalos temporales y espaciales
   integer, parameter :: nx1 = 50, nz1 = 45, nz = 64
   real, parameter :: dx1 = 300., dt1 = 2., dt2 = 1., dt3 = .2
end module dimen

