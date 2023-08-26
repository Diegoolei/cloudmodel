*     Include para las constantes matematicas y fisicas
*                   para precipitacion (7/6/99)

      real pi,gam1p8,gam2p8,gam3p8,gam4p8,G
      parameter (pi=3.1415926,gam1p8=.9134,gam2p8=1.6765,
     &           gam3p8=4.6941742,gam4p8=17.837862,G=9.8)

      real Rd,Rv,Kapa,Eps,T0,P00,Kair,Cwl,Cwv,Cwi,Lvl0,Lsl0,Lvs0,
     &     Cp,Cv
      parameter (Rd=287.04,Rv=461.05,Kapa=.2857,Eps=.622646,
     &           T0=273.15,P00=101300,Kair=2.40e-2,Cwl=4218.,
     &           Cwv=1839.,Cwi=2106.,Lvl0=2.500e6,Lsl0=79.7,
     &           Lvs0=2.835e6,Cp=1003.,Cv=716.)

      real elvs0,esvs0
      parameter(elvs0=610.78,esvs0=610.918)
*$$
      real Dv0,Vis0,rhow,rhocri,rhonie,rhogra,N0got,N0llu,N0nie
     &     ,N0gra,Av0,Vtnie0,Efcol,Efcolgn

      parameter (Dv0=2.11e-5,Vis0=1.718e-5,rhow=1000.,rhocri=900.
     &          ,rhonie=100.,rhogra=500.,N0got=2.9e24,N0llu=4912189.
     &          ,N0nie=1.66e5,N0gra=310.  !(mod 20/11/99)
*     &          ,N0nie=1.66e5,N0gra=1000.  !(mod 7/6/99)
     &          ,Av0=1455.,Vtnie0=.5,Efcol=.8,Efcolgn=.7)  !(mod 7/6/99)

      dimension Tvis(210:320)
      dimension Tlvl(210:320)
      dimension Tlsl(210:320)
      dimension Tlvs(210:320)
      dimension Telvs(210:320)
      dimension Tesvs(210:320)
      dimension Eautcn(210:320)
      dimension Eacrcn(210:320)
      
      common /const/Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Eautcn,Eacrcn
      real Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Eautcn,Eacrcn








