*     include para las velocidades
      integer t
      real dvx,dvy,dvz,diver
      real dprex,dprey,dprez
      real vel0,vel1,vel2,vel3
      real presi,presix,presiy,presiz
      real facx,facy,facz
      integer i,j,k
      real prom1,prom,kkk,presprom,nnn

      facx=.05
      facy=.05
      facz=.05

      prom=.3/6.*(dt3/.2)
      prom1=1.-prom*6.
      kkk=.01
      nnn=(nx1+2)**2.*(nz1+1)
