*     include de las velocidades y las presiones

      dimension U3(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension V3(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension W3(-3:nx1+3,-3:nx1+3,-2:nz1+2)      
      dimension Pres3(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      common /vp3/ Pres3,U3,V3,W3
      real Pres3,U3,V3,W3
