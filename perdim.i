*     Include para las perturbaciones de las variables dinamicas

      dimension U1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension V1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension W1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension U2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension V2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension W2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      common /UVW/U1,V1,W1,U2,V2,W2
      real U2,V2,W2,U1,V1,W1
      
      dimension Pres1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Pres2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      common /pre/ Pres1,Pres2
      real Pres1,Pres2

      dimension Tempa1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Titaa1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Titaa2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      common /temps/ Tempa1,Titaa1,Titaa2
      real Tempa1,Titaa1,Titaa2

      dimension Ktur(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Fcalo(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      common /caltur/ Ktur,Fcalo
      real Ktur,Fcalo
