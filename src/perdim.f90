!     Include para las perturbaciones de las variables dinamicas

module perdim
      USE dimen
      real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: U1,V1,W1,U2,V2,W2
      real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Pres1,Pres2
      real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Tempa1,Titaa1,Titaa2
      real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Ktur,Fcalo

      common /UVW/U1,V1,W1,U2,V2,W2
      common /pre/ Pres1,Pres2
      common /temps/ Tempa1,Titaa1,Titaa2
      common /caltur/ Ktur,Fcalo
end module perdim