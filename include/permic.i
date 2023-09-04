!     Include para las perturbaciones de las variables microfisicas

      dimension Qvap1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qvap2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qgot1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qgot2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension aer1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension aer2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qllu1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qllu2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qcri1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qcri2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qnie1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qnie2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qgra1(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Qgra2(-3:nx1+3,-3:nx1+3,-2:nz1+2)
      dimension Av(-3:2*nz1+5),Vtnie(-3:2*nz1+5),Vtgra0(-3:2*nz1+5)
      common /microf/Qvap1,Qvap2,Qgot1,Qgot2,aer1,aer2,Qllu1,Qllu2,
     &       Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2
      real Qvap1,Qvap2,Qgot1,Qgot2,aer1,aer2,Qllu1,Qllu2,
     &     Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2

      common /velter/Av,Vtnie,Vtgra0
      real Av,Vtnie,Vtgra0
