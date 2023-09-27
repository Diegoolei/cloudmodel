module permic
      USE dimen
      real, dimension(-3:nx1+3,-3:nx1+3,-2:nz1+2) :: Qvap1,Qvap2,Qgot1,Qgot2,aer1,aer2,Qllu1,Qllu2,Qcri1,Qcri2,&
      Qnie1,Qnie2,Qgra1,Qgra2
      real, dimension(-3:2*nz1+5) :: Av,Vtnie,Vtgra0

      common /microf/Qvap1,Qvap2,Qgot1,Qgot2,aer1,aer2,Qllu1,Qllu2,Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2
      common /velter/Av,Vtnie,Vtgra0
end module permic