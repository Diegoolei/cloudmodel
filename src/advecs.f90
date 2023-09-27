!     Include para los terminos de adveccion 29/12/98
module advecs
      USE dimen
      real, dimension(-2:nx1+2,-2:nx1+2) :: advaer1,advaer2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advgot1,advgot2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advllu1,advllu2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advcri1,advcri2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advnie1,advnie2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advgra1,advgra2
      real, dimension(-2:nx1+2,-2:nx1+2) :: advvap1,advvap2
      common /advaer/ advaer1,advaer2
      common /advgot/ advgot1,advgot2
      common /advllu/ advllu1,advllu2
      common /advcri/ advcri1,advcri2
      common /advnie/ advnie1,advnie2
      common /advgra/ advgra1,advgra2
      common /advvap/ advvap1,advvap2
end module advecs