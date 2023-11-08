!     revisado 12/02/99
!     Esta subrutina calcula los terminos inomogeneos para las velocidades
subroutine inomo(i,j,k,dden0z)
   USE cant01
   USE constants
   USE perdim
   USE permic
   USE const
   USE estbas
   USE turbvar
   implicit none
   include 'fuvw.i'
   include 'inomo.i'

   integer, intent(in) :: i,j,k
   real, intent(in) :: dden0z
   
   dvelxx=(U2(i-2,j,k)-U2(i+2,j,k))+8.*(U2(i+1,j,k)-U2(i-1,j,k))
   dvelxy=(U2(i,j-2,k)-U2(i,j+2,k))+8.*(U2(i,j+1,k)-U2(i,j-1,k))
   dvelxz=(U2(i,j,k-2)-U2(i,j,k+2))+8.*(U2(i,j,k+1)-U2(i,j,k-1))
   dvelyx=(V2(i-2,j,k)-V2(i+2,j,k))+8.*(V2(i+1,j,k)-V2(i-1,j,k))
   dvelyy=(V2(i,j-2,k)-V2(i,j+2,k))+8.*(V2(i,j+1,k)-V2(i,j-1,k))
   dvelyz=(V2(i,j,k-2)-V2(i,j,k+2))+8.*(V2(i,j,k+1)-V2(i,j,k-1))
   dvelzx=(W2(i-2,j,k)-W2(i+2,j,k))+8.*(W2(i+1,j,k)-W2(i-1,j,k))
   dvelzy=(W2(i,j-2,k)-W2(i,j+2,k))+8.*(W2(i,j+1,k)-W2(i,j-1,k))
   dvelzz=(W2(i,j,k-2)-W2(i,j,k+2))+8.*(W2(i,j,k+1)-W2(i,j,k-1))

   diverx=(U2(i,j,k)*dvelxx+V2(i,j,k)*dvelxy)+W2(i,j,k)*dvelxz
   divery=(U2(i,j,k)*dvelyx+V2(i,j,k)*dvelyy)+W2(i,j,k)*dvelyz
   diverz=(U2(i,j,k)*dvelzx+V2(i,j,k)*dvelzy)+W2(i,j,k)*dvelzz

   a1=(KM1*DD(1,1)+KM2*DD(1,2))+KM3*DD(1,3)
   a2=KMM*((D1(1)+D2(1))+D3(1))
   a3=KMM*DD(1,3)*dden0z
   turbulx=cteturb*((a1+a2)+a3)


   a1=(KM1*DD(2,1)+KM2*DD(2,2))+KM3*DD(2,3)
   a2=KMM*((D1(2)+D2(2))+D3(2))
   a3=KMM*DD(2,3)*dden0z
   turbuly=cteturb*((a1+a2)+a3)

   a1=(KM1*DD(3,1)+KM2*DD(3,2))+KM3*DD(3,3)
   a2=KMM*((D1(3)+D2(3))+D3(3))
   a3=KMM*DD(3,3)*dden0z
   turbulz=cteturb*((a1+a2)+a3)

   !$$
   grave=G*(Titaa1(i,j,k)/Tita0(k)+(AA*Qvap1(i,j,k)-&
      Qgot1(i,j,k)-Qllu1(i,j,k)-Qcri1(i,j,k)-&
      Qnie1(i,j,k)-Qgra1(i,j,k))/Den0(k))
   !      grave=G*(Titaa1(i,j,k)/Tita0(k))

   fu(i,j,k)=turbulx/dx8-diverx/dx12
   fv(i,j,k)=turbuly/dx8-divery/dx12
   fw(i,j,k)=turbulz/dx8-diverz/dx12+grave

   !     agregado para la P (23/8/97)
   laplap=(Pres2(i+1,j,k)+Pres2(i-1,j,k))+(Pres2(i,j+1,k)+&
      Pres2(i,j-1,k))+Pres2(i,j,k+1)+Pres2(i,j,k-1)-&
      6.*Pres2(i,j,k)
   fp(i,j,k)=cteturb*KMM/dx2*laplap

   !##
   !      if (i.eq.2.and. (j.eq.5.or.j.eq.28) .and. k.eq.1) then
   !	write(*,*) 'in',i,j,k
   !        write(*,'(5g16.8)') cteturb,fu(i,j,k),fv(i,j,k),fw(i,j,k),
   !     &	fp(i,j,k)
   !        write(*,'(4g16.8)') fw(i,j,k),turbulz,diverz,grave
   !        write(*,'(7g16.8)') G,Titaa1(i,j,k),Tita0(k),AA,Qvap1(i,j,k),
   !     &         Qgot1(i,j,k),Den0(k)
   !      endif

   !##

   return
end


