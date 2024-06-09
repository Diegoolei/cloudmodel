!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!> TEMPE01
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


subroutine tempot(i,j,k,dden0z,Fcal)
   !> Fcalo es el calor liberado por cambio de fase, por unidad de masa de aire
   USE cant01
   USE dimen
   USE perdim
   USE const
   USE estbas
   USE turbvar
   USE tempe01
   implicit none

   real, intent(in) :: dden0z,Fcal
   integer, intent(in) :: i,j,k

   adv(1)=(U2(i,j,k)+UU(k))*(Titaa1(i+1,j,k)-Titaa1(i-1,j,k))
   adv(2)=(V2(i,j,k)+VV(k))*(Titaa1(i,j+1,k)-Titaa1(i,j-1,k))
   adv(3)=W2(i,j,k)*(Titaa1(i,j,k+1)-Titaa1(i,j,k-1))

   advec=-((adv(1)+adv(2))+adv(3))

   verti=-W2(i,j,k)*(Tita0(k+1)-Tita0(k-1))

   calor=Fcal*Tita0(k)/(Temp0(k)*Cp)

   dtita(1)=Titaa1(i+1,j,k)-Titaa1(i-1,j,k)
   dtita(2)=Titaa1(i,j+1,k)-Titaa1(i,j-1,k)
   dtita(3)=Titaa1(i,j,k+1)-Titaa1(i,j,k-1)

   escal=(dtita(1)*KM1+dtita(2)*KM2)+dtita(3)*KM3

   lapla=(Titaa1(i+1,j,k)+Titaa1(i-1,j,k))+(Titaa1(i,j+1,k)+&
      Titaa1(i,j-1,k))+Titaa1(i,j,k+1)+Titaa1(i,j,k-1)-&
      6*Titaa1(i,j,k)
   lapla=lapla+(Tita0(k+1)+Tita0(k-1)-2.*Tita0(k))

   turden=dden0z*(Tita0(k+1)-Tita0(k-1))

   turbul=3.*cteturb/dx8*(escal+KMM*(4.*lapla+turden))

   Titaa2(i,j,k)=dt1*((advec+verti)/dx2+turbul+calor)+&
      Titaa1(i,j,k)

   !     control de locura
   if (abs(Titaa2(i,j,k)).gt.30) then
      stop
   endif

   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> CORREC01
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de gotitas
 !     revisada 28/01/99
subroutine corgot
   USE dimen
   USE permic
   USE lmngot
   USE corgot_vars
   implicit none

   neg1=0.
   pos1=0.
   do n=ngot(1),ngot(2)
      do l=lgot(1),lgot(2)
         do m=mgot(1),mgot(2)
            if (Qgot2(l,m,n).lt.0.) then
               neg1=neg1+Qgot2(l,m,n)
               Qgot2(l,m,n)=0
            else
               pos1=pos1+Qgot2(l,m,n)
            endif
         end do
      end do
   end do

   if(pos1.le.-neg1) then
      do l=lgot(1),lgot(2)
         do m=mgot(1),mgot(2)
            do n=ngot(1),ngot(2)
               Qgot2(l,m,n)=0.
            end do
         end do
      end do

      if (-neg1.gt.1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do l=lgot(1),lgot(2)
         do m=mgot(1),mgot(2)
            do n=ngot(1),ngot(2)
               if(Qgot2(l,m,n).gt.0) then
                  Qgot2(l,m,n)=Qgot2(l,m,n)*(1.+aux1)
               endif
            end do
         end do
      end do

   endif

   return
end

 !********************************************************************

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de gotas
 !     revisada 8/01/99
subroutine corllu
   USE dimen
   USE permic
   USE lmnllu
   USE corgot_vars
   implicit none

   neg1=0.
   pos1=0.
   do n=nllu(1),nllu(2)
      do l=lllu(1),lllu(2)
         do m=mllu(1),mllu(2)
            if (Qllu2(l,m,n).lt.0.) then
               neg1=neg1+Qllu2(l,m,n)
               Qllu2(l,m,n)=0
            else
               pos1=pos1+Qllu2(l,m,n)
            endif
         end do
      end do
   end do
   if(pos1.le.-neg1) then
      do l=lllu(1),lllu(2)
         do m=mllu(1),mllu(2)
            do n=nllu(1),nllu(2)
               Qllu2(l,m,n)=0.
            end do
         end do
      end do

      if (-neg1.gt.1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do l=lllu(1),lllu(2)
         do m=mllu(1),mllu(2)
            do n=nllu(1),nllu(2)
               if(Qllu2(l,m,n).gt.0) then
                  Qllu2(l,m,n)=Qllu2(l,m,n)*(1.+aux1)
               endif
            end do
         end do
      end do

   endif

   return
end

 !********************************************************************
 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de cristales
 !     revisada 8/01/99
subroutine corcri
   USE dimen
   USE permic
   USE lmncri
   USE corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do n=ncri(1),ncri(2)
      do l=lcri(1),lcri(2)
         do m=mcri(1),mcri(2)
            if (Qcri2(l,m,n).lt.0.) then
               neg1=neg1+Qcri2(l,m,n)
               Qcri2(l,m,n)=0
            else
               pos1=pos1+Qcri2(l,m,n)
            endif
         end do
      end do
   end do

   if(pos1.le.-neg1) then
      do l=lcri(1),lcri(2)
         do m=mcri(1),mcri(2)
            do n=ncri(1),ncri(2)
               Qcri2(l,m,n)=0.
            end do
         end do
      end do

      if (-neg1.gt.1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do l=lcri(1),lcri(2)
         do m=mcri(1),mcri(2)
            do n=ncri(1),ncri(2)
               if(Qcri2(l,m,n).gt.0) then
                  Qcri2(l,m,n)=Qcri2(l,m,n)*(1.+aux1)
               endif
            end do
         end do
      end do

   endif

   return
end

 !********************************************************************
 !$$
 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de nieve
 !     revisada 11/01/99
subroutine cornie
   USE dimen
   USE permic
   USE lmnnie
   USE corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do n=nnie(1),nnie(2)
      do l=lnie(1),lnie(2)
         do m=mnie(1),mnie(2)
            if (Qnie2(l,m,n).lt.0.) then
               neg1=neg1+Qnie2(l,m,n)
               Qnie2(l,m,n)=0
            else
               pos1=pos1+Qnie2(l,m,n)
            endif
         end do
      end do
   end do

   if(pos1.le.-neg1) then
      do l=lnie(1),lnie(2)
         do m=mnie(1),mnie(2)
            do n=nnie(1),nnie(2)
               Qnie2(l,m,n)=0.
            end do
         end do
      end do

      if (-neg1.gt.1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do l=lnie(1),lnie(2)
         do m=mnie(1),mnie(2)
            do n=nnie(1),nnie(2)
               if(Qnie2(l,m,n).gt.0) then
                  Qnie2(l,m,n)=Qnie2(l,m,n)*(1.+aux1)
               endif
            end do
         end do
      end do

   endif

   return
end

 !********************************************************************
 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de granizos
 !     revisada 2/02/99
subroutine corgra
   USE dimen
   USE permic
   USE lmngra
   USE corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do n=ngra(1),ngra(2)
      do l=lgra(1),lgra(2)
         do m=mgra(1),mgra(2)
            if (Qgra2(l,m,n).lt.0.) then
               neg1=neg1+Qgra2(l,m,n)
               Qgra2(l,m,n)=0
            else
               pos1=pos1+Qgra2(l,m,n)
            endif
         end do
      end do
   end do

   if(pos1.le.-neg1) then
      do l=lgra(1),lgra(2)
         do m=mgra(1),mgra(2)
            do n=ngra(1),ngra(2)
               Qgra2(l,m,n)=0.
            end do
         end do
      end do

      if (-neg1.gt.1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do l=lgra(1),lgra(2)
         do m=mgra(1),mgra(2)
            do n=ngra(1),ngra(2)
               if(Qgra2(l,m,n).gt.0) then
                  Qgra2(l,m,n)=Qgra2(l,m,n)*(1.+aux1)
               endif
            end do
         end do
      end do

   endif

   return
end

 !********************************************************************

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de vapor
 !     revisada 8/01/99
subroutine corvap(Qvapneg)
   USE dimen
   USE permic
   USE estbas
   USE corvap_vars
   implicit none

   real(8), intent(in) :: Qvapneg

   do k=1,nz1
      dq=Qvapneg*Qvaprel(k)/nx1**2.
      do i=1,nx1
         do j=1,nx1
            Qvap2(i,j,k)=Qvap2(i,j,k)+dq
            if (Qvap2(i,j,k)+Qvap0(k).lt.0) Qvap2(i,j,k)=-Qvap0(k)
         end do
      end do
   end do

   return
end

 !********************************************************************

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de aerosoles
 !     revisada 14/09/99
subroutine coraer(aerneg)
   USE dimen
   USE permic
   USE estbas
   USE coraer_vars
   implicit none

   real(8), intent(in) :: aerneg

   do k=1,nz1
      dq=aerneg*aerrel(k)/nx1**2.
      do i=1,nx1
         do j=1,nx1
            aer2(i,j,k)=aer2(i,j,k)+dq
            if (aer2(i,j,k)+aer0(k).lt.0) aer2(i,j,k)=-aer0(k)
         end do
      end do
   end do

   return
end


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !     Revision 28/04/98
subroutine daeros(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE daeros_vars
   implicit none

   integer, intent(in) :: l,m,n

   daer(1)=aer1(l+1,m,n)-aer1(l-1,m,n)
   daer(2)=aer1(l,m+1,n)-aer1(l,m-1,n)
   daer(3)=aer1(l,m,n+1)-aer1(l,m,n-1)


   adv(1)=(((U2(l+1,m,n)+U2(l,m,n))*(aer1(l+1,m,n)+aer1(l,m,n)))-&
      ((U2(l-1,m,n)+U2(l,m,n))*(aer1(l-1,m,n)+aer1(l,m,n))))/4.
   adv(1)=adv(1)+daer(1)/2.*UU(n)

   adv(2)=(((V2(l,m+1,n)+V2(l,m,n))*(aer1(l,m+1,n)+aer1(l,m,n)))&
      -((V2(l,m-1,n)+V2(l,m,n))*(aer1(l,m-1,n)+aer1(l,m,n))))/4.
   adv(2)=adv(2)+daer(2)/2.*VV(n)

   advaer2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (aer1(l,m,n)+aer1(l,m,n+1))/4.

   adv(3)=advaer2(l,m)-advaer1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))

   verti=-((W2(l,m,n+1)+W2(l,m,n))*(aer0(n+1)+aer0(n))-&
      (W2(l,m,n-1)+W2(l,m,n))*(aer0(n-1)+aer0(n)))/4.

   !## agregado
   aux=-((U2(l+1,m,n)-U2(l-1,m,n))+(V2(l,m+1,n)-V2(l,m-1,n)))*&
      aer0(n)/2.

   escal=daer(1)*KM1+daer(2)*KM2+daer(3)*KM3

   lapla=aer1(l+1,m,n)+aer1(l,m+1,n)+aer1(l,m,n+1)+&
      aer1(l-1,m,n)+aer1(l,m-1,n)+aer1(l,m,n-1)-&
      6.*aer1(l,m,n)

   lapla=((aer1(l+1,m,n)+aer1(l-1,m,n))+(aer1(l,m+1,n)+&
      aer1(l,m-1,n)))+aer1(l,m,n-1)+aer1(l,m,n+1)-&
      6.*aer1(l,m,n)


   lapla=lapla+(aer0(n+1)+aer0(n-1)-2.*aer0(n))

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   aer2(l,m,n)=dt1*((advec+verti+aux)/dx1+turbul)+aer1(l,m,n)
   return
end

 !********************************************************************
 !     revisada 28/04/98
subroutine dgotit(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dgotit_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqgot(1)=Qgot1(l+1,m,n)-Qgot1(l-1,m,n)
   dqgot(2)=Qgot1(l,m+1,n)-Qgot1(l,m-1,n)
   dqgot(3)=Qgot1(l,m,n+1)-Qgot1(l,m,n-1)

   adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qgot1(l+1,m,n)+Qgot1(l,m,n))&
      -(U2(l-1,m,n)+U2(l,m,n))*(Qgot1(l-1,m,n)+Qgot1(l,m,n)))/4.
   adv(1)=adv(1)+dqgot(1)/2.*UU(n)

   adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qgot1(l,m+1,n)+Qgot1(l,m,n))&
      -(V2(l,m-1,n)+V2(l,m,n))*(Qgot1(l,m-1,n)+Qgot1(l,m,n)))/4.
   adv(2)=adv(2)+dqgot(2)/2.*VV(n)

   advgot2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qgot1(l,m,n)+Qgot1(l,m,n+1))/4.
   if ((advgot2(l,m)-advgot1(l,m))*dt1/dx1.gt.Qgot1(l,m,n) .and.&
      W2(l,m,n).gt.0) then
      advgot2(l,m)=advgot1(l,m)+Qgot1(l,m,n)*dx1/dt1
   endif
   adv(3)=advgot2(l,m)-advgot1(l,m)


   advec=-(adv(1)+adv(2)+adv(3))/dx1

   escal=dqgot(1)*KM1+dqgot(2)*KM2+dqgot(3)*KM3


   lapla=Qgot1(l+1,m,n)+Qgot1(l,m+1,n)+Qgot1(l,m,n+1)+&
      Qgot1(l-1,m,n)+Qgot1(l,m-1,n)+Qgot1(l,m,n-1)-&
      6.*Qgot1(l,m,n)

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   Qgot2(l,m,n)=dt1*(advec+turbul)+Qgot1(l,m,n)

   return
end
 !********************************************************************
 !     Revision 28/04/98
subroutine dvapor(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dvapor_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqvap(1)=Qvap1(l+1,m,n)-Qvap1(l-1,m,n)
   dqvap(2)=Qvap1(l,m+1,n)-Qvap1(l,m-1,n)
   dqvap(3)=Qvap1(l,m,n+1)-Qvap1(l,m,n-1)

   adv(1)=(((U2(l+1,m,n)+U2(l,m,n))*&
      (Qvap1(l+1,m,n)+Qvap1(l,m,n)))-&
      ((U2(l-1,m,n)+U2(l,m,n))*&
      (Qvap1(l-1,m,n)+Qvap1(l,m,n))))/4.
   adv(1)=adv(1)+dqvap(1)/2.*UU(n)

   adv(2)=(((V2(l,m+1,n)+V2(l,m,n))*&
      (Qvap1(l,m+1,n)+Qvap1(l,m,n)))-&
      ((V2(l,m-1,n)+V2(l,m,n))*&
      (Qvap1(l,m-1,n)+Qvap1(l,m,n))))/4.
   adv(2)=adv(2)+dqvap(2)/2.*VV(n)


   advvap2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qvap1(l,m,n)+Qvap1(l,m,n+1))/4.

   adv(3)=advvap2(l,m)-advvap1(l,m)

   advec=-((adv(1)+adv(2))+adv(3))

   verti=-((W2(l,m,n+1)+W2(l,m,n))*(Qvap0(n+1)+Qvap0(n))-&
      (W2(l,m,n-1)+W2(l,m,n))*(Qvap0(n-1)+Qvap0(n)))/4.

   !## agregado
   aux=-(U2(l+1,m,n)-U2(l-1,m,n)+V2(l,m+1,n)-V2(l,m-1,n))*&
      Qvap0(n)/2.

   escal=dqvap(1)*KM1+dqvap(2)*KM2+dqvap(3)*KM3

   lapla=((Qvap1(l+1,m,n)+Qvap1(l-1,m,n))+(Qvap1(l,m+1,n)+&
      Qvap1(l,m-1,n)))+Qvap1(l,m,n+1)+Qvap1(l,m,n-1)-&
      6.*Qvap1(l,m,n)

   lapla=lapla+(Qvap0(n+1)+Qvap0(n-1)-2.*Qvap0(n))

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   !      Qvap2(l,m,n)=dt1*((advec+verti)/dx1+turbul)+Qvap1(l,m,n)
   Qvap2(l,m,n)=dt1*((advec+verti+aux)/dx1+turbul)+Qvap1(l,m,n)

   aux=dt1/dx1
   aux=aux*(verti+advec)+dt1*turbul+Qvap1(l,m,n)

   return
end
 !**********************************************************
 !     Revision 28/04/98
subroutine dlluvi(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dlluvi_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqllu(1)=Qllu1(l+1,m,n)-Qllu1(l-1,m,n)
   dqllu(2)=Qllu1(l,m+1,n)-Qllu1(l,m-1,n)
   dqllu(3)=Qllu1(l,m,n+1)-Qllu1(l,m,n-1)

   adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qllu1(l+1,m,n)+Qllu1(l,m,n))&
      -(U2(l-1,m,n)+U2(l,m,n))*(Qllu1(l-1,m,n)+Qllu1(l,m,n)))/4.
   adv(1)=adv(1)+dqllu(1)/2.*UU(n)

   adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qllu1(l,m+1,n)+Qllu1(l,m,n))-&
      (V2(l,m-1,n)+V2(l,m,n))*(Qllu1(l,m-1,n)+Qllu1(l,m,n)))/4.
   adv(2)=adv(2)+dqllu(2)/2.*VV(n)

   advllu2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qllu1(l,m,n)+Qllu1(l,m,n+1))/4.

   adv(3)=advllu2(l,m)-advllu1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))

   escal=dqllu(1)*KM1+dqllu(2)*KM2+dqllu(3)*KM3

   lapla=Qllu1(l+1,m,n)+Qllu1(l,m+1,n)+Qllu1(l,m,n+1)+&
      Qllu1(l-1,m,n)+Qllu1(l,m-1,n)+Qllu1(l,m,n-1)-&
      6.*Qllu1(l,m,n)

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   !***  termino de sedimentacion

   Qllus=(Qllu1(l,m,n+1)+Qllu1(l,m,n))/2.
   Qllui=(Qllu1(l,m,n-1)+Qllu1(l,m,n))/2.
   Rms=(Qllus/cteqllu)**.25
   Rmm=(Qllu1(l,m,n)/cteqllu)**.25
   Rmi=(Qllui/cteqllu)**.25
   Vtllus=(Av(2*n+1)*Rms**.8+Av(2*n)*Rmm**.8)/2.
   Vtllui=(Av(2*n-1)*Rmi**.8+Av(2*n)*Rmm**.8)/2.
   if (n.eq.1) then
      Vtllui=Av(2*n)*Rmm**.8
      Qllui=Qllu1(l,m,n)
   endif

   sedim=gam4p8/6.*(Qllus*Vtllus-Qllui*Vtllui)
   !***

   Qllu2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qllu1(l,m,n)
   return
end

 !**********************************************************
 !     Revision 29/12/98
subroutine dcrist(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dcrist_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqcri(1)=Qcri1(l+1,m,n)-Qcri1(l-1,m,n)
   dqcri(2)=Qcri1(l,m+1,n)-Qcri1(l,m-1,n)
   dqcri(3)=Qcri1(l,m,n+1)-Qcri1(l,m,n-1)

   adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qcri1(l+1,m,n)+Qcri1(l,m,n))&
      -(U2(l-1,m,n)+U2(l,m,n))*(Qcri1(l-1,m,n)+Qcri1(l,m,n)))/4.
   adv(1)=adv(1)+dqcri(1)/2.*UU(n)

   adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qcri1(l,m+1,n)+Qcri1(l,m,n))-&
      (V2(l,m-1,n)+V2(l,m,n))*(Qcri1(l,m-1,n)+Qcri1(l,m,n)))/4.
   adv(2)=adv(2)+dqcri(2)/2.*VV(n)

   advcri2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qcri1(l,m,n)+Qcri1(l,m,n+1))/4.

   adv(3)=advcri2(l,m)-advcri1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))/dx1

   escal=dqcri(1)*KM1+dqcri(2)*KM2+dqcri(3)*KM3

   lapla=Qcri1(l+1,m,n)+Qcri1(l,m+1,n)+Qcri1(l,m,n+1)+&
      Qcri1(l-1,m,n)+Qcri1(l,m-1,n)+Qcri1(l,m,n-1)-&
      6.*Qcri1(l,m,n)

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   Qcri2(l,m,n)=dt1*(advec+turbul)+Qcri1(l,m,n)

   return
end
 !**********************************************************
 !$$
 !     Revision 7/06/99
subroutine dnieve(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dnieve_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqnie(1)=Qnie1(l+1,m,n)-Qnie1(l-1,m,n)
   dqnie(2)=Qnie1(l,m+1,n)-Qnie1(l,m-1,n)
   dqnie(3)=Qnie1(l,m,n+1)-Qnie1(l,m,n-1)

   adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qnie1(l+1,m,n)+Qnie1(l,m,n))&
      -(U2(l-1,m,n)+U2(l,m,n))*(Qnie1(l-1,m,n)+Qnie1(l,m,n)))/4.
   adv(1)=adv(1)+dqnie(1)/2.*UU(n)

   adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qnie1(l,m+1,n)+Qnie1(l,m,n))-&
      (V2(l,m-1,n)+V2(l,m,n))*(Qnie1(l,m-1,n)+Qnie1(l,m,n)))/4.
   adv(2)=adv(2)+dqnie(2)/2.*VV(n)

   advnie2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qnie1(l,m,n)+Qnie1(l,m,n+1))/4.

   adv(3)=advnie2(l,m)-advnie1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))

   escal=dqnie(1)*KM1+dqnie(2)*KM2+dqnie(3)*KM3

   lapla=Qnie1(l+1,m,n)+Qnie1(l,m+1,n)+Qnie1(l,m,n+1)+&
      Qnie1(l-1,m,n)+Qnie1(l,m-1,n)+Qnie1(l,m,n-1)-&
      6.*Qnie1(l,m,n)

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   !***  termino de sedimentacion

   Qnies=(Qnie1(l,m,n+1)+Qnie1(l,m,n))/2.
   Qniei=(Qnie1(l,m,n-1)+Qnie1(l,m,n))/2.

   sedim=Vtnie(2*n+1)*Qnies-Vtnie(2*n-1)*Qniei
   !***

   Qnie2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qnie1(l,m,n)

   return
end
 !**********************************************************
 !     Revision 1/02/99
subroutine dgrani(l,m,n)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE advecs
   USE turbvar
   USE dgrani_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqgra(1)=Qgra1(l+1,m,n)-Qgra1(l-1,m,n)
   dqgra(2)=Qgra1(l,m+1,n)-Qgra1(l,m-1,n)
   dqgra(3)=Qgra1(l,m,n+1)-Qgra1(l,m,n-1)

   adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qgra1(l+1,m,n)+Qgra1(l,m,n))&
      -(U2(l-1,m,n)+U2(l,m,n))*(Qgra1(l-1,m,n)+Qgra1(l,m,n)))/4.
   adv(1)=adv(1)+dqgra(1)/2.*UU(n)

   adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qgra1(l,m+1,n)+Qgra1(l,m,n))-&
      (V2(l,m-1,n)+V2(l,m,n))*(Qgra1(l,m-1,n)+Qgra1(l,m,n)))/4.
   adv(2)=adv(2)+dqgra(2)/2.*VV(n)

   advgra2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*&
      (Qgra1(l,m,n)+Qgra1(l,m,n+1))/4.

   adv(3)=advgra2(l,m)-advgra1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))

   escal=dqgra(1)*KM1+dqgra(2)*KM2+dqgra(3)*KM3

   lapla=Qgra1(l+1,m,n)+Qgra1(l,m+1,n)+Qgra1(l,m,n+1)+&
      Qgra1(l-1,m,n)+Qgra1(l,m-1,n)+Qgra1(l,m,n-1)-&
      6.*Qgra1(l,m,n)

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   !***  termino de sedimentacion
   Qgras=(Qgra1(l,m,n+1)+Qgra1(l,m,n))/2.
   Qgrai=(Qgra1(l,m,n-1)+Qgra1(l,m,n))/2.
   Rms=(Qgras/cteqgra)**.25
   Rmm=(Qgra1(l,m,n)/cteqgra)**.25
   Rmi=(Qgrai/cteqgra)**.25
   Vtgras=(Vtgra0(2*n+1)*Rms**.8+Vtgra0(2*n)*Rmm**.8)/2.
   Vtgrai=(Vtgra0(2*n-1)*Rmi**.8+Vtgra0(2*n)*Rmm**.8)/2.
   if (n.eq.1) then
      Vtgrai=Vtgra0(2*n)*Rmm**.8
      Qgrai=Qgra1(l,m,n)
   endif

   sedim=gam4p8/6.*(Qgras*Vtgras-Qgrai*Vtgrai)
   !***


   Qgra2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qgra1(l,m,n)


   return
end


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> VELPRE01
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


 !> Calcula la evolucion del la presion y las velocidades con un paso de tiempo menor lt3
 !> Las cantidades 1 son las presentes en el paso grande y las 2 son las del paso futuro, las 3 son auxiliares
 !> Le resta la perturbacion promedio

subroutine velpre
   USE cant01
   USE dimen
   USE perdim
   USE const
   USE estbas
   USE velpre01
   USE p3v3
   USE fuvw
   implicit none

   call velpre01_init()

   do i=0,nx1+1
      do j=0,nx1+1
         do k=0,nz1
            U2(i,j,k)=U1(i,j,k)
            V2(i,j,k)=V1(i,j,k)
            W2(i,j,k)=W1(i,j,k)
            Pres2(i,j,k)=Pres1(i,j,k)
         end do
      end do
   end do

   do t=1,lt3
      presprom=0.
      do k=1,nz1-1
         presi=-Cp*Tita0(k)*(1.+.61*Qvap0(k)/Den0(k))
         vel0=Tita0(k)*(Den0(k)+.61*Qvap0(k))
         vel1=Tita0(k-1)*(Den0(k-1)+.61*Qvap0(k-1))
         vel2=Tita0(k+1)*(Den0(k+1)+.61*Qvap0(k+1))
         vel3=cc2(k)/presi/vel0
         do i=1,nx1
            do j=1,nx1

               dprex=Pres2(i+1,j,k)-Pres2(i-1,j,k)
               dprey=Pres2(i,j+1,k)-Pres2(i,j-1,k)
               dprez=Pres2(i,j,k+1)-Pres2(i,j,k-1)

               presix=presi*dprex/dx2
               presiy=presi*dprey/dx2
               presiz=presi*dprez/dx2

               U3(i,j,k)=dt3*(presix+fu(i,j,k))+U2(i,j,k)
               V3(i,j,k)=dt3*(presiy+fv(i,j,k))+V2(i,j,k)
               W3(i,j,k)=dt3*(presiz+fw(i,j,k))+W2(i,j,k)


               dvx=vel0*(U2(i+1,j,k)-U2(i-1,j,k))
               dvy=vel0*(V2(i,j+1,k)-V2(i,j-1,k))
               if (k.eq.1) then
                  !      dvz=tiene 80% de (W2(2)-W2(1) y 20% de (W2(1)-W2(0)
                  dvz=(.8*vel2*W2(i,j,k+1)-.8*vel1*W2(i,j,k))*2.
               else
                  dvz=vel2*W2(i,j,k+1)-vel1*W2(i,j,k-1)
               endif

               diver=vel3*((dvx+dvy)+dvz)/dx2

               !      modificado para agrega turbulencia en la P 23/8/97
               Pres3(i,j,k)=dt3*(diver+fp(i,j,k))+Pres2(i,j,k)
            end do
         end do
      end do

      !*      redefiniciones y contornos
      do i=1,nx1
         do j=1,nx1
            Pres3(i,j,0)=Pres3(i,j,1)
            Pres3(i,j,nz1)=Pres3(i,j,nz1-1)
         end do
      end do
      do i=1,nx1
         do k=0,nz1
            Pres3(i,0,k)=Pres3(i,1,k)
            Pres3(i,nx1+1,k)=Pres3(i,nx1,k)
            Pres3(0,i,k)=Pres3(1,i,k)
            Pres3(nx1+1,i,k)=Pres3(nx1,i,k)
         end do
      end do

      do i=1,nx1
         do j=1,nx1
            do k=1,nz1-1
               if (k.eq.1) then
                  U2(i,j,k)=U3(i,j,k)-kkk*&
                     (2.*U3(i,j,k)-U3(i,j,k+1))
                  V2(i,j,k)=V3(i,j,k)-kkk*&
                     (2.*V3(i,j,k)-V3(i,j,k+1))
                  W2(i,j,k)=W3(i,j,k)-kkk*&
                     (2.*W3(i,j,k)-W3(i,j,k+1))
               else
                  U2(i,j,k)=U3(i,j,k)
                  V2(i,j,k)=V3(i,j,k)
                  W2(i,j,k)=W3(i,j,k)
               endif
               Pres2(i,j,k)=prom1*Pres3(i,j,k)+prom*(&
                  ((Pres3(i+1,j,k)+ Pres3(i-1,j,k))+&
                  (Pres3(i,j+1,k)+Pres3(i,j-1,k)))+&
                  Pres3(i,j,k+1)+Pres3(i,j,k-1))
               presprom=Pres2(i,j,k)+presprom
            end do

            U2(i,j,0)=0
            V2(i,j,0)=0
            W2(i,j,0)=0
            Pres2(i,j,0)=Pres2(i,j,1)
            U2(i,j,nz1)=U2(i,j,nz1-1)
            V2(i,j,nz1)=V2(i,j,nz1-1)
            W2(i,j,nz1)=W2(i,j,nz1-1)
            Pres2(i,j,nz1)=Pres2(i,j,nz1-1)
         end do
      end do
      do i=1,nx1
         do k=0,nz1
            U2(0,i,k)=U2(1,i,k)
            V2(0,i,k)=V2(1,i,k)
            W2(0,i,k)=W2(1,i,k)
            Pres2(0,i,k)=Pres2(1,i,k)
            U2(nx1+1,i,k)=U2(nx1,i,k)
            V2(nx1+1,i,k)=V2(nx1,i,k)
            W2(nx1+1,i,k)=W2(nx1,i,k)
            Pres2(nx1+1,i,k)=Pres2(nx1,i,k)
            U2(i,0,k)=U2(i,1,k)
            V2(i,0,k)=V2(i,1,k)
            W2(i,0,k)=W2(i,1,k)
            Pres2(i,0,k)=Pres2(i,1,k)
            U2(i,nx1+1,k)=U2(i,nx1,k)
            V2(i,nx1+1,k)=V2(i,nx1,k)
            W2(i,nx1+1,k)=W2(i,nx1,k)
            Pres2(i,nx1+1,k)=Pres2(i,nx1,k)
         end do
      end do

      presprom=presprom/nnn
      do i=0,nx1+1
         do j=0,nx1+1
            do k=0,nz1
               Pres2(i,j,k)=Pres2(i,j,k)-presprom
            end do
         end do
      end do

      if (t.eq.lt3/2) then
         do i=0,nx1+1
            do j=0,nx1+1
               do k=0,nz1
                  U1(i,j,k)=U2(i,j,k)
                  V1(i,j,k)=V2(i,j,k)
                  W1(i,j,k)=W2(i,j,k)
                  Pres1(i,j,k)=Pres2(i,j,k)
               end do
            end do
         end do
      endif

   end do

   !**********************************************************
   !*    suavizado

   call filtro(Pres1,.15,.15,.1)

   call filtro(Pres2,.15,.15,.1)

   call filtro(U1,facx,facy,facz)
   call filtro(U2,facx,facy,facz)
   call filtro(V1,facx,facy,facz)
   call filtro(V2,facx,facy,facz)
   call filtro(W1,facx,facy,facz)
   call filtro(W2,facx,facy,facz)

   do i=1,nx1
      do j=1,nx1
         Pres1(i,j,0)=Pres1(i,j,1)
         Pres2(i,j,0)=Pres2(i,j,1)
      end do
   end do
   !**********************************************************

   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> TURBUA01
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine turbu1(kk)
   !> Esta subrutina calcula los Dnm para cada plano Z
   USE cant01
   USE dimen
   USE perdim
   USE const
   USE turbvar1
   USE turbu1_vars
   implicit none

   integer, intent(in) :: kk

   k=kk+1
   do i=0,nx1+1
      do j=0,nx1+1

         if (kk.eq.1) then
            do n=1,2
               do m=1,n
                  D(n,m,i,j,1)=0.
               end do
            end do
            D(3,1,i,j,1)=U2(i,j,1)
            D(3,2,i,j,1)=V2(i,j,1)
            D(3,3,i,j,1)=W2(i,j,1)*2./3.
            do n=1,3
               do m=1,n
                  D(m,n,i,j,1)=D(n,m,i,j,1)
               end do
            end do
            do lx=-1,1
               do ly=-1,1
                  do lz=-1,1
                     ldis=abs(lx)+abs(ly)+abs(lz)
                     if (ldis.le.1) then
                        vel(1,lx,ly,lz)=U2(lx+i,ly+j,lz+1)
                        vel(2,lx,ly,lz)=V2(lx+i,ly+j,lz+1)
                        vel(3,lx,ly,lz)=W2(lx+i,ly+j,lz+1)
                     endif
                  end do
               end do
            end do
            !     calculo de Dij
            do n=1,3
               dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
               dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
               dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
            end do
            do n=1,3
               do m=1,n
                  D(n,m,i,j,2)=(dv(n,m)+dv(m,n))
                  D(m,n,i,j,2)=D(n,m,i,j,2)
                  if (n.eq.m) D(n,n,i,j,2)=2./3.*D(n,n,i,j,2)
               end do
            end do
         else
            do n=1,3
               do m=1,3
                  do lz=1,2
                     D(n,m,i,j,lz)=D(n,m,i,j,lz+1)
                  end do
               end do
            end do
         endif
         !*********************************************************

         !     Lectura de las velocidades necesarias
         do lx=-1,1
            do ly=-1,1
               do lz=-1,1
                  ldis=abs(lx)+abs(ly)+abs(lz)
                  if (ldis.le.1) then
                     vel(1,lx,ly,lz)=U2(lx+i,ly+j,lz+k)
                     vel(2,lx,ly,lz)=V2(lx+i,ly+j,lz+k)
                     vel(3,lx,ly,lz)=W2(lx+i,ly+j,lz+k)
                  endif
               end do
            end do
         end do
         !     calculo de Dij
         do n=1,3
            dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
            dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
            dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
         end do
         do n=1,3
            do m=1,n
               D(n,m,i,j,3)=(dv(n,m)+dv(m,n))
               D(m,n,i,j,3)=D(n,m,i,j,3)
               if (n.eq.m) D(n,n,i,j,3)=2./3.*D(n,n,i,j,3)
            end do
         end do
      end do
   end do

   return
end


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> TURBUB70
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !      Esta subrutina calcula las cantidades referida a los terminos de
 !      turbulencia: K, DK, DDij
 !      Dij viene de turbu1
 !     En realidad para tener los valores de las 4 cantidades falta
 !     multiplicarlas por:
 !                      KMM : cteturb*dx1/2
 !                      DK  : cteturb/4
 !                      Dij : 1/dx2 (simetrico)
 !                      DDij: 1/dx2**2
 !     revisado : 28/04/97

subroutine turbu2(i,j)
   USE dimen
   USE turbvar
   USE turbvar1
   USE turbu2_vars
   implicit none
   integer, intent(in) :: i,j
   real aux

   !     calculo de KM
   do lx=-1,1
      do ly=-1,1
         do lz=-1,1
            ldis=abs(lx)+abs(ly)+abs(lz)
            if (ldis.le.2) then
               call suma(aux,&
                  D(1,1,i+lx,j+ly,2+lz)**2.,&
                  D(2,2,i+lx,j+ly,2+lz)**2.,&
                  D(3,3,i+lx,j+ly,2+lz)**2.)
               sum=aux/2.
               call suma(aux,&
                  D(1,2,i+lx,j+ly,2+lz)**2.,&
                  D(1,3,i+lx,j+ly,2+lz)**2.,&
                  D(2,3,i+lx,j+ly,2+lz)**2.)
               sum=sum+aux

               KM(lx,ly,lz)=sum**.5
            endif
         end do
      end do
   end do



   !     calculo de las derivadas (sin la distancia abajo)
   KM1=KM(1,0,0)-KM(-1,0,0)
   KM2=KM(0,1,0)-KM(0,-1,0)
   KM3=KM(0,0,1)-KM(0,0,-1)
   do n=1,3
      D1(n)=D(n,1,i+1,j,2)-D(n,1,i-1,j,2)
      D2(n)=D(n,2,i,j+1,2)-D(n,2,i,j-1,2)
      D3(n)=D(n,3,i,j,3)-D(n,3,i,j,1)
   end do
   KMM=KM(0,0,0)
   do n=1,3
      do m=1,3
         DD(n,m)=D(n,m,i,j,2)
      end do
   end do

   return
end

 !*********************************************************************
subroutine suma(sum,a1,a2,a3)
   implicit none
   real a1,a2,a3,sum,aux
   integer j

   do j=1,2
      if (a1.gt.a2) then
         aux=a1
         a1=a2
         a2=aux
      endif
      if (a2.gt.a3) then
         aux=a2
         a2=a3
         a3=aux
      endif
   end do

   sum=(a1+a2)+a3
   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> NUCLEA91
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine nuclea(Qvap,Qliq,Naer,TT,rhoa,e1,esl,ess,rl,rs,Lvl,Lvs,l,m,n,Naux,auxl,auxs)
   USE cant01
   USE dimen
   USE const
   USE nuclea61
   implicit none

   real, intent(in) ::  Naer, rhoa, rs, Lvl, Lvs
   real, intent(inout) :: Qliq, Qvap, TT, e1, esl, ess, rl
   real, intent(inout) :: Naux, auxl, auxs   !     Numero de aesosoles
   integer, intent(in) :: l, m, n

   B=Lvl/Rv
   auxl=0.
   auxs=0.
   Naux=0.

   Rcri=5e-5
   if(TT.lt.T0) Rcri=Rcri-4e-5*(T0-TT)/40.
   if(Rcri.lt.1e-5) Rcri=1e-5

   !     nucleacion sobre cristales
   Tc=T0-TT
   if (Tc .gt. 0 .and. rs .gt. 0 .and. Naer.gt.0) then

      mcri=pi*Rcri**3./10.*rhocri

      Naux=Acri*exp(Bcri*Tc)
      if (Naux .gt. .9) Naux=.9
      Naux=Naux*Naer*1e6   ! en m3
      auxs=Naux*mcri

      if (auxs .gt. (Qvap-ess/Rv/TT*.95)) then
         auxs=(Qvap-ess/Rv/TT)*.95
         Naux=auxs/mcri     ! en m3
      endif
      Qvap=Qvap-auxs
      TT1=TT
      TT2=TT+(auxs*Lvs)/(Cp*rhoa)
      TT=TT2
      e1=Qvap*Rv*TT
      esl=esl*exp(B*(TT2-TT1)/TT2/TT1)
      ess=ess*exp(B*(TT2-TT1)/TT2/TT1)
      Naux=-Naux/1e6        ! en cm3

   endif

   !     nucleacion sobre gotitas
   if (e1.gt.esl) then
      s=0
      hhh=0
      xxx=0
      TT1=TT
      Ti=TT
      ei=e1
      esli=esl
      auxl=0.
      caux=B*esli/Ti
      !TODO: Check this loop
10    continue
      F0=Lvl/Rv*(esl/TT1-ei/Ti)+Cp*rhoa*(TT1-Ti)
      F0p=Cp*rhoa+B/TT1**2.*caux
      TT2=TT1-F0/F0p
      auxl=(ei/Ti-esl/TT2)/Rv
      Qliq1=Qliq+auxl
      e1=esl

      if (Qliq1.lt.0) then
         Qliq1=0.
         auxl=-Qliq
         TT2=esl/(e1/TT1-auxl*Rv)
         e1=(Qvap-auxl)*Rv*TT2
         hhh=1
         if (s.eq.1) then
            stop
         endif
      endif

      s=1
      esl=esl*exp(B*(TT2-TT1)/TT2/TT1)
      TT1=TT2
      rl=abs(e1-esl)/esl
      xxx=xxx+1
      if (rl.gt.1e-3 .and. hhh.eq.0) goto 10

      Qvap=Qvap-auxl
      Qliq=Qliq+auxl

      !*     control de mocos
      if (auxl.lt.0 .or. auxs.lt.0) then
         stop
      endif

      !*
      !      variacion en los aerosoles
      !      considerando que las nuevas gotitas tienen un radio Rgotmin

      Naux=Naux-auxl/(4./3.*pi*rhow*Rgotmin**3.)/1e6

      TT=TT2

   endif

   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> INOMO60
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !     revisado 12/02/99
 !     Esta subrutina calcula los terminos inomogeneos para las velocidades
subroutine inomo(i,j,k,dden0z)
   USE cant01
   USE dimen
   USE perdim
   USE permic
   USE const
   USE estbas
   USE turbvar
   USE fuvw
   USE inomo_var
   implicit none

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

   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> FILTRO01
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !     Esta subrutina filtra componentes de alta frecuencia espacial.
 !     El valor de la variable del punto j se filtra con los valores
 !     extrapolados linalmente de los puntos j-3 y j-1 y similares,
 !     pasando un polinomio de grado 4.

subroutine filtro(varia1,facx,facy,facz)
   USE dimen
   USE filtro01
   implicit none
   character*50 text
   REAL, DIMENSION(-3:NX1+3,-3:NX1+3,-2:NZ1+2), intent(inout) :: varia1
   real, intent(in) :: facx,facy,facz
   fact=1.-(facx+facy+facz)

   if (fact.lt.0.25) then
      stop
   endif

   !**********************************************************
   !     Redefiniciones y contornos

   do i=0,nx1+1
      do j=0,nx1+1
         do k=0,nz1
            varia2(i,j,k)=varia1(i,j,k)
         end do
      end do
   end do

   do k=0,nz1
      do i=0,nx1
         varia2(i,-1,k)=varia2(i,1,k)
         varia2(i,-2,k)=varia2(i,1,k)
         varia2(i,nx1+2,k)=varia2(i,nx1,k)
         varia2(i,nx1+3,k)=varia2(i,nx1,k)
         varia2(-1,i,k)=varia2(1,i,k)
         varia2(-2,i,k)=varia2(1,i,k)
         varia2(nx1+2,i,k)=varia2(nx1,i,k)
         varia2(nx1+3,i,k)=varia2(nx1,i,k)
      end do
   end do

   do i=1,nx1
      do j=1,nx1
         varia2(i,j,-1)=varia2(i,j,0)
         varia2(i,j,-2)=varia2(i,j,0)
         varia2(i,j,nz1+1)=varia2(i,j,nz1)
         varia2(i,j,nz1+2)=varia2(i,j,nz1)
      end do
   end do

   !**********************************************************
   !     Filtro

   do i=1,nx1
      do j=1,nx1
         do k=1,nz1-1
            varx=(9.*(varia2(i-1,j,k)+varia2(i+1,j,k))-&
               (varia2(i-3,j,k)+varia2(i+3,j,k)))/16.
            vary=(9.*(varia2(i,j-1,k)+varia2(i,j+1,k))-&
               (varia2(i,j-3,k)+varia2(i,j+3,k)))/16.
            varz=(9.*(varia2(i,j,k-1)+varia2(i,j,k+1))-&
               (varia2(i,j,k-3)+varia2(i,j,k+3)))/16.

            varia1(i,j,k)=((facx*varx+facy*vary)+facz*varz)+&
               fact*varia2(i,j,k)
         end do
      end do
   end do

   do k=1,nz1-1
      do i=1,nx1
         varia1(i,0,k)=varia1(i,1,k)
         varia1(i,nx1+1,k)=varia1(i,nx1,k)
         varia1(0,i,k)=varia1(1,i,k)
         varia1(nx1+1,i,k)=varia1(nx1,i,k)
      end do
   end do

   do i=1,nx1
      do j=1,nx1
         varia1(i,j,nz1)=varia1(i,j,nz1-1)
      end do
   end do
   !**********************************************************

   return
end