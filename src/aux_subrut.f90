subroutine tempot(i,j,k,dden0z,Fcal)
   !> heat_force es el calor liberado por cambio de fase, por unidad de masa de aire
   use cant01
   use dimen
   use dinamic_var_perturbation
   use constants
   use estbas
   use turbvar
   use tempe01
   implicit none

   real, intent(in) :: dden0z,Fcal
   integer, intent(in) :: i,j,k

   adv(1)=(u_perturbed(i,j,k)+UU(k))*(thermal_property_1(i+1,j,k)-thermal_property_1(i-1,j,k))
   adv(2)=(v_perturbed(i,j,k)+VV(k))*(thermal_property_1(i,j+1,k)-thermal_property_1(i,j-1,k))
   adv(3)=w_perturbed(i,j,k)*(thermal_property_1(i,j,k+1)-thermal_property_1(i,j,k-1))

   advec=-((adv(1)+adv(2))+adv(3))

   verti=-w_perturbed(i,j,k)*(Tita0(k+1)-Tita0(k-1))

   calor=Fcal*Tita0(k)/(Temp0(k)*Cp)

   dtita(1)=thermal_property_1(i+1,j,k)-thermal_property_1(i-1,j,k)
   dtita(2)=thermal_property_1(i,j+1,k)-thermal_property_1(i,j-1,k)
   dtita(3)=thermal_property_1(i,j,k+1)-thermal_property_1(i,j,k-1)

   escal=(dtita(1)*KM1+dtita(2)*KM2)+dtita(3)*KM3

   lapla=(thermal_property_1(i+1,j,k)+thermal_property_1(i-1,j,k))+(thermal_property_1(i,j+1,k)+&
      thermal_property_1(i,j-1,k))+thermal_property_1(i,j,k+1)+thermal_property_1(i,j,k-1)-&
      6*thermal_property_1(i,j,k)
   lapla=lapla+(Tita0(k+1)+Tita0(k-1)-2.*Tita0(k))

   turden=dden0z*(Tita0(k+1)-Tita0(k-1))

   turbul=3.*cteturb/dx8*(escal+KMM*(4.*lapla+turden))

   thermal_property_2(i,j,k)=dt1*((advec+verti)/dx2+turbul+calor)+&
      thermal_property_1(i,j,k)

   !     control de locura
   if (abs(thermal_property_2(i,j,k)) > 30) then
      stop
   endif

   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de gotitas
subroutine corgot
   use dimen
   use permic
   use lmngot
   use corgot_vars
   implicit none

   neg1=0.
   pos1=0.
   do concurrent(n=ngot(1):ngot(2), l=lgot(1):lgot(2), m=mgot(1):mgot(2))
      if (Qgot2(l,m,n) < 0.) then
         neg1=neg1+Qgot2(l,m,n)
         Qgot2(l,m,n)=0
      else
         pos1=pos1+Qgot2(l,m,n)
      endif
   end do

   if(pos1 <= -neg1) then
      do concurrent(l=lgot(1):lgot(2), m=mgot(1):mgot(2), n=ngot(1):ngot(2))
         Qgot2(l,m,n)=0.
      end do
      if (-neg1 > 1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do concurrent(l=lgot(1):lgot(2), m=mgot(1):mgot(2), n=ngot(1):ngot(2))
         if(Qgot2(l,m,n) > 0) then
            Qgot2(l,m,n)=Qgot2(l,m,n)*(1.+aux1)
         endif
      end do
   endif

   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de gotas
subroutine corllu
   use dimen
   use permic
   use lmnllu
   use corgot_vars
   implicit none

   neg1=0.
   pos1=0.
   do concurrent(n=nllu(1):nllu(2), l=lllu(1):lllu(2), m=mllu(1):mllu(2))
      if (Qllu2(l,m,n) < 0.) then
         neg1=neg1+Qllu2(l,m,n)
         Qllu2(l,m,n)=0
      else
         pos1=pos1+Qllu2(l,m,n)
      endif
   end do

   if(pos1 <= -neg1) then
      do concurrent(l=lllu(1):lllu(2), m=mllu(1):mllu(2), n=nllu(1):nllu(2))
         Qllu2(l,m,n)=0.
      end do
      if (-neg1 > 1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do concurrent(l=lllu(1):lllu(2), m=mllu(1):mllu(2), n=nllu(1):nllu(2))
         if(Qllu2(l,m,n) > 0) then
            Qllu2(l,m,n)=Qllu2(l,m,n)*(1.+aux1)
         endif
      end do
   endif
   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de cristales
subroutine corcri
   use dimen
   use permic
   use lmncri
   use corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do concurrent(n=ncri(1):ncri(2), l=lcri(1):lcri(2), m=mcri(1):mcri(2))
      if (Qcri2(l,m,n) < 0.) then
         neg1=neg1+Qcri2(l,m,n)
         Qcri2(l,m,n)=0
      else
         pos1=pos1+Qcri2(l,m,n)
      endif
   end do

   if(pos1 <= -neg1) then
      do concurrent(l=lcri(1):lcri(2), m=mcri(1):mcri(2), n=ncri(1):ncri(2))
         Qcri2(l,m,n)=0.
      end do
      if (-neg1 > 1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do concurrent(l=lcri(1):lcri(2), m=mcri(1):mcri(2), n=ncri(1):ncri(2))
         if(Qcri2(l,m,n) > 0) then
            Qcri2(l,m,n)=Qcri2(l,m,n)*(1.+aux1)
         endif
      end do
   endif

   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de nieve
subroutine cornie
   use dimen
   use permic
   use lmnnie
   use corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do concurrent(n=nnie(1):nnie(2), l=lnie(1):lnie(2), m=mnie(1):mnie(2))
      if (Qnie2(l,m,n) < 0.) then
         neg1=neg1+Qnie2(l,m,n)
         Qnie2(l,m,n)=0
      else
         pos1=pos1+Qnie2(l,m,n)
      endif
   end do

   if(pos1 <= -neg1) then
      do concurrent(l=lnie(1):lnie(2), m=mnie(1):mnie(2), n=nnie(1):nnie(2))
         Qnie2(l,m,n)=0.
      end do

      if (-neg1 > 1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do concurrent(l=lnie(1):lnie(2), m=mnie(1):mnie(2), n=nnie(1):nnie(2))
         if(Qnie2(l,m,n) > 0) then
            Qnie2(l,m,n)=Qnie2(l,m,n)*(1.+aux1)
         endif
      end do
   endif
   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de granizos
subroutine corgra
   use dimen
   use permic
   use lmngra
   use corgot_vars
   implicit none
   neg1=0.
   pos1=0.
   do concurrent(n=ngra(1):ngra(2), l=lgra(1):lgra(2), m=mgra(1):mgra(2))
      if (Qgra2(l,m,n) < 0.) then
         neg1=neg1+Qgra2(l,m,n)
         Qgra2(l,m,n)=0
      else
         pos1=pos1+Qgra2(l,m,n)
      endif
   end do

   if(pos1 <= -neg1) then
      do concurrent(l=lgra(1):lgra(2), m=mgra(1):mgra(2), n=ngra(1):ngra(2))
         Qgra2(l,m,n)=0.
      end do

      if (-neg1 > 1e-3) then
         stop
      endif
   else
      aux1=neg1/pos1
      do concurrent(l=lgra(1):lgra(2), m=mgra(1):mgra(2), n=ngra(1):ngra(2))
         if(Qgra2(l,m,n) > 0) then
            Qgra2(l,m,n)=Qgra2(l,m,n)*(1.+aux1)
         endif
      end do
   endif
   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de vapor
subroutine corvap(Qvapneg)
   use dimen
   use permic
   use estbas
   use corvap_vars
   implicit none

   real(8), intent(in) :: Qvapneg

   do concurrent(k=1:nz1)
      dq=Qvapneg*Qvaprel(k)/nx1**2.
      do concurrent(i=1:nx1, j=1:nx1)
         Qvap2(i,j,k)=Qvap2(i,j,k)+dq
         if (Qvap2(i,j,k)+Qvap0(k) < 0) Qvap2(i,j,k)=-Qvap0(k)
      end do
   end do

   return
end

 !     Esta subrutina corrige los lugares en donde la dinamica da
 !     negativa la cantidad de aerosoles
subroutine coraer(aerneg)
   use dimen
   use permic
   use estbas
   use coraer_vars
   implicit none

   real(8), intent(in) :: aerneg

   do concurrent(k=1:nz1)
      dq=aerneg*aerrel(k)/nx1**2.
      do concurrent(i=1:nx1, j=1:nx1)
         aer2(i,j,k)=aer2(i,j,k)+dq
         if (aer2(i,j,k)+aer0(k) < 0) aer2(i,j,k)=-aer0(k)
      end do
   end do
   return
end

subroutine daeros(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use daeros_vars
   implicit none

   integer, intent(in) :: l,m,n

   daer(1)=aer1(l+1,m,n)-aer1(l-1,m,n)
   daer(2)=aer1(l,m+1,n)-aer1(l,m-1,n)
   daer(3)=aer1(l,m,n+1)-aer1(l,m,n-1)


   adv(1)=(((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(aer1(l+1,m,n)+aer1(l,m,n)))-&
      ((u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(aer1(l-1,m,n)+aer1(l,m,n))))/4.
   adv(1)=adv(1)+daer(1)/2.*UU(n)

   adv(2)=(((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(aer1(l,m+1,n)+aer1(l,m,n)))&
      -((v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(aer1(l,m-1,n)+aer1(l,m,n))))/4.
   adv(2)=adv(2)+daer(2)/2.*VV(n)

   advaer2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
      (aer1(l,m,n)+aer1(l,m,n+1))/4.

   adv(3)=advaer2(l,m)-advaer1(l,m)

   advec=-(adv(1)+adv(2)+adv(3))

   verti=-((w_perturbed(l,m,n+1)+w_perturbed(l,m,n))*(aer0(n+1)+aer0(n))-&
      (w_perturbed(l,m,n-1)+w_perturbed(l,m,n))*(aer0(n-1)+aer0(n)))/4.

   aux=-((u_perturbed(l+1,m,n)-u_perturbed(l-1,m,n))+(v_perturbed(l,m+1,n)-v_perturbed(l,m-1,n)))*&
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

subroutine dgotit(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dgotit_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqgot(1)=Qgot1(l+1,m,n)-Qgot1(l-1,m,n)
   dqgot(2)=Qgot1(l,m+1,n)-Qgot1(l,m-1,n)
   dqgot(3)=Qgot1(l,m,n+1)-Qgot1(l,m,n-1)

   adv(1)=((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(Qgot1(l+1,m,n)+Qgot1(l,m,n))&
      -(u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(Qgot1(l-1,m,n)+Qgot1(l,m,n)))/4.
   adv(1)=adv(1)+dqgot(1)/2.*UU(n)

   adv(2)=((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(Qgot1(l,m+1,n)+Qgot1(l,m,n))&
      -(v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(Qgot1(l,m-1,n)+Qgot1(l,m,n)))/4.
   adv(2)=adv(2)+dqgot(2)/2.*VV(n)

   advgot2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
      (Qgot1(l,m,n)+Qgot1(l,m,n+1))/4.
   if ((advgot2(l,m)-advgot1(l,m))*dt1/dx1 > Qgot1(l,m,n) .and.&
      w_perturbed(l,m,n) > 0) then
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

subroutine dvapor(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dvapor_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqvap(1)=Qvap1(l+1,m,n)-Qvap1(l-1,m,n)
   dqvap(2)=Qvap1(l,m+1,n)-Qvap1(l,m-1,n)
   dqvap(3)=Qvap1(l,m,n+1)-Qvap1(l,m,n-1)

   adv(1)=(((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*&
      (Qvap1(l+1,m,n)+Qvap1(l,m,n)))-&
      ((u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*&
      (Qvap1(l-1,m,n)+Qvap1(l,m,n))))/4.
   adv(1)=adv(1)+dqvap(1)/2.*UU(n)

   adv(2)=(((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*&
      (Qvap1(l,m+1,n)+Qvap1(l,m,n)))-&
      ((v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*&
      (Qvap1(l,m-1,n)+Qvap1(l,m,n))))/4.
   adv(2)=adv(2)+dqvap(2)/2.*VV(n)

   advvap2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
      (Qvap1(l,m,n)+Qvap1(l,m,n+1))/4.

   adv(3)=advvap2(l,m)-advvap1(l,m)

   advec=-((adv(1)+adv(2))+adv(3))

   verti=-((w_perturbed(l,m,n+1)+w_perturbed(l,m,n))*(Qvap0(n+1)+Qvap0(n))-&
      (w_perturbed(l,m,n-1)+w_perturbed(l,m,n))*(Qvap0(n-1)+Qvap0(n)))/4.

   aux=-(u_perturbed(l+1,m,n)-u_perturbed(l-1,m,n)+v_perturbed(l,m+1,n)-v_perturbed(l,m-1,n))*&
      Qvap0(n)/2.

   escal=dqvap(1)*KM1+dqvap(2)*KM2+dqvap(3)*KM3

   lapla=((Qvap1(l+1,m,n)+Qvap1(l-1,m,n))+(Qvap1(l,m+1,n)+&
      Qvap1(l,m-1,n)))+Qvap1(l,m,n+1)+Qvap1(l,m,n-1)-&
      6.*Qvap1(l,m,n)

   lapla=lapla+(Qvap0(n+1)+Qvap0(n-1)-2.*Qvap0(n))

   turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

   Qvap2(l,m,n)=dt1*((advec+verti+aux)/dx1+turbul)+Qvap1(l,m,n)

   aux=dt1/dx1
   aux=aux*(verti+advec)+dt1*turbul+Qvap1(l,m,n)

   return
end

subroutine dlluvi(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dlluvi_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqllu(1)=Qllu1(l+1,m,n)-Qllu1(l-1,m,n)
   dqllu(2)=Qllu1(l,m+1,n)-Qllu1(l,m-1,n)
   dqllu(3)=Qllu1(l,m,n+1)-Qllu1(l,m,n-1)

   adv(1)=((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(Qllu1(l+1,m,n)+Qllu1(l,m,n))&
      -(u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(Qllu1(l-1,m,n)+Qllu1(l,m,n)))/4.
   adv(1)=adv(1)+dqllu(1)/2.*UU(n)

   adv(2)=((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(Qllu1(l,m+1,n)+Qllu1(l,m,n))-&
      (v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(Qllu1(l,m-1,n)+Qllu1(l,m,n)))/4.
   adv(2)=adv(2)+dqllu(2)/2.*VV(n)

   advllu2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
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
   if (n == 1) then
      Vtllui=Av(2*n)*Rmm**.8
      Qllui=Qllu1(l,m,n)
   endif

   sedim=gam4p8/6.*(Qllus*Vtllus-Qllui*Vtllui)

   Qllu2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qllu1(l,m,n)
   return
end

subroutine dcrist(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dcrist_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqcri(1)=Qcri1(l+1,m,n)-Qcri1(l-1,m,n)
   dqcri(2)=Qcri1(l,m+1,n)-Qcri1(l,m-1,n)
   dqcri(3)=Qcri1(l,m,n+1)-Qcri1(l,m,n-1)

   adv(1)=((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(Qcri1(l+1,m,n)+Qcri1(l,m,n))&
      -(u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(Qcri1(l-1,m,n)+Qcri1(l,m,n)))/4.
   adv(1)=adv(1)+dqcri(1)/2.*UU(n)

   adv(2)=((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(Qcri1(l,m+1,n)+Qcri1(l,m,n))-&
      (v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(Qcri1(l,m-1,n)+Qcri1(l,m,n)))/4.
   adv(2)=adv(2)+dqcri(2)/2.*VV(n)

   advcri2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
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

subroutine dnieve(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dnieve_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqnie(1)=Qnie1(l+1,m,n)-Qnie1(l-1,m,n)
   dqnie(2)=Qnie1(l,m+1,n)-Qnie1(l,m-1,n)
   dqnie(3)=Qnie1(l,m,n+1)-Qnie1(l,m,n-1)

   adv(1)=((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(Qnie1(l+1,m,n)+Qnie1(l,m,n))&
      -(u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(Qnie1(l-1,m,n)+Qnie1(l,m,n)))/4.
   adv(1)=adv(1)+dqnie(1)/2.*UU(n)

   adv(2)=((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(Qnie1(l,m+1,n)+Qnie1(l,m,n))-&
      (v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(Qnie1(l,m-1,n)+Qnie1(l,m,n)))/4.
   adv(2)=adv(2)+dqnie(2)/2.*VV(n)

   advnie2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
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
   Qnie2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qnie1(l,m,n)

   return
end

subroutine dgrani(l,m,n)
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use advecs
   use turbvar
   use dgrani_vars
   implicit none

   integer, intent(in) :: l,m,n

   dqgra(1)=Qgra1(l+1,m,n)-Qgra1(l-1,m,n)
   dqgra(2)=Qgra1(l,m+1,n)-Qgra1(l,m-1,n)
   dqgra(3)=Qgra1(l,m,n+1)-Qgra1(l,m,n-1)

   adv(1)=((u_perturbed(l+1,m,n)+u_perturbed(l,m,n))*(Qgra1(l+1,m,n)+Qgra1(l,m,n))&
      -(u_perturbed(l-1,m,n)+u_perturbed(l,m,n))*(Qgra1(l-1,m,n)+Qgra1(l,m,n)))/4.
   adv(1)=adv(1)+dqgra(1)/2.*UU(n)

   adv(2)=((v_perturbed(l,m+1,n)+v_perturbed(l,m,n))*(Qgra1(l,m+1,n)+Qgra1(l,m,n))-&
      (v_perturbed(l,m-1,n)+v_perturbed(l,m,n))*(Qgra1(l,m-1,n)+Qgra1(l,m,n)))/4.
   adv(2)=adv(2)+dqgra(2)/2.*VV(n)

   advgra2(l,m)=(w_perturbed(l,m,n)+w_perturbed(l,m,n+1))*&
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
   if (n == 1) then
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
!############### calculo de la velocidad y la presion ################
subroutine speed_pressure()
   use cant01
   use dimen
   use dinamic_var_perturbation
   use constants
   use estbas
   use velpre01
   use p3v3
   use fuvw
   implicit none

   call velpre01_init()

   do concurrent (i=0:nx1+1, j=0:nx1+1, k=0:nz1)
      u_perturbed(i,j,k)=u_original(i,j,k)
      v_perturbed(i,j,k)=v_original(i,j,k)
      w_perturbed(i,j,k)=w_original(i,j,k)
      pressure_perturbed(i,j,k)=pressure_original(i,j,k)
   end do

   do concurrent(t=1:lt3)
      do concurrent(k=1:nz1-1)
         presi=-Cp*Tita0(k)*(1.+.61*Qvap0(k)/Den0(k))
         vel0=Tita0(k)*(Den0(k)+.61*Qvap0(k))
         vel1=Tita0(k-1)*(Den0(k-1)+.61*Qvap0(k-1))
         vel2=Tita0(k+1)*(Den0(k+1)+.61*Qvap0(k+1))
         vel3=cc2(k)/presi/vel0
         do concurrent (i=1:nx1, j=1:nx1)
            dprex=pressure_perturbed(i+1,j,k)-pressure_perturbed(i-1,j,k)
            dprey=pressure_perturbed(i,j+1,k)-pressure_perturbed(i,j-1,k)
            dprez=pressure_perturbed(i,j,k+1)-pressure_perturbed(i,j,k-1)

            presix=presi*dprex/dx2
            presiy=presi*dprey/dx2
            presiz=presi*dprez/dx2

            U3(i,j,k)=dt3*(presix+fu(i,j,k))+u_perturbed(i,j,k)
            V3(i,j,k)=dt3*(presiy+fv(i,j,k))+v_perturbed(i,j,k)
            W3(i,j,k)=dt3*(presiz+fw(i,j,k))+w_perturbed(i,j,k)


            dvx=vel0*(u_perturbed(i+1,j,k)-u_perturbed(i-1,j,k))
            dvy=vel0*(v_perturbed(i,j+1,k)-v_perturbed(i,j-1,k))
            if (k == 1) then
               !      dvz=tiene 80% de (w_perturbed(2)-w_perturbed(1) y 20% de (w_perturbed(1)-w_perturbed(0)
               dvz=(.8*vel2*w_perturbed(i,j,k+1)-.8*vel1*w_perturbed(i,j,k))*2.
            else
               dvz=vel2*w_perturbed(i,j,k+1)-vel1*w_perturbed(i,j,k-1)
            endif

            diver=vel3*((dvx+dvy)+dvz)/dx2

            !      modificado para agrega turbulencia en la P 23/8/97
            Pres3(i,j,k)=dt3*(diver+fp(i,j,k))+pressure_perturbed(i,j,k)
         end do
      end do

      !*      redefiniciones y contornos
      do concurrent(i=1:nx1, j=1:nx1)
         Pres3(i,j,0)=Pres3(i,j,1)
         Pres3(i,j,nz1)=Pres3(i,j,nz1-1)
      end do
      do concurrent(i=1:nx1, k=0:nz1)
         Pres3(i,0,k)=Pres3(i,1,k)
         Pres3(i,nx1+1,k)=Pres3(i,nx1,k)
         Pres3(0,i,k)=Pres3(1,i,k)
         Pres3(nx1+1,i,k)=Pres3(nx1,i,k)
      end do

      presprom=0.
      do concurrent(i=1:nx1, j=1:nx1)
         do k=1,nz1-1
            if (k == 1) then
               u_perturbed(i,j,k)=U3(i,j,k)-kkk*&
                  (2.*U3(i,j,k)-U3(i,j,k+1))
               v_perturbed(i,j,k)=V3(i,j,k)-kkk*&
                  (2.*V3(i,j,k)-V3(i,j,k+1))
               w_perturbed(i,j,k)=W3(i,j,k)-kkk*&
                  (2.*W3(i,j,k)-W3(i,j,k+1))
            else
               u_perturbed(i,j,k)=U3(i,j,k)
               v_perturbed(i,j,k)=V3(i,j,k)
               w_perturbed(i,j,k)=W3(i,j,k)
            endif
            pressure_perturbed(i,j,k)=prom1*Pres3(i,j,k)+prom*(&
               ((Pres3(i+1,j,k)+ Pres3(i-1,j,k))+&
               (Pres3(i,j+1,k)+Pres3(i,j-1,k)))+&
               Pres3(i,j,k+1)+Pres3(i,j,k-1))
            presprom=pressure_perturbed(i,j,k)+presprom
         end do

         u_perturbed(i,j,0)=0
         v_perturbed(i,j,0)=0
         w_perturbed(i,j,0)=0
         pressure_perturbed(i,j,0)=pressure_perturbed(i,j,1)
         u_perturbed(i,j,nz1)=u_perturbed(i,j,nz1-1)
         v_perturbed(i,j,nz1)=v_perturbed(i,j,nz1-1)
         w_perturbed(i,j,nz1)=w_perturbed(i,j,nz1-1)
         pressure_perturbed(i,j,nz1)=pressure_perturbed(i,j,nz1-1)
      end do
      do concurrent(i=1:nx1, k=0:nz1)
         u_perturbed(0,i,k)=u_perturbed(1,i,k)
         v_perturbed(0,i,k)=v_perturbed(1,i,k)
         w_perturbed(0,i,k)=w_perturbed(1,i,k)
         pressure_perturbed(0,i,k)=pressure_perturbed(1,i,k)
         u_perturbed(nx1+1,i,k)=u_perturbed(nx1,i,k)
         v_perturbed(nx1+1,i,k)=v_perturbed(nx1,i,k)
         w_perturbed(nx1+1,i,k)=w_perturbed(nx1,i,k)
         pressure_perturbed(nx1+1,i,k)=pressure_perturbed(nx1,i,k)
         u_perturbed(i,0,k)=u_perturbed(i,1,k)
         v_perturbed(i,0,k)=v_perturbed(i,1,k)
         w_perturbed(i,0,k)=w_perturbed(i,1,k)
         pressure_perturbed(i,0,k)=pressure_perturbed(i,1,k)
         u_perturbed(i,nx1+1,k)=u_perturbed(i,nx1,k)
         v_perturbed(i,nx1+1,k)=v_perturbed(i,nx1,k)
         w_perturbed(i,nx1+1,k)=w_perturbed(i,nx1,k)
         pressure_perturbed(i,nx1+1,k)=pressure_perturbed(i,nx1,k)
      end do

      presprom=presprom/nnn
      do concurrent(i=0:nx1+1, j=0:nx1+1, k=0:nz1)
         pressure_perturbed(i,j,k)=pressure_perturbed(i,j,k)-presprom
      end do

      if (t == lt3/2) then
         do concurrent(i=0:nx1+1, j=0:nx1+1, k=0:nz1)
            u_original(i,j,k)=u_perturbed(i,j,k)
            v_original(i,j,k)=v_perturbed(i,j,k)
            w_original(i,j,k)=w_perturbed(i,j,k)
            pressure_original(i,j,k)=pressure_perturbed(i,j,k)
         end do
      endif

   end do

   !**********************************************************
   !*    suavizado

   call filtro(pressure_original,.15,.15,.1)

   call filtro(pressure_perturbed,.15,.15,.1)

   call filtro(u_original,facx,facy,facz)
   call filtro(u_perturbed,facx,facy,facz)
   call filtro(v_original,facx,facy,facz)
   call filtro(v_perturbed,facx,facy,facz)
   call filtro(w_original,facx,facy,facz)
   call filtro(w_perturbed,facx,facy,facz)

   do concurrent(i=1:nx1, j=1:nx1)
      pressure_original(i,j,0)=pressure_original(i,j,1)
      pressure_perturbed(i,j,0)=pressure_perturbed(i,j,1)
   end do
   !**********************************************************

   return
end subroutine speed_pressure

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !> TURBUA01
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine turbu1(kk)
   !> Esta subrutina calcula los Dnm para cada plano Z
   use cant01
   use dimen
   use dinamic_var_perturbation
   use constants
   use turbvar1
   use turbu1_vars
   implicit none

   integer, intent(in) :: kk

   k=kk+1
   do concurrent(i=0:nx1+1, j=0:nx1+1)
      if (kk == 1) then
         do concurrent(n=1:2)
            do concurrent(m=1:n)
               D(n,m,i,j,1)=0.
            end do
         end do
         D(3,1,i,j,1)=u_perturbed(i,j,1)
         D(3,2,i,j,1)=v_perturbed(i,j,1)
         D(3,3,i,j,1)=w_perturbed(i,j,1)*2./3.
         do concurrent(n=1:3)
            do concurrent(m=1:n)
               D(m,n,i,j,1)=D(n,m,i,j,1)
            end do
         end do
         do concurrent(lx=-1:1, ly=-1:1, lz=-1:1)
            ldis=abs(lx)+abs(ly)+abs(lz)
            if (ldis <= 1) then
               vel(1,lx,ly,lz)=u_perturbed(lx+i,ly+j,lz+1)
               vel(2,lx,ly,lz)=v_perturbed(lx+i,ly+j,lz+1)
               vel(3,lx,ly,lz)=w_perturbed(lx+i,ly+j,lz+1)
            endif
         end do
         !     calculo de Dij
         do concurrent(n=1:3)
            dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
            dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
            dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
         end do
         do concurrent(n=1:3)
            do concurrent(m=1:n)
               D(n,m,i,j,2)=(dv(n,m)+dv(m,n))
               D(m,n,i,j,2)=D(n,m,i,j,2)
               if (n == m) D(n,n,i,j,2)=2./3.*D(n,n,i,j,2)
            end do
         end do
      else
         do concurrent(n=1:3, m=1:3, lz=1:2)
            D(n,m,i,j,lz)=D(n,m,i,j,lz+1)
         end do
      endif
      !*********************************************************

      !     Lectura de las velocidades necesarias
      do concurrent(lx=-1:1, ly=-1:1, lz=-1:1)
         ldis=abs(lx)+abs(ly)+abs(lz)
         if (ldis <= 1) then
            vel(1,lx,ly,lz)=u_perturbed(lx+i,ly+j,lz+k)
            vel(2,lx,ly,lz)=v_perturbed(lx+i,ly+j,lz+k)
            vel(3,lx,ly,lz)=w_perturbed(lx+i,ly+j,lz+k)
         endif
      end do
      !     calculo de Dij
      do concurrent(n=1:3)
         dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
         dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
         dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
      end do
      do concurrent(n=1:3)
         do concurrent(m=1:n)
            D(n,m,i,j,3)=(dv(n,m)+dv(m,n))
            D(m,n,i,j,3)=D(n,m,i,j,3)
            if (n == m) D(n,n,i,j,3)=2./3.*D(n,n,i,j,3)
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
   use dimen
   use turbvar
   use turbvar1
   use turbu2_vars
   implicit none
   integer, intent(in) :: i,j
   real aux

   !     calculo de KM
   do lx=-1,1
      do ly=-1,1
         do lz=-1,1
            ldis=abs(lx)+abs(ly)+abs(lz)
            if (ldis <= 2) then
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
   do concurrent(n=1:3)
      D1(n)=D(n,1,i+1,j,2)-D(n,1,i-1,j,2)
      D2(n)=D(n,2,i,j+1,2)-D(n,2,i,j-1,2)
      D3(n)=D(n,3,i,j,3)-D(n,3,i,j,1)
   end do
   KMM=KM(0,0,0)
   do concurrent(n=1:3, m=1:3)
      DD(n,m)=D(n,m,i,j,2)
   end do

   return
end

 !*********************************************************************
subroutine suma(sum,a1,a2,a3)
   implicit none
   real a1,a2,a3,sum,aux
   integer j

   do concurrent(j=1:2)
      if (a1 > a2) then
         aux=a1
         a1=a2
         a2=aux
      endif
      if (a2 > a3) then
         aux=a2
         a2=a3
         a3=aux
      endif
   end do

   sum=(a1+a2)+a3
   return
end

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! > NUCLEA91
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine nuclea(Qvap,Qliq,Naer,TT,rhoa,e1,esl,ess,rl,rs,Lvl,Lvs,Naux,auxl,auxs)
   use cant01
   use dimen
   use constants
   use nuclea61
   implicit none

   real, intent(in) ::  Naer, rhoa, rs, Lvl, Lvs
   real, intent(inout) :: Qliq, Qvap, TT, e1, esl, ess, rl
   real, intent(inout) :: Naux, auxl, auxs   !     Numero de aesosoles

   B=Lvl/Rv
   auxl=0.
   auxs=0.
   Naux=0.

   Rcri=5e-5
   if(TT < T0) Rcri=Rcri-4e-5*(T0-TT)/40.
   if(Rcri < 1e-5) Rcri=1e-5

   !     nucleacion sobre cristales
   Tc=T0-TT
   if (Tc  >  0 .and. rs  >  0 .and. Naer > 0) then

      mcri=pi*Rcri**3./10.*rhocri

      Naux=Acri*exp(Bcri*Tc)
      if (Naux  >  .9) Naux=.9
      Naux=Naux*Naer*1e6   ! en m3
      auxs=Naux*mcri

      if (auxs  >  (Qvap-ess/Rv/TT*.95)) then
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
   if (e1 > esl) then
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

      if (Qliq1 < 0) then
         Qliq1=0.
         auxl=-Qliq
         TT2=esl/(e1/TT1-auxl*Rv)
         e1=(Qvap-auxl)*Rv*TT2
         hhh=1
         if (s == 1) then
            stop
         endif
      endif

      s=1
      esl=esl*exp(B*(TT2-TT1)/TT2/TT1)
      TT1=TT2
      rl=abs(e1-esl)/esl
      xxx=xxx+1
      if (rl > 1e-3 .and. hhh == 0) goto 10

      Qvap=Qvap-auxl
      Qliq=Qliq+auxl

      !*     control de mocos
      if (auxl < 0 .or. auxs < 0) then
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
   use cant01
   use dimen
   use dinamic_var_perturbation
   use permic
   use constants
   use estbas
   use turbvar
   use fuvw
   use inomo_var
   implicit none

   integer, intent(in) :: i,j,k
   real, intent(in) :: dden0z

   dvelxx=(u_perturbed(i-2,j,k)-u_perturbed(i+2,j,k))+8.*(u_perturbed(i+1,j,k)-u_perturbed(i-1,j,k))
   dvelxy=(u_perturbed(i,j-2,k)-u_perturbed(i,j+2,k))+8.*(u_perturbed(i,j+1,k)-u_perturbed(i,j-1,k))
   dvelxz=(u_perturbed(i,j,k-2)-u_perturbed(i,j,k+2))+8.*(u_perturbed(i,j,k+1)-u_perturbed(i,j,k-1))
   dvelyx=(v_perturbed(i-2,j,k)-v_perturbed(i+2,j,k))+8.*(v_perturbed(i+1,j,k)-v_perturbed(i-1,j,k))
   dvelyy=(v_perturbed(i,j-2,k)-v_perturbed(i,j+2,k))+8.*(v_perturbed(i,j+1,k)-v_perturbed(i,j-1,k))
   dvelyz=(v_perturbed(i,j,k-2)-v_perturbed(i,j,k+2))+8.*(v_perturbed(i,j,k+1)-v_perturbed(i,j,k-1))
   dvelzx=(w_perturbed(i-2,j,k)-w_perturbed(i+2,j,k))+8.*(w_perturbed(i+1,j,k)-w_perturbed(i-1,j,k))
   dvelzy=(w_perturbed(i,j-2,k)-w_perturbed(i,j+2,k))+8.*(w_perturbed(i,j+1,k)-w_perturbed(i,j-1,k))
   dvelzz=(w_perturbed(i,j,k-2)-w_perturbed(i,j,k+2))+8.*(w_perturbed(i,j,k+1)-w_perturbed(i,j,k-1))

   diverx=(u_perturbed(i,j,k)*dvelxx+v_perturbed(i,j,k)*dvelxy)+w_perturbed(i,j,k)*dvelxz
   divery=(u_perturbed(i,j,k)*dvelyx+v_perturbed(i,j,k)*dvelyy)+w_perturbed(i,j,k)*dvelyz
   diverz=(u_perturbed(i,j,k)*dvelzx+v_perturbed(i,j,k)*dvelzy)+w_perturbed(i,j,k)*dvelzz

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
   grave=G*(thermal_property_1(i,j,k)/Tita0(k)+(AA*Qvap1(i,j,k)-&
      Qgot1(i,j,k)-Qllu1(i,j,k)-Qcri1(i,j,k)-&
      Qnie1(i,j,k)-Qgra1(i,j,k))/Den0(k))
   !      grave=G*(thermal_property_1(i,j,k)/Tita0(k))

   fu(i,j,k)=turbulx/dx8-diverx/dx12
   fv(i,j,k)=turbuly/dx8-divery/dx12
   fw(i,j,k)=turbulz/dx8-diverz/dx12+grave

   !     agregado para la P (23/8/97)
   laplap=(pressure_perturbed(i+1,j,k)+pressure_perturbed(i-1,j,k))+(pressure_perturbed(i,j+1,k)+&
      pressure_perturbed(i,j-1,k))+pressure_perturbed(i,j,k+1)+pressure_perturbed(i,j,k-1)-&
      6.*pressure_perturbed(i,j,k)
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
   !filtro para thermal_property_1 Qvap1
   use dimen
   use filtro01
   implicit none
   character*50 text
   REAL, DIMENSION(-3:NX1+3,-3:NX1+3,-2:NZ1+2), intent(inout) :: varia1
   real, intent(in) :: facx,facy,facz
   fact=1.-(facx+facy+facz)

   if (fact < 0.25) then
      stop
   endif

   !**********************************************************
   !     Redefiniciones y contornos

   do concurrent(i=0:nx1+1, j=0:nx1+1, k=0:nz1)
      varia2(i,j,k)=varia1(i,j,k)
   end do

   do concurrent(k=0:nz1, i=0:nx1)
      varia2(i,-1,k)=varia2(i,1,k)
      varia2(i,-2,k)=varia2(i,1,k)
      varia2(i,nx1+2,k)=varia2(i,nx1,k)
      varia2(i,nx1+3,k)=varia2(i,nx1,k)
      varia2(-1,i,k)=varia2(1,i,k)
      varia2(-2,i,k)=varia2(1,i,k)
      varia2(nx1+2,i,k)=varia2(nx1,i,k)
      varia2(nx1+3,i,k)=varia2(nx1,i,k)
   end do

   do concurrent(i=1:nx1, j=1:nx1)
      varia2(i,j,-1)=varia2(i,j,0)
      varia2(i,j,-2)=varia2(i,j,0)
      varia2(i,j,nz1+1)=varia2(i,j,nz1)
      varia2(i,j,nz1+2)=varia2(i,j,nz1)
   end do

   !**********************************************************
   !     Filtro

   do concurrent(i=1:nx1, j=1:nx1, k=1:nz1-1)
      varx=(9.*(varia2(i-1,j,k)+varia2(i+1,j,k))-&
         (varia2(i-3,j,k)+varia2(i+3,j,k)))/16.
      vary=(9.*(varia2(i,j-1,k)+varia2(i,j+1,k))-&
         (varia2(i,j-3,k)+varia2(i,j+3,k)))/16.
      varz=(9.*(varia2(i,j,k-1)+varia2(i,j,k+1))-&
         (varia2(i,j,k-3)+varia2(i,j,k+3)))/16.

      varia1(i,j,k)=((facx*varx+facy*vary)+facz*varz)+&
         fact*varia2(i,j,k)
   end do

   do concurrent(k=1:nz1-1, i=1:nx1)
      varia1(i,0,k)=varia1(i,1,k)
      varia1(i,nx1+1,k)=varia1(i,nx1,k)
      varia1(0,i,k)=varia1(1,i,k)
      varia1(nx1+1,i,k)=varia1(nx1,i,k)
   end do

   do concurrent(i=1:nx1, j=1:nx1)
      varia1(i,j,nz1)=varia1(i,j,nz1-1)
   end do
   !**********************************************************
   return
end
