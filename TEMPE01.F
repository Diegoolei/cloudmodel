*234567890123456789012345678901234567890123456789012345678901234567890
*     revisado 28/04/98
*     Fcalo es el calor liberado por cambio de fase, por unidad de masa
*     de aire

      subroutine tempot(i,j,k,dden0z,Fcal)

      implicit none

      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'tempe01.i'

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

	lapla=(Titaa1(i+1,j,k)+Titaa1(i-1,j,k))+(Titaa1(i,j+1,k)+
     &        Titaa1(i,j-1,k))+Titaa1(i,j,k+1)+Titaa1(i,j,k-1)-
     &        6*Titaa1(i,j,k)
	lapla=lapla+(Tita0(k+1)+Tita0(k-1)-2.*Tita0(k))

	turden=dden0z*(Tita0(k+1)-Tita0(k-1))

	turbul=3.*cteturb/dx8*(escal+KMM*(4.*lapla+turden))

	Titaa2(i,j,k)=dt1*((advec+verti)/dx2+turbul+calor)+
     &               Titaa1(i,j,k)

*     control de locura
      if (abs(Titaa2(i,j,k)).gt.30) then
       write(*,*) 'tita loca',i,j,k,Titaa2(i,j,k),advec,
     &    verti,turbul,calor,Titaa1(i,j,k)
       write(*,*) adv(3),W2(i,j,k),Titaa1(i,j,k+1)
     &     ,Titaa1(i,j,k-1),Titaa1(i,j,k)
       write(*,*) adv(1),U2(i+1,j,k),Titaa1(i+1,j,k)
     &       , U2(i-1,j,k),Titaa1(i-1,j,k)
       write(*,*) adv(2),V2(i,j+1,k),Titaa1(i,j+1,k)
     &       , V2(i,j-1,k),Titaa1(i,j-1,k)
      stop
      endif
       
      return
      end
