*     Revision 28/04/98
      subroutine daeros(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'cant01.i'
      include 'advecs.i'
      include 'turbvar.i'
      include 'daeros.i'


      daer(1)=aer1(l+1,m,n)-aer1(l-1,m,n)
      daer(2)=aer1(l,m+1,n)-aer1(l,m-1,n)
      daer(3)=aer1(l,m,n+1)-aer1(l,m,n-1)


      adv(1)=(((U2(l+1,m,n)+U2(l,m,n))*(aer1(l+1,m,n)+aer1(l,m,n)))-
     &      ((U2(l-1,m,n)+U2(l,m,n))*(aer1(l-1,m,n)+aer1(l,m,n))))/4.
      adv(1)=adv(1)+daer(1)/2.*UU(n)

      adv(2)=(((V2(l,m+1,n)+V2(l,m,n))*(aer1(l,m+1,n)+aer1(l,m,n)))
     &      -((V2(l,m-1,n)+V2(l,m,n))*(aer1(l,m-1,n)+aer1(l,m,n))))/4.
      adv(2)=adv(2)+daer(2)/2.*VV(n)

      advaer2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (aer1(l,m,n)+aer1(l,m,n+1))/4.

      adv(3)=advaer2(l,m)-advaer1(l,m)
	    
      advec=-(adv(1)+adv(2)+adv(3))
      
      verti=-((W2(l,m,n+1)+W2(l,m,n))*(aer0(n+1)+aer0(n))-
     &      (W2(l,m,n-1)+W2(l,m,n))*(aer0(n-1)+aer0(n)))/4.

*## agregado
       aux=-((U2(l+1,m,n)-U2(l-1,m,n))+(V2(l,m+1,n)-V2(l,m-1,n)))*
     &     aer0(n)/2.
      
      escal=daer(1)*KM1+daer(2)*KM2+daer(3)*KM3    

      lapla=aer1(l+1,m,n)+aer1(l,m+1,n)+aer1(l,m,n+1)+
     &      aer1(l-1,m,n)+aer1(l,m-1,n)+aer1(l,m,n-1)-
     &      6.*aer1(l,m,n)

      lapla=((aer1(l+1,m,n)+aer1(l-1,m,n))+(aer1(l,m+1,n)+
     &      aer1(l,m-1,n)))+aer1(l,m,n-1)+aer1(l,m,n+1)-
     &      6.*aer1(l,m,n)


      lapla=lapla+(aer0(n+1)+aer0(n-1)-2.*aer0(n))

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)
      
      aer2(l,m,n)=dt1*((advec+verti+aux)/dx1+turbul)+aer1(l,m,n)

*      if(((l.eq.11.and.m.eq.3).or.(l.eq.3.and.m.eq.11)).and.n.eq.4)then
*	write(*,'(a5,i3,4g16.8)')'daer',l,aer2(l,m,n),aer1(l,m,n),aux
*	write(*,'(5g16.8)') turbul,advec,verti,adv(1),adv(2),adv(3)
*	write(*,'(4g16.8)') turbul,escal,KMM,lapla
*	write(*,'(2g16.8)') advaer2(l,m),advaer1(l,m)
*	write(*,'(3g16.8)') W2(l,m,n-1),W2(l,m,n),W2(l,m,n+1)
*	write(*,'(3g16.8)') advaer2(l,m),aer1(l,m,n-1),aer1(l,m,n+1)
*        stop
*      endif


      return
      end
      
*********************************************************************
*     revisada 28/04/98
      subroutine dgotit(l,m,n)
      
      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'perdim.i'
      include 'permic.i'
      include 'estbas.i'
      include 'advecs.i'      
      include 'cant01.i'      
      include 'turbvar.i'
      include 'dgotit.i'      
     
      dqgot(1)=Qgot1(l+1,m,n)-Qgot1(l-1,m,n)
      dqgot(2)=Qgot1(l,m+1,n)-Qgot1(l,m-1,n)
      dqgot(3)=Qgot1(l,m,n+1)-Qgot1(l,m,n-1)
      
      adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qgot1(l+1,m,n)+Qgot1(l,m,n))
     &      -(U2(l-1,m,n)+U2(l,m,n))*(Qgot1(l-1,m,n)+Qgot1(l,m,n)))/4.
      adv(1)=adv(1)+dqgot(1)/2.*UU(n)

      adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qgot1(l,m+1,n)+Qgot1(l,m,n))
     &      -(V2(l,m-1,n)+V2(l,m,n))*(Qgot1(l,m-1,n)+Qgot1(l,m,n)))/4.
      adv(2)=adv(2)+dqgot(2)/2.*VV(n)

      advgot2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qgot1(l,m,n)+Qgot1(l,m,n+1))/4.
      if ((advgot2(l,m)-advgot1(l,m))*dt1/dx1.gt.Qgot1(l,m,n) .and.
     &     W2(l,m,n).gt.0) then
	advgot2(l,m)=advgot1(l,m)+Qgot1(l,m,n)*dx1/dt1
      endif     
      adv(3)=advgot2(l,m)-advgot1(l,m)
      

      advec=-(adv(1)+adv(2)+adv(3))/dx1
      
      escal=dqgot(1)*KM1+dqgot(2)*KM2+dqgot(3)*KM3

   
      lapla=Qgot1(l+1,m,n)+Qgot1(l,m+1,n)+Qgot1(l,m,n+1)+
     &      Qgot1(l-1,m,n)+Qgot1(l,m-1,n)+Qgot1(l,m,n-1)-
     &      6.*Qgot1(l,m,n)
      
      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)
      
      Qgot2(l,m,n)=dt1*(advec+turbul)+Qgot1(l,m,n)

*##
      if(l.eq.13 .and. m.eq.16 .and. n.eq.8) then
	write(*,*) 'qgot',advgot2(l,m),advgot1(l,m),Qgot2(l,m,n)
     &             ,Qgot1(l,m,n)
      endif

      return
      end
*********************************************************************
*     Revision 28/04/98
      subroutine dvapor(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'advecs.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'dvapor.i'

      dqvap(1)=Qvap1(l+1,m,n)-Qvap1(l-1,m,n)
      dqvap(2)=Qvap1(l,m+1,n)-Qvap1(l,m-1,n)
      dqvap(3)=Qvap1(l,m,n+1)-Qvap1(l,m,n-1)

      adv(1)=(((U2(l+1,m,n)+U2(l,m,n))*(Qvap1(l+1,m,n)+Qvap1(l,m,n)))-
     &      ((U2(l-1,m,n)+U2(l,m,n))*(Qvap1(l-1,m,n)+Qvap1(l,m,n))))/4.
      adv(1)=adv(1)+dqvap(1)/2.*UU(n)

      adv(2)=(((V2(l,m+1,n)+V2(l,m,n))*(Qvap1(l,m+1,n)+Qvap1(l,m,n)))-
     &      ((V2(l,m-1,n)+V2(l,m,n))*(Qvap1(l,m-1,n)+Qvap1(l,m,n))))/4.
      adv(2)=adv(2)+dqvap(2)/2.*VV(n)


      advvap2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qvap1(l,m,n)+Qvap1(l,m,n+1))/4.

      adv(3)=advvap2(l,m)-advvap1(l,m)

      advec=-((adv(1)+adv(2))+adv(3))
      
      verti=-((W2(l,m,n+1)+W2(l,m,n))*(Qvap0(n+1)+Qvap0(n))-
     &      (W2(l,m,n-1)+W2(l,m,n))*(Qvap0(n-1)+Qvap0(n)))/4.

*## agregado
       aux=-(U2(l+1,m,n)-U2(l-1,m,n)+V2(l,m+1,n)-V2(l,m-1,n))*
     &     Qvap0(n)/2.
      
      escal=dqvap(1)*KM1+dqvap(2)*KM2+dqvap(3)*KM3    

      lapla=((Qvap1(l+1,m,n)+Qvap1(l-1,m,n))+(Qvap1(l,m+1,n)+
     &      Qvap1(l,m-1,n)))+Qvap1(l,m,n+1)+Qvap1(l,m,n-1)-
     &      6.*Qvap1(l,m,n)

      lapla=lapla+(Qvap0(n+1)+Qvap0(n-1)-2.*Qvap0(n))

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)
      
*      Qvap2(l,m,n)=dt1*((advec+verti)/dx1+turbul)+Qvap1(l,m,n)
      Qvap2(l,m,n)=dt1*((advec+verti+aux)/dx1+turbul)+Qvap1(l,m,n)

      aux=dt1/dx1
      aux=aux*(verti+advec)+dt1*turbul+Qvap1(l,m,n)

*      if(l.eq.14 .and. m.eq.15 .and. n.eq.21) then
*	write(*,*) 'dvapo',Qvap2(l,m,n),Qvap1(l,m,n),aux,turbul
*	write(*,*) advec,verti,adv(1),adv(2),adv(3),advvap2(l,m)
*	write(*,*) advvap1(l,m),W2(l,m,n-1),W2(l,m,n),W2(l,m,n+1)
*	write(*,*) advvap2(l,m),Qvap1(l,m,n-1),Qvap1(l,m,n+1)
*	write(*,*) Qvap0(n-1),Qvap0(n),Qvap0(n+1)
*        stop
*      endif

      return
      end
***********************************************************
*     Revision 28/04/98
      subroutine dlluvi(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'advecs.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'dlluvi.i'

      dqllu(1)=Qllu1(l+1,m,n)-Qllu1(l-1,m,n)
      dqllu(2)=Qllu1(l,m+1,n)-Qllu1(l,m-1,n)
      dqllu(3)=Qllu1(l,m,n+1)-Qllu1(l,m,n-1)

      adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qllu1(l+1,m,n)+Qllu1(l,m,n))
     &      -(U2(l-1,m,n)+U2(l,m,n))*(Qllu1(l-1,m,n)+Qllu1(l,m,n)))/4.
      adv(1)=adv(1)+dqllu(1)/2.*UU(n)

      adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qllu1(l,m+1,n)+Qllu1(l,m,n))-
     &     (V2(l,m-1,n)+V2(l,m,n))*(Qllu1(l,m-1,n)+Qllu1(l,m,n)))/4.
      adv(2)=adv(2)+dqllu(2)/2.*VV(n)

      advllu2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qllu1(l,m,n)+Qllu1(l,m,n+1))/4.

      adv(3)=advllu2(l,m)-advllu1(l,m)

      advec=-(adv(1)+adv(2)+adv(3))

      escal=dqllu(1)*KM1+dqllu(2)*KM2+dqllu(3)*KM3    

      lapla=Qllu1(l+1,m,n)+Qllu1(l,m+1,n)+Qllu1(l,m,n+1)+
     &      Qllu1(l-1,m,n)+Qllu1(l,m-1,n)+Qllu1(l,m,n-1)-
     &      6.*Qllu1(l,m,n)

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

****  termino de sedimentacion

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
****
      
      Qllu2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qllu1(l,m,n)

*##
*      if (l.eq.16 .and. m.eq.16 .and. n.gt.4 .and. n.le. 7) then
*        write(*,*) 'dina',Qllu2(l,m,n),Qllu1(l,m,n),advec,sedim,turbul
*        write(*,*) gam4p8,Qllus,Vtllus,Qllui,Vtllui
*       write(*,*) Av(2*n+1),Rms,Av(2*n),Rmm
*       write(*,*) Qllu1(l,m,n),Qllu1(l,m,n+1),cteqllu
*      endif


      return
      end

***********************************************************
*     Revision 29/12/98
      subroutine dcrist(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'advecs.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'dcrist.i'

      dqcri(1)=Qcri1(l+1,m,n)-Qcri1(l-1,m,n)
      dqcri(2)=Qcri1(l,m+1,n)-Qcri1(l,m-1,n)
      dqcri(3)=Qcri1(l,m,n+1)-Qcri1(l,m,n-1)

      adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qcri1(l+1,m,n)+Qcri1(l,m,n))
     &      -(U2(l-1,m,n)+U2(l,m,n))*(Qcri1(l-1,m,n)+Qcri1(l,m,n)))/4.
      adv(1)=adv(1)+dqcri(1)/2.*UU(n)

      adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qcri1(l,m+1,n)+Qcri1(l,m,n))-
     &     (V2(l,m-1,n)+V2(l,m,n))*(Qcri1(l,m-1,n)+Qcri1(l,m,n)))/4.
      adv(2)=adv(2)+dqcri(2)/2.*VV(n)

      advcri2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qcri1(l,m,n)+Qcri1(l,m,n+1))/4.

      adv(3)=advcri2(l,m)-advcri1(l,m)

      advec=-(adv(1)+adv(2)+adv(3))/dx1

      escal=dqcri(1)*KM1+dqcri(2)*KM2+dqcri(3)*KM3    

      lapla=Qcri1(l+1,m,n)+Qcri1(l,m+1,n)+Qcri1(l,m,n+1)+
     &      Qcri1(l-1,m,n)+Qcri1(l,m-1,n)+Qcri1(l,m,n-1)-
     &      6.*Qcri1(l,m,n)

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)
      
      Qcri2(l,m,n)=dt1*(advec+turbul)+Qcri1(l,m,n)


*##
*      if(l.eq.17.and.m.eq.19.and.n.eq.17) then
*        write(*,*) 'crisd',Qcri1(l,m,n),Qcri2(l,m,n),turbul
*        write(*,*) advec,adv(1),adv(2),adv(3)
*        write(*,*) W2(l,m,n-1),W2(l,m,n),W2(l,m,n+1)
*        write(*,*) Qcri1(l,m,n-1),Qcri1(l,m,n),Qcri1(l,m,n+1)	
*      endif

      return
      end
***********************************************************
*$$
*     Revision 7/06/99
      subroutine dnieve(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'advecs.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'dnieve.i'

      dqnie(1)=Qnie1(l+1,m,n)-Qnie1(l-1,m,n)
      dqnie(2)=Qnie1(l,m+1,n)-Qnie1(l,m-1,n)
      dqnie(3)=Qnie1(l,m,n+1)-Qnie1(l,m,n-1)

      adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qnie1(l+1,m,n)+Qnie1(l,m,n))
     &      -(U2(l-1,m,n)+U2(l,m,n))*(Qnie1(l-1,m,n)+Qnie1(l,m,n)))/4.
      adv(1)=adv(1)+dqnie(1)/2.*UU(n)

      adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qnie1(l,m+1,n)+Qnie1(l,m,n))-
     &     (V2(l,m-1,n)+V2(l,m,n))*(Qnie1(l,m-1,n)+Qnie1(l,m,n)))/4.
      adv(2)=adv(2)+dqnie(2)/2.*VV(n)

      advnie2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qnie1(l,m,n)+Qnie1(l,m,n+1))/4.

      adv(3)=advnie2(l,m)-advnie1(l,m)

      advec=-(adv(1)+adv(2)+adv(3))

      escal=dqnie(1)*KM1+dqnie(2)*KM2+dqnie(3)*KM3    

      lapla=Qnie1(l+1,m,n)+Qnie1(l,m+1,n)+Qnie1(l,m,n+1)+
     &      Qnie1(l-1,m,n)+Qnie1(l,m-1,n)+Qnie1(l,m,n-1)-
     &      6.*Qnie1(l,m,n)

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

****  termino de sedimentacion

      Qnies=(Qnie1(l,m,n+1)+Qnie1(l,m,n))/2.
      Qniei=(Qnie1(l,m,n-1)+Qnie1(l,m,n))/2.

      sedim=Vtnie(2*n+1)*Qnies-Vtnie(2*n-1)*Qniei
****
      
      Qnie2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qnie1(l,m,n)


*##
*      if(l.eq.32.and.m.eq.34.and.n.eq.1) then
*        write(*,*) 'nied',Qnie1(l,m,n),Qnie2(l,m,n),turbul
*        write(*,*) sedim,advec,adv(1),adv(2),adv(3)
*        write(*,*) Qnies,Qniei,Qnie1(l,m,n-1),Qnie1(l,m,n+1)
*      endif

      return
      end
***********************************************************
*     Revision 1/02/99
      subroutine dgrani(l,m,n)

      implicit none
      include 'dimen.i'
      include 'const.i'
      include 'estbas.i'
      include 'perdim.i'
      include 'permic.i'
      include 'advecs.i'
      include 'cant01.i'
      include 'turbvar.i'
      include 'dgrani.i'

      dqgra(1)=Qgra1(l+1,m,n)-Qgra1(l-1,m,n)
      dqgra(2)=Qgra1(l,m+1,n)-Qgra1(l,m-1,n)
      dqgra(3)=Qgra1(l,m,n+1)-Qgra1(l,m,n-1)

      adv(1)=((U2(l+1,m,n)+U2(l,m,n))*(Qgra1(l+1,m,n)+Qgra1(l,m,n))
     &      -(U2(l-1,m,n)+U2(l,m,n))*(Qgra1(l-1,m,n)+Qgra1(l,m,n)))/4.
      adv(1)=adv(1)+dqgra(1)/2.*UU(n)

      adv(2)=((V2(l,m+1,n)+V2(l,m,n))*(Qgra1(l,m+1,n)+Qgra1(l,m,n))-
     &     (V2(l,m-1,n)+V2(l,m,n))*(Qgra1(l,m-1,n)+Qgra1(l,m,n)))/4.
      adv(2)=adv(2)+dqgra(2)/2.*VV(n)

      advgra2(l,m)=(W2(l,m,n)+W2(l,m,n+1))*
     &             (Qgra1(l,m,n)+Qgra1(l,m,n+1))/4.

      adv(3)=advgra2(l,m)-advgra1(l,m)

      advec=-(adv(1)+adv(2)+adv(3))

      escal=dqgra(1)*KM1+dqgra(2)*KM2+dqgra(3)*KM3    

      lapla=Qgra1(l+1,m,n)+Qgra1(l,m+1,n)+Qgra1(l,m,n+1)+
     &      Qgra1(l-1,m,n)+Qgra1(l,m-1,n)+Qgra1(l,m,n-1)-
     &      6.*Qgra1(l,m,n)

      turbul=cteturb*(escal/dx8+KMM/dx2*lapla)

****  termino de sedimentacion
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
****


      Qgra2(l,m,n)=dt1*((advec+sedim)/dx1+turbul)+Qgra1(l,m,n)


      return
      end
