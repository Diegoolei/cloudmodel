!     revisado 13/04/99
!23456789*123456789*123456789*123456789*123456789*123456789*123456789*
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

!##
!      if (l.eq.16.and.m.eq.17.and.n.ge.12.and.n.le.14) then
!        write(*,*) 'nuclcri',l,m,n,auxs,Naux,mcri,Naer
!        write(*,*) Naux,Acri,Bcri,Tc,rs
!        write(*,*) auxs,(Qvap-ess/Rv/TT*.95),Qvap
!        stop
!      endif


!        write(*,*) 'vapor?',auxs,Qvap-ess/Rv/TT

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

!      if (l.eq.16.and.(m.eq.15.or.m.eq.18).and. n.eq.10) then
!        write(*,*) 'anucl',m,e1,esl,TT,Qvap
!      endif


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
            write(*,*) 'moco en nuclea',l,m,n
            stop
         endif
      endif

      s=1
      esl=esl*exp(B*(TT2-TT1)/TT2/TT1)
      TT1=TT2
      rl=abs(e1-esl)/esl
      xxx=xxx+1
      if (rl.gt.1e-3  .and. hhh.eq.0) goto 10

      Qvap=Qvap-auxl
      Qliq=Qliq+auxl

!*     control de mocos
      if (auxl.lt.0 .or. auxs.lt.0) then
         write(*,*) 'Mocazo en nuclea',auxl,auxs,l,m,n
         stop
      endif

!*
!      variacion en los aerosoles
!      considerando que las nuevas gotitas tienen un radio Rgotmin

      Naux=Naux-auxl/(4./3.*pi*rhow*Rgotmin**3.)/1e6

      TT=TT2

   endif

!      if (auxs.gt.0.or. auxl.gt.0) then
!        write(*,*) 'nucl',l,m,n,auxs,auxl,TT
!	write(*,*) Naux,Qvap,e1,esl,ess,rl,rs
!        pause
!      endif

   return
end
