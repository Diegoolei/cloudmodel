!23456789*123456789*123456789*123456789*123456789*123456789*123456789*
!     revisado 28/04/98
!     Esta subrutina calcula la evolucion del la presion y las
!     velocidades con un paso de tiempo menor lt3
!     Las cantidades 1 son las presentes en el paso grande
!     y las 2 son las del paso futuro, las 3 son auxiliares
!     Le resta la perturbacion promedio
subroutine velpre
   USE cant01
   USE dimen
   implicit none
   include 'const.i'
   include 'perdim.i'
   include 'estbas.i'
   include 'p3v3.i'
   include 'fuvw.i'
   include 'velpre01.i'


   write(*,*) cc2(1)

   do 10 i=0,nx1+1
      do 10 j=0,nx1+1
         do 10 k=0,nz1
            U2(i,j,k)=U1(i,j,k)
            V2(i,j,k)=V1(i,j,k)
            W2(i,j,k)=W1(i,j,k)
            Pres2(i,j,k)=Pres1(i,j,k)
10 continue

   do 20 t=1,lt3
      presprom=0.
      do 30 k=1,nz1-1
         presi=-Cp*Tita0(k)*(1.+.61*Qvap0(k)/Den0(k))
         vel0=Tita0(k)*(Den0(k)+.61*Qvap0(k))
         vel1=Tita0(k-1)*(Den0(k-1)+.61*Qvap0(k-1))
         vel2=Tita0(k+1)*(Den0(k+1)+.61*Qvap0(k+1))
         vel3=cc2(k)/presi/vel0
         do 30 i=1,nx1
            do 30 j=1,nx1

               dprex=Pres2(i+1,j,k)-Pres2(i-1,j,k)
               dprey=Pres2(i,j+1,k)-Pres2(i,j-1,k)
               dprez=Pres2(i,j,k+1)-Pres2(i,j,k-1)

               presix=presi*dprex/dx2
               presiy=presi*dprey/dx2
               presiz=presi*dprez/dx2

               U3(i,j,k)=dt3*(presix+fu(i,j,k))+U2(i,j,k)
               V3(i,j,k)=dt3*(presiy+fv(i,j,k))+V2(i,j,k)
               W3(i,j,k)=dt3*(presiz+fw(i,j,k))+W2(i,j,k)

!##
!      if (i.eq.2.and.(j.eq.5.or.j.eq.28).and.k.eq.1) then
!       write(*,501) t,U3(i,j,k),dt3,presix,fu(i,j,k),U2(i,j,k)
!       write(*,501) t,V3(i,j,k),dt3,presiy,fv(i,j,k),V2(i,j,k)
!       write(*,501) t,W3(i,j,k),dt3,presiz,fw(i,j,k),W2(i,j,k)
!        write(*,501) t,U3(i,j,k),dt3,Pres2(i+1,j,k),Pres2(i-1,j,k)
!        write(*,501) t,V3(i,j,k),dt3,Pres2(i,j+1,k),Pres2(i,j-1,k)
!      endif
!##

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

!##
!      if (i.eq.16.and.(j.eq.8.or.j.eq.25).and.k.eq.1) then
!        write(*,501) t,Pres3(i,j,k),dt3,diver,fp(i,j,k),Pres2(i,j,k)
!        write(*,'(a4,3f16.8)') 'div',dvx,dvy,dvz
!        write(*,'(a5,3f16.8)') 'divz',(.8*vel2*W2(i,j,k+1)
!     &        -.8*vel1*W2(i,j,k))*2.,
!        write(*,'(4g16.8)')  vel2,W2(i,j,k+1),vel1,W2(i,j,k)
!      endif
!##


30    continue

!#
!      do 1010 i=1,nx1
!      do 1010 j=1,nx1/2
!      do 1010 k=1,nz1-1
!        if((U3(i,j,k).ne.U3(i,nx1-j+1,k)).or.(V3(i,j,k).ne.
!     &      -V3(i,nx1-j+1,k)).or.(Pres3(i,j,k).ne.
!     &      Pres3(i,nx1-j+1,k).or.W3(i,j,k).ne.
!     &      W3(i,nx1-j+1,k))) then
!           write(*,*) 'no iguales en el ciclo',t,i,j,nx1-j+1,k
!           write(*,505)  U3(i,j,k),U3(i,nx1-j+1,k)
!           write(*,505)  V3(i,j,k),V3(i,nx1-j+1,k)
!           write(*,505)  W3(i,j,k),W3(i,nx1-j+1,k)
!           write(*,505)  Pres3(i,j,k),Pres3(i,nx1-j+1,k)
!           stop
!        endif
! 1010 continue
!#




!*      redefiniciones y contornos
      do 35 i=1,nx1
         do 35 j=1,nx1
            Pres3(i,j,0)=Pres3(i,j,1)
            Pres3(i,j,nz1)=Pres3(i,j,nz1-1)
35    continue
      do 37 i=1,nx1
         do 37 k=0,nz1
            Pres3(i,0,k)=Pres3(i,1,k)
            Pres3(i,nx1+1,k)=Pres3(i,nx1,k)
            Pres3(0,i,k)=Pres3(1,i,k)
            Pres3(nx1+1,i,k)=Pres3(nx1,i,k)
37    continue

      do 40 i=1,nx1
         do 40 j=1,nx1
            do 50 k=1,nz1-1
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

!#
!      if (i.eq.16.and.(j.eq.8.or.j.eq.25).and.k.eq.1) then
!        write(*,502) 'p',t,Pres2(i,j,k),presprom
!        write(*,510) Pres3(i,j,k),Pres3(i+1,j,k),
!     &              Pres3(i-1,j,k),Pres3(i,j+1,k),Pres3(i,j-1,k),
!     &              Pres3(i,j,k+1),Pres3(i,j,k-1)
!      endif
!#


50          continue
            U2(i,j,0)=0
            V2(i,j,0)=0
            W2(i,j,0)=0
            Pres2(i,j,0)=Pres2(i,j,1)
            U2(i,j,nz1)=U2(i,j,nz1-1)
            V2(i,j,nz1)=V2(i,j,nz1-1)
            W2(i,j,nz1)=W2(i,j,nz1-1)
            Pres2(i,j,nz1)=Pres2(i,j,nz1-1)
40    continue
      do 60 i=1,nx1
         do 60 k=0,nz1
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
60    continue

      presprom=presprom/nnn
      do 70 i=0,nx1+1
         do 70 j=0,nx1+1
            do 70 k=0,nz1
               Pres2(i,j,k)=Pres2(i,j,k)-presprom

!#
!      if (i.eq.16.and.(j.eq.8.or.j.eq.25).and.k.eq.1) then
!        write(*,502) 'q',t,Pres2(i,j,k),presprom
!      endif
!#

70    continue

      if (t.eq.lt3/2) then
         do 80 i=0,nx1+1
            do 80 j=0,nx1+1
               do 80 k=0,nz1
                  U1(i,j,k)=U2(i,j,k)
                  V1(i,j,k)=V2(i,j,k)
                  W1(i,j,k)=W2(i,j,k)
                  Pres1(i,j,k)=Pres2(i,j,k)
80       continue
      endif

20 continue

!**********************************************************
!*    suavizado

   write(*,500) Pres2(9,2,2),Pres2(9,31,2)&
      ,Pres2(9,2,2),Pres2(9,31,2)


!#
!      do 1001 i=1,nx1
!      do 1001 j=1,nx1/2
!      do 1001 k=1,nz1
!        if((U2(i,j,k).ne.U2(i,nx1-j+1,k)).or.(V2(i,j,k).ne.
!     &      -V2(i,nx1-j+1,k)).or.(Pres2(i,j,k).ne.
!     &      Pres2(i,nx1-j+1,k).or.W2(i,j,k).ne.
!     &      W2(i,nx1-j+1,k))) then
!           write(*,*) 'no iguales antes del filtro',i,j,k
!           write(*,505)  U2(i,j,k),U2(i,nx1-j+1,k)
!           write(*,505)  V2(i,j,k),V2(i,nx1-j+1,k)
!           write(*,505)  W2(i,j,k),W2(i,nx1-j+1,k)
!           write(*,505)  Pres2(i,j,k),Pres2(i,nx1-j+1,k)
!           stop
!        endif
! 1001 continue
!#



   call filtro(Pres1,.15,.15,.1)

   call filtro(Pres2,.15,.15,.1)

   call filtro(U1,facx,facy,facz)
   call filtro(U2,facx,facy,facz)
   call filtro(V1,facx,facy,facz)
   call filtro(V2,facx,facy,facz)
   call filtro(W1,facx,facy,facz)
   call filtro(W2,facx,facy,facz)

   write(*,500) Pres1(9,2,1),Pres1(9,31,1)&
      ,Pres2(9,2,1),Pres2(9,31,1)

!#
!      do 1000 i=1,nx1
!      do 1000 j=1,nx1/2
!      do 1000 k=1,nz1
!        if((U2(i,j,k).ne.U2(i,nx1-j+1,k)).or.(V2(i,j,k).ne.
!     &      -V2(i,nx1-j+1,k)).or.(Pres2(i,j,k).ne.
!     &      Pres2(i,nx1-j+1,k).or.W2(i,j,k).ne.
!     &      W2(i,nx1-j+1,k))) then
!           write(*,*) 'no iguales despues del filtro',i,j,k
!           write(*,505)  U2(i,j,k),U2(i,nx1-j+1,k)
!           write(*,505)  V2(i,j,k),V2(i,nx1-j+1,k)
!           write(*,505)  W2(i,j,k),W2(i,nx1-j+1,k)
!           write(*,505)  Pres2(i,j,k),Pres2(i,nx1-j+1,k)
!           stop
!        endif
! 1000 continue
!#


500 format(4g16.8)
501 format(i3,5g16.8)
502 format(a3,i3,4g16.8)
505 format(2g16.8)
510 format(7g16.8)

   do 200 i=1,nx1
      do 200 j=1,nx1
         Pres1(i,j,0)=Pres1(i,j,1)
         Pres2(i,j,0)=Pres2(i,j,1)
200 continue
!**********************************************************

   return
end
