!     Esta subrutina filtra componentes de alta frecuencia espacial.
!     El valor de la variable del punto j se filtra con los valores
!     extrapolados linalmente de los puntos j-3 y j-1 y similares,
!     pasando un polinomio de grado 4.

subroutine filtro(varia1,facx,facy,facz)
   USE dimen
   implicit none
   character*50 text
   include 'filtro01.i'

   !REAL, DIMENSION(-3:NX1+3,-3:NX1+3,-2:NZ1+2), intent(in) :: varia1
   !real, intent(in) :: facx,facy,facz
   fact=1.-(facx+facy+facz)

   if (fact.lt.0.25) then
      write(*,*) '              Moco, fact= ',fact
      text = 'El valor de la variable no debe ser menor a 0.25'
      write(*,*) text
      stop
   endif

!**********************************************************
!     Redefiniciones y contornos

   do 10 i=0,nx1+1
      do 10 j=0,nx1+1
         do 10 k=0,nz1
            varia2(i,j,k)=varia1(i,j,k)
10 continue

   do 20 k=0,nz1
      do 20 i=0,nx1
         varia2(i,-1,k)=varia2(i,1,k)
         varia2(i,-2,k)=varia2(i,1,k)
         varia2(i,nx1+2,k)=varia2(i,nx1,k)
         varia2(i,nx1+3,k)=varia2(i,nx1,k)
         varia2(-1,i,k)=varia2(1,i,k)
         varia2(-2,i,k)=varia2(1,i,k)
         varia2(nx1+2,i,k)=varia2(nx1,i,k)
         varia2(nx1+3,i,k)=varia2(nx1,i,k)
20 continue

   do 30 i=1,nx1
      do 30 j=1,nx1
         varia2(i,j,-1)=varia2(i,j,0)
         varia2(i,j,-2)=varia2(i,j,0)
         varia2(i,j,nz1+1)=varia2(i,j,nz1)
         varia2(i,j,nz1+2)=varia2(i,j,nz1)
30 continue

!**********************************************************
!     Filtro

   do 100 i=1,nx1
      do 100 j=1,nx1
         do 100 k=1,nz1-1
            varx=(9.*(varia2(i-1,j,k)+varia2(i+1,j,k))-&
               (varia2(i-3,j,k)+varia2(i+3,j,k)))/16.
            vary=(9.*(varia2(i,j-1,k)+varia2(i,j+1,k))-&
               (varia2(i,j-3,k)+varia2(i,j+3,k)))/16.
            varz=(9.*(varia2(i,j,k-1)+varia2(i,j,k+1))-&
               (varia2(i,j,k-3)+varia2(i,j,k+3)))/16.

            varia1(i,j,k)=((facx*varx+facy*vary)+facz*varz)+&
               fact*varia2(i,j,k)

!#
!      if (i.eq.9 .and.(j.eq.2.or.j.eq.31).and.k.eq.1)then
!      write(*,500) 'fil',j,varia1(i,j,k),varx,vary,varz,varia2(i,j,k)
!      write(*,501) 'zet',varia2(i,j,k-1),varia2(i,j,k+1),
!     &             varia2(i,j,k-3),varia2(i,j,k+3)
!      endif
!#


100 continue

   do 110 k=1,nz1-1
      do 110 i=1,nx1
         varia1(i,0,k)=varia1(i,1,k)
         varia1(i,nx1+1,k)=varia1(i,nx1,k)
         varia1(0,i,k)=varia1(1,i,k)
         varia1(nx1+1,i,k)=varia1(nx1,i,k)
110 continue

   do 120 i=1,nx1
      do 120 j=1,nx1
         varia1(i,j,nz1)=varia1(i,j,nz1-1)
120 continue

500 format(a4,i3,5g16.8)
501 format(a4,4g16.8)
!**********************************************************

   return
end
