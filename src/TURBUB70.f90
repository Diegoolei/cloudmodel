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

subroutine turbu2(i,j,k)
   USE constants
   USE turbvar
   USE turbvar1
   implicit none
   include 'turbu2.i'
   real aux

!     calculo de KM
   do 20 lx=-1,1
      do 20 ly=-1,1
         do 20 lz=-1,1
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

20 continue



!     calculo de las derivadas (sin la distancia abajo)
   KM1=KM(1,0,0)-KM(-1,0,0)
   KM2=KM(0,1,0)-KM(0,-1,0)
   KM3=KM(0,0,1)-KM(0,0,-1)
   do 100 n=1,3
      D1(n)=D(n,1,i+1,j,2)-D(n,1,i-1,j,2)
      D2(n)=D(n,2,i,j+1,2)-D(n,2,i,j-1,2)
      D3(n)=D(n,3,i,j,3)-D(n,3,i,j,1)
100 continue
   KMM=KM(0,0,0)
   do 110 n=1,3
      do 110 m=1,3
         DD(n,m)=D(n,m,i,j,2)
110 continue

   return
end

!*********************************************************************
subroutine suma(sum,a1,a2,a3)
   implicit none
   real a1,a2,a3,sum,aux
   integer j

   do 10 j=1,2
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
10 continue

   sum=(a1+a2)+a3
   return
end
