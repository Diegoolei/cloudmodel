!     Esta subrutina calcula los Dnm para cada plano Z
subroutine turbu1(kk)
   USE cant01
   USE dimen
   USE perdim
   USE const
   USE turbvar1
   USE turbu1_vars
   implicit none

   integer, intent(in) :: kk

   k=kk+1
   do 10 i=0,nx1+1
      do 10 j=0,nx1+1

         if (kk.eq.1) then
            do 20 n=1,2
               do 20 m=1,n
                  D(n,m,i,j,1)=0.
20          continue
            D(3,1,i,j,1)=U2(i,j,1)
            D(3,2,i,j,1)=V2(i,j,1)
            D(3,3,i,j,1)=W2(i,j,1)*2./3.
            do 25 n=1,3
               do 25 m=1,n
                  D(m,n,i,j,1)=D(n,m,i,j,1)
25          continue
            do 30 lx=-1,1
               do 30 ly=-1,1
                  do 30 lz=-1,1
                     ldis=abs(lx)+abs(ly)+abs(lz)
                     if (ldis.le.1) then
                        vel(1,lx,ly,lz)=U2(lx+i,ly+j,lz+1)
                        vel(2,lx,ly,lz)=V2(lx+i,ly+j,lz+1)
                        vel(3,lx,ly,lz)=W2(lx+i,ly+j,lz+1)
                     endif
30          continue
!     calculo de Dij
            do 40 n=1,3
               dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
               dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
               dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
40          continue
            do 50 n=1,3
               do 50 m=1,n
                  D(n,m,i,j,2)=(dv(n,m)+dv(m,n))
                  D(m,n,i,j,2)=D(n,m,i,j,2)
                  if (n.eq.m) D(n,n,i,j,2)=2./3.*D(n,n,i,j,2)
50          continue
         else
            do 60 n=1,3
               do 60 m=1,3
                  do 60 lz=1,2
                     D(n,m,i,j,lz)=D(n,m,i,j,lz+1)
60          continue
         endif
!*********************************************************

!     Lectura de las velocidades necesarias
         do 100 lx=-1,1
            do 100 ly=-1,1
               do 100 lz=-1,1
                  ldis=abs(lx)+abs(ly)+abs(lz)
                  if (ldis.le.1) then
                     vel(1,lx,ly,lz)=U2(lx+i,ly+j,lz+k)
                     vel(2,lx,ly,lz)=V2(lx+i,ly+j,lz+k)
                     vel(3,lx,ly,lz)=W2(lx+i,ly+j,lz+k)
                  endif
100      continue
!     calculo de Dij
         do 130 n=1,3
            dv(n,1)=vel(n,1,0,0)-vel(n,-1,0,0)
            dv(n,2)=vel(n,0,1,0)-vel(n,0,-1,0)
            dv(n,3)=vel(n,0,0,1)-vel(n,0,0,-1)
130      continue
         do 140 n=1,3
            do 140 m=1,n
               D(n,m,i,j,3)=(dv(n,m)+dv(m,n))
               D(m,n,i,j,3)=D(n,m,i,j,3)
               if (n.eq.m) D(n,n,i,j,3)=2./3.*D(n,n,i,j,3)
140      continue


10 continue

   return
end
