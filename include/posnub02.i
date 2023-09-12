!     desplazamientos horizontales a partir de la velocidad media de la nube
!     calculo la altura media de la nube, la velocidad media de la nube
!     es tomada como la velocidad del aire sin perturbar a esa altura
!                  (1/3/2000)

      if (ngot(2).ge.1 .or. ncri(2).gt.1) then

         impx=0.
         impy=0.
         zmed=0.
         Qaguat=0.
         spos=1
         Xnub(tte)=Xnub(tte-1)
         Ynub(tte)=Ynub(tte-1)

         laux1=min(lgot(1),lcri(1))
         laux2=max(lgot(2),lcri(2))
         maux1=min(mgot(1),mcri(1))
         maux2=max(mgot(2),mcri(2))
         naux2=max(ngot(2),ncri(2))

         do 1300 k=1,naux2
            do 1300 i=laux1,laux2
               do 1300 j=maux1,maux2

                  Qagua=Qgot1(i,j,k)+Qcri1(i,j,k)+Qllu1(i,j,k)+&
                       Qnie1(i,j,k)+Qgra1(i,j,k)
                  Qagua=Qgot1(i,j,k)+Qllu1(i,j,k)+&
                       (Qcri1(i,j,k)+Qnie1(i,j,k)+Qgra1(i,j,k))/1000
                  zmed=zmed+k*Qagua
                  Qaguat=Qaguat+Qagua

 1300    continue

         if (Qaguat.gt.1e-3) then
            zmed=zmed/Qaguat
            Xnub(tte)=Xnub(tte)+UU(nint(zmed))*lte
            Ynub(tte)=Ynub(tte)+VV(nint(zmed))*lte

!          write(*,*) 'posnub'
!          write(*,*) zmed,Qaguat,Xnub(tte),UU(nint(zmed))
!     &              ,Ynub(tte),VV(nint(zmed)),lte,tte
!               pause

         endif

      endif










