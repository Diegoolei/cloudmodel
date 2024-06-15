module model_aux
contains
   subroutine vapor_advection()
      !######################## Adveccion de vapores #######################
      USE advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advaer1
      USE permic, only: aer1, Qvap1, Qllu1, Qnie1, Qgra1
      USE perdim, only: W2
      USE dimen, only: nx1
      implicit none
      integer :: i, j
      do i=0,nx1+1
         do j=0,nx1+1
            advvap1(i,j)=W2(i,j,1)*(Qvap1(i,j,1)+Qvap1(i,j,0))/4.
            advgot1(i,j)=0.
            advllu1(i,j)=W2(i,j,1)*Qllu1(i,j,1)
            if (W2(i,j,1).gt.0) advllu1(i,j)=0.
            advaer1(i,j)=W2(i,j,1)*(aer1(i,j,1)+aer1(i,j,0))/4.
            if(W2(i,j,1).lt.0) advaer1(i,j)=advaer1(i,j)*1.5
            advcri1(i,j)=0.
            advnie1(i,j)=W2(i,j,1)*Qnie1(i,j,1)
            if (W2(i,j,1).gt.0) advnie1(i,j)=0.
            advgra1(i,j)=W2(i,j,1)*Qgra1(i,j,1)
            if (W2(i,j,1).gt.0) advgra1(i,j)=0.
         end do
      end do
   end subroutine vapor_advection

   subroutine dinamics()
      !########### calculo de la dinamica y de la termodinamica ############
      USE dimen, only: dt1, dx1, nx1, nz1
      USE model_var, only: dden0z, ener1, s, llluneg, lcrineg, lnieneg, lgraneg,&
         lvapneg, laerneg, Qvapneg, aerneg
      USE estbas, only: Den0, Qvap0, aer0
      USE perdim, only: W2, U2, V2, Fcalo
      USE advecs, only: advllu1, advaer1, advnie1, advgra1, advvap1, advgot1,&
         advcri1, advvap2, advgot2, advllu2, advcri2, advnie2, advgra2, advaer1,&
         advaer2
      USE permic, only: Qgot2, Qllu2, Qcri2, Qnie2, Qgra2, Qvap2, aer2
      USE lmngot, only: lgot, mgot, ngot
      USE lmnllu, only: lllu, mllu, nllu
      USE lmncri, only: lcri, mcri, ncri
      USE lmnnie, only: lnie, mnie, nnie
      USE lmngra, only: lgra, mgra, ngra
      implicit none
      integer :: i, j, k, l, m, n

      s=0
      Qvapneg=0.
      lvapneg=0
      aerneg=0.
      laerneg=0
      llluneg=0
      lcrineg=0.
      lnieneg=0.
      lgraneg=0.
      ener1=0.
      lgot(1)=nx1
      lgot(2)=0
      mgot(1)=nx1
      mgot(2)=0
      ngot(1)=nz1
      ngot(2)=0
      lllu(1)=nx1
      lllu(2)=0
      mllu(1)=nx1
      mllu(2)=0
      nllu(1)=nz1
      nllu(2)=0
      lcri(1)=nx1
      lcri(2)=0
      mcri(1)=nx1
      mcri(2)=0
      ncri(1)=nz1
      ncri(2)=0
      lnie(1)=nx1
      lnie(2)=0
      mnie(1)=nx1
      mnie(2)=0
      nnie(1)=nz1
      nnie(2)=0
      lgra(1)=nx1
      lgra(2)=0
      mgra(1)=nx1
      mgra(2)=0
      ngra(1)=nz1
      ngra(2)=0
      do k=1,nz1-1
         n=k
         dden0z=(Den0(k+1)-Den0(k-1))/Den0(k)
         call turbu1(n)
         do i=1,nx1
            l=i
            do j=1,nx1
               m=j
               !calculo del coeficiente de turbulencia y derivadas
               call turbu2(l,m,n)

               !calculo de las inhomogeneidades para las velocidades
               call inomo(l,m,n,dden0z)

               !calculo de la energia cinetica
               ener1=.5*Den0(k)*(U2(i,j,k)**2.+V2(i,j,k)**2.+W2(i,j,k)**2.)+ener1

               !calculo de la temperatura potencial
               call tempot(l,m,n,dden0z,Fcalo(i,j,k))
               Fcalo(i,j,k)=0.

               !dinamica del vapor y de las gotitas
               call dvapor(l,m,n)
               advvap1(i,j)=advvap2(i,j)
               call dgotit(l,m,n)
               advgot1(i,j)=advgot2(i,j)
               call dlluvi(l,m,n)
               advllu1(i,j)=advllu2(i,j)
               call dcrist(l,m,n)
               advcri1(i,j)=advcri2(i,j)
               call dnieve(l,m,n)
               advnie1(i,j)=advnie2(i,j)
               call dgrani(l,m,n)
               advgra1(i,j)=advgra2(i,j)
               call daeros(l,m,n)
               advaer1(i,j)=advaer2(i,j)

               !limites de la nube
               if(Qgot2(i,j,k).ne.0) then
                  if (i.lt.lgot(1)) lgot(1)=i
                  if (i.gt.lgot(2)) lgot(2)=i
                  if (j.lt.mgot(1)) mgot(1)=j
                  if (j.gt.mgot(2)) mgot(2)=j
                  if (k.lt.ngot(1)) ngot(1)=k
                  if (k.gt.ngot(2)) ngot(2)=k
                  s=s+1
               endif

               !limites de la lluvia
               if(Qllu2(i,j,k).ne.0) then
                  if (i.lt.lllu(1)) lllu(1)=i
                  if (i.gt.lllu(2)) lllu(2)=i
                  if (j.lt.mllu(1)) mllu(1)=j
                  if (j.gt.mllu(2)) mllu(2)=j
                  if (k.lt.nllu(1)) nllu(1)=k
                  if (k.gt.nllu(2)) nllu(2)=k
                  llluneg=1
               endif

               !limites de los cristales
               if(Qcri2(i,j,k).ne.0) then
                  if (i.lt.lcri(1)) lcri(1)=i
                  if (i.gt.lcri(2)) lcri(2)=i
                  if (j.lt.mcri(1)) mcri(1)=j
                  if (j.gt.mcri(2)) mcri(2)=j
                  if (k.lt.ncri(1)) ncri(1)=k
                  if (k.gt.ncri(2)) ncri(2)=k
                  lcrineg=1
               endif

               !limites de la nieve
               if(Qnie2(i,j,k).ne.0) then
                  if (i.lt.lnie(1)) lnie(1)=i
                  if (i.gt.lnie(2)) lnie(2)=i
                  if (j.lt.mnie(1)) mnie(1)=j
                  if (j.gt.mnie(2)) mnie(2)=j
                  if (k.lt.nnie(1)) nnie(1)=k
                  if (k.gt.nnie(2)) nnie(2)=k
                  lnieneg=1
               endif

               !limites del granizo
               if(Qgra2(i,j,k).ne.0) then
                  if (i.lt.lgra(1)) lgra(1)=i
                  if (i.gt.lgra(2)) lgra(2)=i
                  if (j.lt.mgra(1)) mgra(1)=j
                  if (j.gt.mgra(2)) mgra(2)=j
                  if (k.lt.ngra(1)) ngra(1)=k
                  if (k.gt.ngra(2)) ngra(2)=k
                  lgraneg=1
               endif

               if(Qvap0(k)+Qvap2(i,j,k).lt.0) then
                  Qvapneg=Qvapneg+Qvap0(k)+Qvap2(i,j,k)
                  lvapneg=1
               endif

               if(aer0(k)+aer2(i,j,k).lt.0) then
                  aerneg=aerneg+aer0(k)+aer2(i,j,k)
                  laerneg=1
               endif
            end do
         end do
      end do

   end subroutine dinamics

   subroutine negative_correction
      USE model_var, only: s, llluneg, lcrineg, lnieneg, lgraneg, lvapneg, laerneg,&
         Qvapneg, aerneg
      implicit none
      if(s.ge.1) call corgot
      if (llluneg.eq.1) call corllu
      if (lcrineg.eq.1) call corcri
      if (lnieneg.eq.1) call cornie
      if (lgraneg.eq.1) call corgra
      if (lvapneg.eq.1) call corvap(Qvapneg)
      if (laerneg.eq.1) call coraer(aerneg)
   end subroutine negative_correction
end module model_aux
