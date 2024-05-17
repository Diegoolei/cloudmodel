!      include para variables del modelo

      integer tt,t1,t2,n,m,l,i,j,k,lll,s,iT,tte
!$$
      integer lvapneg,llluneg,lcrineg,laerneg,lnieneg,lgraneg,yy

      real T,P
      real Dv,Lvl,Lvs,Lsl,Vis,Qvap,Qliq,densi,nu
      real Lsl00
      real Eaucn,Eaccn,Eacng
!$$
      real*8 qgotaux,qvapaux,qlluaux,qcriaux,qnieaux,qgraaux,aeraux
      real*8 auxx,auxy,auxz
      real*8 Taux,Qvapneg,aerneg
      real Naer,dqgot,dqcri,daer,daer2
      real Fcal
      real elvs,esvs,e1,rl,rs,dden0z
      real aux,aux1,aux2,aux3,aux4
      real*8 ener,ener1,ener2,ener3,ener4,ener5,qv,qg,daitot

      real*8 vapt1,vapt2,vapt3,vapt4
      real*8 gott1,gott2,gott3,gott4
      real*8 aert1,aert2,aert3,aert4
      real*8 totnuc,totmic

      real*8 impx,impy,Qagua,Qaguat
      real*8 Xnub(5000),Ynub(5000)
      real   posxx,posyy,zmed
      integer posx(-3:5000),posy(-3:5000),spos

      character*50 bre
      character*12 nombre
      character*2 tie

      integer laux1,laux2,maux1,maux2,naux2
      integer umax,umin,vmax,vmin,wmax,wmin,titamax,titamin
!$$
      integer qvapmax,qvapmin,qgotmax,qllumax,qcrimax,qniemax,qgramax
      integer aermax
      integer lumax,mumax,numax,lumin,mumin,numin
      integer lvmax,mvmax,nvmax,lvmin,mvmin,nvmin
      integer lwmax,mwmax,nwmax,lwmin,mwmin,nwmin
      integer ltitamax,mtitamax,ntitamax,ltitamin,mtitamin
      integer ntitamin,lqvapmax,mqvapmax,nqvapmax,lqvapmin
      integer mqvapmin,nqvapmin,lqgotmax,mqgotmax,nqgotmax
      integer lqllumax,mqllumax,nqllumax,laermax,maermax,naermax
!$$
      integer lqcrimax,mqcrimax,nqcrimax
      integer lqniemax,mqniemax,nqniemax
      integer lqgramax,mqgramax,nqgramax
      integer qgottot,qllutot,qcritot,qnietot,qgratot

      real cks,turbu,lapla      

      bre='01020304050607080910111213141516171819202122232425'
      tie='31'      
      
      ini=0                            !inicio por vez primera= 0
      t1=0                             !paso a inicio (si ini=0->t1=0)
      ltt = 45. * 60.
!      ltt=25.*3.*60.                      !tiempo total de simulacion
!      ltt=15.*6.*60.                       !tiempo total de simulacion
!      ltt=3.*60.                          !tiempo total de simulacion
!      ltt=2.*2.                           !tiempo total de simulacion
!      ltg=3*60.                           !tiempo de grabacion
!      ltg=6.*60.                          !tiempo de grabacion
      ltg= 3. * 60.                              !tiempo de grabacion
      lte= 45. * 60.                             !tiempo de grabacion estadistica
!      lte=2.                              !tiempo de grabacion estadistica
!      ltb=3*60.                           ! tiempo de grabacion de backup
      ltb= 45. * 60.     


      ctur=0.5
      
      pro1=1.-2e-2*(dt1/5.)
      pro2=(1.-pro1)/6.
      pro3=1.-2e-2*(dt1/5.)
      pro4=(1.-pro1)/4.
            
      lt1=nint(ltt/dt1)
      lt2=nint(dt1/dt2)
      lt3=2*nint(dt1/dt3)
      cteturb=ctur/2.**.5

      cks=cteturb*2.

      dx2=2.*dx1
      dx8=8.*dx1
      dx12=12.*dx1
      AA=1./Eps-1.
      ikapa=1./Kapa

      cteqgot=160./3**6.*pi*rhow*N0got
      cteqllu=8.*pi*rhow*N0llu
      cteqnie=.6*pi*rhonie*N0nie
      cteqgra=8.*pi*rhogra*N0gra


      spos=0
      posxx=0.
      posyy=0.
      posx(0)=0
      posy(0)=0
      tte=0











