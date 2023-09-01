*     calculo de algunas cantidades de interes (12/02/99)
        
        umax=0.
          lumax=0
          mumax=0
          numax=0
        vmax=0.
          lvmax=0
          mvmax=0
          nvmax=0
        wmax=0.
          lwmax=0
          mwmax=0
          nwmax=0
        titamax=0.
          ltitamax=0
          mtitamax=0
          ntitamax=0
        qvapmax=0.
          lqvapmax=0
          mqvapmax=0
          nqvapmax=0
        qgotmax=0.
          lqgotmax=0
          mqgotmax=0
          nqgotmax=0
        qllumax=0.
          lqllumax=0
          mqllumax=0
          nqllumax=0
        qcrimax=0.
          lqcrimax=0
          mqcrimax=0
          nqcrimax=0
*$$
        qniemax=0.
          lqniemax=0
          mqniemax=0
          nqniemax=0
        qgramax=0.
          lqgramax=0
          mqgramax=0
          nqgramax=0
        aermax=0.
          laermax=0
          maermax=0
          naermax=0

        umin=0.
          lumin=0
          mumin=0
          numin=0
        vmin=0.
          lvmin=0
          mvmin=0
          nvmin=0
        wmin=0.
          lwmin=0
          mwmin=0
          nwmin=0
        titamin=0.
          ltitamin=0
          mtitamin=0
          ntitamin=0
        qvapmin=0.
          lqvapmin=0
          mqvapmin=0
          nqvapmin=0

        qgottot=0.
        qllutot=0.
        qcritot=0.
        qnietot=0.
        qgratot=0.

        do 700 k=1,nz1
        do 700 i=1,nx1
        do 700 j=1,nx1
          
          if (umax.lt.U1(i,j,k)*100) then
            umax=U1(i,j,k)*100
            lumax=i
            mumax=j
            numax=k
          endif

          if (umin.gt.U1(i,j,k)*100) then
            umin=U1(i,j,k)*100
            lumin=i
            mumin=j
            numin=k
          endif

          if (vmax.lt.V1(i,j,k)*100) then
            vmax=V1(i,j,k)*100
            lvmax=i
            mvmax=j
            nvmax=k
          endif

          if (vmin.gt.V1(i,j,k)*100) then
            vmin=V1(i,j,k)*100
            lvmin=i
            mvmin=j
            nvmin=k
          endif

          if (wmax.lt.W1(i,j,k)*100) then
            wmax=W1(i,j,k)*100
            lwmax=i
            mwmax=j
            nwmax=k
          endif

          if (wmin.gt.W1(i,j,k)*100) then
            wmin=W1(i,j,k)*100
            lwmin=i
            mwmin=j
            nwmin=k
          endif

          if (titamax.lt.Titaa1(i,j,k)*1000) then
            titamax=Titaa1(i,j,k)*1000
            ltitamax=i
            mtitamax=j
            ntitamax=k
          endif
          
          if (titamin.gt.Titaa1(i,j,k)*1000) then
            titamin=Titaa1(i,j,k)*1000
            ltitamin=i
            mtitamin=j
            ntitamin=k
          endif
          
          if (qvapmax.lt.Qvap1(i,j,k)*1e6) then
            qvapmax=Qvap1(i,j,k)*1e6
            lqvapmax=i
            mqvapmax=j
            nqvapmax=k
          endif
          
          if (qvapmin.gt.Qvap1(i,j,k)*1e6) then
            qvapmin=Qvap1(i,j,k)*1e6
            lqvapmin=i
            mqvapmin=j
            nqvapmin=k
          endif
          
          if (qgotmax.lt.Qgot1(i,j,k)*1e6) then
            qgotmax=Qgot1(i,j,k)*1e6
            lqgotmax=i
            mqgotmax=j
            nqgotmax=k
          endif
            qgottot=qgottot+Qgot1(i,j,k)*1e6
          
          if (qllumax.lt.Qllu1(i,j,k)*1e6) then
            qllumax=Qllu1(i,j,k)*1e6
            lqllumax=i
            mqllumax=j
            nqllumax=k
          endif
            qllutot=qllutot+Qllu1(i,j,k)*1e6

          if (qcrimax.lt.Qcri1(i,j,k)*1e6) then
            qcrimax=Qcri1(i,j,k)*1e6
            lqcrimax=i
            mqcrimax=j
            nqcrimax=k


*            write(*,*) qcrimax,Qcri1(i,j,k),i,j,k
*            pause

          endif
            qcritot=qcritot+Qcri1(i,j,k)*1e6

*$$          
          if (qniemax.lt.Qnie1(i,j,k)*1e6) then
            qniemax=Qnie1(i,j,k)*1e6
            lqniemax=i
            mqniemax=j
            nqniemax=k
          endif
            qnietot=qnietot+Qnie1(i,j,k)*1e6

          if (qgramax.lt.Qgra1(i,j,k)*1e6) then
            qgramax=Qgra1(i,j,k)*1e6
            lqgramax=i
            mqgramax=j
            nqgramax=k
          endif
            qgratot=qgratot+Qgra1(i,j,k)*1e6

          if (aermax.lt.aer1(i,j,k)/1000) then
            aermax=aer1(i,j,k)/1000
            laermax=i
            maermax=j
            naermax=k
          endif
  700   continue


        qgotmax=0.
        qllumax=0.
        qcrimax=0.
        qniemax=0.
        qgramax=0.
        do 719 i=-1,1
        do 719 j=-1,1
        do 719 k=-1,1
          qgotmax=qgotmax+1e5*Qgot1(lqgotmax+i,mqgotmax+j,nqgotmax+k)
          qllumax=qllumax+1e5*Qllu1(lqllumax+i,mqllumax+j,nqllumax+k)
          qcrimax=qcrimax+1e5*Qcri1(lqcrimax+i,mqcrimax+j,nqcrimax+k)
          qniemax=qniemax+1e5*Qnie1(lqniemax+i,mqniemax+j,nqniemax+k)
          qgramax=qgramax+1e5*Qgra1(lqgramax+i,mqgramax+j,nqgramax+k)


*        write(*,*) qcrimax,Qcri1(lqcrimax+i,mqcrimax+j,nqcrimax+k),
*     &             lqcrimax+i,mqcrimax+j,nqcrimax+k


  719  continue

       qgotmax=qgotmax/27.
       qllumax=qllumax/27.
       qcrimax=qcrimax/27.
       qniemax=qniemax/27.
       qgramax=qgramax/27.
       umax=umax/10
       umin=umin/10
       vmax=vmax/10
       vmin=vmin/10
       wmax=wmax/10
       wmin=wmin/10
       titamax=titamax/10
       titamin=titamin/10
       qvapmax=qvapmax/10
       qvapmin=qvapmin/10
       qgottot=qgottot/1000
       qllutot=qllutot/1000
       qcritot=qcritot/1000
       qnietot=qnietot/1000
       qgratot=qgratot/1000

*$$       
      write(30,710) umax,umin,vmax,vmin,wmax,wmin,titamax,titamin
     &            ,qvapmax,qvapmin,qgotmax,qllumax,qcrimax,qniemax
     &            ,qgramax,aermax
     &            ,lumax,mumax,numax,lumin,mumin,numin
     &            ,lvmax,mvmax,nvmax,lvmin,mvmin,nvmin
     &            ,lwmax,mwmax,nwmax,lwmin,mwmin,nwmin
     &            ,ltitamax,mtitamax,ntitamax,ltitamin,mtitamin
     &            ,ntitamin,lqvapmax,mqvapmax,nqvapmax,lqvapmin
     &            ,mqvapmin,nqvapmin,lqgotmax,mqgotmax,nqgotmax
     &            ,lqllumax,mqllumax,nqllumax
     &            ,lqcrimax,mqcrimax,nqcrimax
     &            ,lqniemax,mqniemax,nqniemax
     &            ,lqgramax,mqgramax,nqgramax
     &            ,laermax,maermax,naermax


      write(33,715) qgottot,qllutot,qcritot,qnietot,qgratot

 710  format(16i5,48i4)
 715  format(5i9)

