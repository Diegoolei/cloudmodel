*     movimiento de la nube (4/01/99)
*     Redefine el valor de todas las variables (deberian 
*      ser las que se graban solamente)

*     calculo de la posicion media de la nube, es decir de las gotitas
*     el centro esta inicialmente en nx1/2+.5, nx1/2+.5
*     posxx y posyy son siempre en modulo menores que dx1
*     En posx y posy se guarda para cada tt la posicion en 
*     puntos de red

      if (spos .eq. 1) then
        posxx=Xnub(tte)
        posyy=Ynub(tte)
      else
        posxx=0.
        posyy=0.
      endif

      posx(tte)=posx(tte-1)
      posy(tte)=posy(tte-1)
      
**    corrimiento en x

      if (posxx.gt.dx1) then
        posx(tte)=posx(tte)+1
        Xnub(tte)=Xnub(tte)-dx1

*##
        write(*,*) 'corri en x pos'

        do 1500 k=0,nz1+1
          do 1501 j=0,nx1+1
          do 1502 i=1,nx1+1
            U1(i-1,j,k)=U1(i,j,k)
            V1(i-1,j,k)=V1(i,j,k)
            W1(i-1,j,k)=W1(i,j,k)
            Pres1(i-1,j,k)=Pres1(i,j,k)
            U2(i-1,j,k)=U2(i,j,k)
            V2(i-1,j,k)=V2(i,j,k)
            W2(i-1,j,k)=W2(i,j,k)
            Pres2(i-1,j,k)=Pres2(i,j,k)
            Titaa1(i-1,j,k)=Titaa1(i,j,k)
            Qvap1(i-1,j,k)=Qvap1(i,j,k)
            Qgot1(i-1,j,k)=Qgot1(i,j,k)
            Qllu1(i-1,j,k)=Qllu1(i,j,k)
            Qcri1(i-1,j,k)=Qcri1(i,j,k)
            Aer1(i-1,j,k)=Aer1(i,j,k)
            Fcalo(i-1,j,k)=Fcalo(i,j,k)
 1502     continue
          i=nx1+1
            U1(i,j,k)=U1(i-1,j,k)
            V1(i,j,k)=V1(i-1,j,k)
            W1(i,j,k)=W1(i-1,j,k)
            Pres1(i,j,k)=Pres1(i-1,j,k)
            U2(i,j,k)=U2(i-1,j,k)
            V2(i,j,k)=V2(i-1,j,k)
            W2(i,j,k)=W2(i-1,j,k)
            Pres2(i,j,k)=Pres2(i-1,j,k)
            Titaa1(i,j,k)=Titaa1(i-1,j,k)
            Qvap1(i,j,k)=Qvap1(i-1,j,k)
            Qgot1(i,j,k)=Qgot1(i-1,j,k)
            Qllu1(i,j,k)=Qllu1(i-1,j,k)
            Qcri1(i,j,k)=Qcri1(i-1,j,k)
            Aer1(i,j,k)=Aer1(i-1,j,k)
            Fcalo(i,j,k)=0.
 1501     continue
          i=nx1
          j=0
            U1(i,j,k)=(U1(i-1,j,k)+U1(i,j+1,k))/2.
            V1(i,j,k)=(V1(i-1,j,k)+V1(i,j+1,k))/2.
            W1(i,j,k)=(W1(i-1,j,k)+W1(i,j+1,k))/2.
            Pres1(i,j,k)=(Pres1(i-1,j,k)+Pres1(i,j+1,k))/2.
            U2(i,j,k)=(U2(i-1,j,k)+U2(i,j+1,k))/2.
            V2(i,j,k)=(V2(i-1,j,k)+V2(i,j+1,k))/2.
            W2(i,j,k)=(W2(i-1,j,k)+W2(i,j+1,k))/2.
            Pres2(i,j,k)=(Pres2(i-1,j,k)+Pres2(i,j+1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i-1,j,k)+Titaa1(i,j+1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i-1,j,k)+Qvap1(i,j+1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i-1,j,k)+Qgot1(i,j+1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i-1,j,k)+Qllu1(i,j+1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i-1,j,k)+Qcri1(i,j+1,k))/2.
            Aer1(i,j,k)=(Aer1(i-1,j,k)+Aer1(i,j+1,k))/2.
            Fcalo(i,j,k)=0.
          j=nx1+1
            U1(i,j,k)=(U1(i-1,j,k)+U1(i,j-1,k))/2.
            V1(i,j,k)=(V1(i-1,j,k)+V1(i,j-1,k))/2.
            W1(i,j,k)=(W1(i-1,j,k)+W1(i,j-1,k))/2.
            Pres1(i,j,k)=(Pres1(i-1,j,k)+Pres1(i,j-1,k))/2.
            U2(i,j,k)=(U2(i-1,j,k)+U2(i,j-1,k))/2.
            V2(i,j,k)=(V2(i-1,j,k)+V2(i,j-1,k))/2.
            W2(i,j,k)=(W2(i-1,j,k)+W2(i,j-1,k))/2.
            Pres2(i,j,k)=(Pres2(i-1,j,k)+Pres2(i,j-1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i-1,j,k)+Titaa1(i,j-1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i-1,j,k)+Qvap1(i,j-1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i-1,j,k)+Qgot1(i,j-1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i-1,j,k)+Qllu1(i,j-1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i-1,j,k)+Qcri1(i,j-1,k))/2.
            Aer1(i,j,k)=(Aer1(i-1,j,k)+Aer1(i,j-1,k))/2.
            Fcalo(i,j,k)=0.
 1500   continue
      endif
 
      if (posxx.lt.-dx1) then
        posx(tte)=posx(tte)-1
        Xnub(tte)=Xnub(tte)+dx1

*##
        write(*,*) 'corri en x neg'

        do 1510 k=0,nz1+1
          do 1511 j=0,nx1+1
          do 1512 i=nx1,0,-1
            U1(i+1,j,k)=U1(i,j,k)
            V1(i+1,j,k)=V1(i,j,k)
            W1(i+1,j,k)=W1(i,j,k)
            Pres1(i+1,j,k)=Pres1(i,j,k)
            U2(i+1,j,k)=U2(i,j,k)
            V2(i+1,j,k)=V2(i,j,k)
            W2(i+1,j,k)=W2(i,j,k)
            Pres2(i+1,j,k)=Pres2(i,j,k)
            Titaa1(i+1,j,k)=Titaa1(i,j,k)
            Qvap1(i+1,j,k)=Qvap1(i,j,k)
            Qgot1(i+1,j,k)=Qgot1(i,j,k)
            Qllu1(i+1,j,k)=Qllu1(i,j,k)
            Qcri1(i+1,j,k)=Qcri1(i,j,k)
            Aer1(i+1,j,k)=Aer1(i,j,k)
            Fcalo(i+1,j,k)=Fcalo(i,j,k)
 1512     continue
          i=0
            U1(i,j,k)=U1(i+1,j,k)
            V1(i,j,k)=V1(i+1,j,k)
            W1(i,j,k)=W1(i+1,j,k)
            Pres1(i,j,k)=Pres1(i+1,j,k)
            U2(i,j,k)=U2(i+1,j,k)
            V2(i,j,k)=V2(i+1,j,k)
            W2(i,j,k)=W2(i+1,j,k)
            Pres2(i,j,k)=Pres2(i+1,j,k)
            Titaa1(i,j,k)=Titaa1(i+1,j,k)
            Qvap1(i,j,k)=Qvap1(i+1,j,k)
            Qgot1(i,j,k)=Qgot1(i+1,j,k)
            Qllu1(i,j,k)=Qllu1(i+1,j,k)
            Qcri1(i,j,k)=Qcri1(i+1,j,k)
            Aer1(i,j,k)=Aer1(i+1,j,k)
            Fcalo(i,j,k)=0.
 1511     continue
          i=1
          j=0
            U1(i,j,k)=(U1(i+1,j,k)+U1(i,j+1,k))/2.
            V1(i,j,k)=(V1(i+1,j,k)+V1(i,j+1,k))/2.
            W1(i,j,k)=(W1(i+1,j,k)+W1(i,j+1,k))/2.
            Pres1(i,j,k)=(Pres1(i+1,j,k)+Pres1(i,j+1,k))/2.
            U2(i,j,k)=(U2(i+1,j,k)+U2(i,j+1,k))/2.
            V2(i,j,k)=(V2(i+1,j,k)+V2(i,j+1,k))/2.
            W2(i,j,k)=(W2(i+1,j,k)+W2(i,j+1,k))/2.
            Pres2(i,j,k)=(Pres2(i+1,j,k)+Pres2(i,j+1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i+1,j,k)+Titaa1(i,j+1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i+1,j,k)+Qvap1(i,j+1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i+1,j,k)+Qgot1(i,j+1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i+1,j,k)+Qllu1(i,j+1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i+1,j,k)+Qcri1(i,j+1,k))/2.
            Aer1(i,j,k)=(Aer1(i+1,j,k)+Aer1(i,j+1,k))/2.
            Fcalo(i,j,k)=0.
          j=nx1+1
            U1(i,j,k)=(U1(i+1,j,k)+U1(i,j-1,k))/2.
            V1(i,j,k)=(V1(i+1,j,k)+V1(i,j-1,k))/2.
            W1(i,j,k)=(W1(i+1,j,k)+W1(i,j-1,k))/2.
            Pres1(i,j,k)=(Pres1(i+1,j,k)+Pres1(i,j-1,k))/2.
            U2(i,j,k)=(U2(i+1,j,k)+U2(i,j-1,k))/2.
            V2(i,j,k)=(V2(i+1,j,k)+V2(i,j-1,k))/2.
            W2(i,j,k)=(W2(i+1,j,k)+W2(i,j-1,k))/2.
            Pres2(i,j,k)=(Pres2(i+1,j,k)+Pres2(i,j-1,k))/2.
            Titaa1(i,j,k)=(Titaa1(i+1,j,k)+Titaa1(i,j-1,k))/2.
            Qvap1(i,j,k)=(Qvap1(i+1,j,k)+Qvap1(i,j-1,k))/2.
            Qgot1(i,j,k)=(Qgot1(i+1,j,k)+Qgot1(i,j-1,k))/2.
            Qllu1(i,j,k)=(Qllu1(i+1,j,k)+Qllu1(i,j-1,k))/2.
            Qcri1(i,j,k)=(Qcri1(i+1,j,k)+Qcri1(i,j-1,k))/2.
            Aer1(i,j,k)=(Aer1(i+1,j,k)+Aer1(i,j-1,k))/2.
            Fcalo(i,j,k)=0.
 1510   continue
      endif

**    corrimiento en y

      if (posyy.gt.dx1) then
        posy(tte)=posy(tte)+1
        Ynub(tte)=Ynub(tte)-dx1

*##
        write(*,*) 'corri en y pos'

        do 1520 k=0,nz1+1
          do 1521 i=0,nx1+1
          do 1522 j=1,nx1+1
            U1(i,j-1,k)=U1(i,j,k)
            V1(i,j-1,k)=V1(i,j,k)
            W1(i,j-1,k)=W1(i,j,k)
            Pres1(i,j-1,k)=Pres1(i,j,k)
            U2(i,j-1,k)=U2(i,j,k)
            V2(i,j-1,k)=V2(i,j,k)
            W2(i,j-1,k)=W2(i,j,k)
            Pres2(i,j-1,k)=Pres2(i,j,k)
            Titaa1(i,j-1,k)=Titaa1(i,j,k)
            Qvap1(i,j-1,k)=Qvap1(i,j,k)
            Qgot1(i,j-1,k)=Qgot1(i,j,k)
            Qllu1(i,j-1,k)=Qllu1(i,j,k)
            Qcri1(i,j-1,k)=Qcri1(i,j,k)
            Aer1(i,j-1,k)=Aer1(i,j,k)
            Fcalo(i,j-1,k)=Fcalo(i,j,k)
 1522     continue
          j=nx1+1
            U1(i,j,k)=U1(i,j-1,k)
            V1(i,j,k)=V1(i,j-1,k)
            W1(i,j,k)=W1(i,j-1,k)
            Pres1(i,j,k)=Pres1(i,j-1,k)
            U2(i,j,k)=U2(i,j-1,k)
            V2(i,j,k)=V2(i,j-1,k)
            W2(i,j,k)=W2(i,j-1,k)
            Pres2(i,j,k)=Pres2(i,j-1,k)
            Titaa1(i,j,k)=Titaa1(i,j-1,k)
            Qvap1(i,j,k)=Qvap1(i,j-1,k)
            Qgot1(i,j,k)=Qgot1(i,j-1,k)
            Qllu1(i,j,k)=Qllu1(i,j-1,k)
            Qcri1(i,j,k)=Qcri1(i,j-1,k)
            Aer1(i,j,k)=Aer1(i,j-1,k)
            Fcalo(i,j,k)=0.
 1521     continue
          j=nx1
          i=0
            U1(i,j,k)=(U1(i,j-1,k)+U1(i+1,j,k))/2.
            V1(i,j,k)=(V1(i,j-1,k)+V1(i+1,j,k))/2.
            W1(i,j,k)=(W1(i,j-1,k)+W1(i+1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j-1,k)+Pres1(i+1,j,k))/2.
            U2(i,j,k)=(U2(i,j-1,k)+U2(i+1,j,k))/2.
            V2(i,j,k)=(V2(i,j-1,k)+V2(i+1,j,k))/2.
            W2(i,j,k)=(W2(i,j-1,k)+W2(i+1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j-1,k)+Pres2(i+1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j-1,k)+Titaa1(i+1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j-1,k)+Qvap1(i+1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j-1,k)+Qgot1(i+1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j-1,k)+Qllu1(i+1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j-1,k)+Qcri1(i+1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j-1,k)+Aer1(i+1,j,k))/2.
            Fcalo(i,j,k)=0.
          i=nx1+1
            U1(i,j,k)=(U1(i,j-1,k)+U1(i-1,j,k))/2.
            V1(i,j,k)=(V1(i,j-1,k)+V1(i-1,j,k))/2.
            W1(i,j,k)=(W1(i,j-1,k)+W1(i-1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j-1,k)+Pres1(i-1,j,k))/2.
            U2(i,j,k)=(U2(i,j-1,k)+U2(i-1,j,k))/2.
            V2(i,j,k)=(V2(i,j-1,k)+V2(i-1,j,k))/2.
            W2(i,j,k)=(W2(i,j-1,k)+W2(i-1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j-1,k)+Pres2(i-1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j-1,k)+Titaa1(i-1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j-1,k)+Qvap1(i-1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j-1,k)+Qgot1(i-1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j-1,k)+Qllu1(i-1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j-1,k)+Qcri1(i-1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j-1,k)+Aer1(i-1,j,k))/2.
            Fcalo(i,j,k)=0.
 1520   continue
      endif
 
      if (posyy.lt.-dx1) then
        posy(tte)=posy(tte)-1
        Xnub(tte)=Xnub(tte)+dx1

*##
        write(*,*) 'corri en y neg'

        do 1530 k=0,nz1+1
          do 1531 i=0,nx1+1
          do 1532 j=nx1,0,-1
            U1(i,j+1,k)=U1(i,j,k)
            V1(i,j+1,k)=V1(i,j,k)
            W1(i,j+1,k)=W1(i,j,k)
            Pres1(i,j+1,k)=Pres1(i,j,k)
            U2(i,j+1,k)=U2(i,j,k)
            V2(i,j+1,k)=V2(i,j,k)
            W2(i,j+1,k)=W2(i,j,k)
            Pres2(i,j+1,k)=Pres2(i,j,k)
            Titaa1(i,j+1,k)=Titaa1(i,j,k)
            Qvap1(i,j+1,k)=Qvap1(i,j,k)
            Qgot1(i,j+1,k)=Qgot1(i,j,k)
            Qllu1(i,j+1,k)=Qllu1(i,j,k)
            Qcri1(i,j+1,k)=Qcri1(i,j,k)
            Aer1(i,j+1,k)=Aer1(i,j,k)
            Fcalo(i,j+1,k)=Fcalo(i,j,k)
 1532     continue
          j=0
            U1(i,j,k)=U1(i,j-1,k)
            V1(i,j,k)=V1(i,j-1,k)
            W1(i,j,k)=W1(i,j-1,k)
            Pres1(i,j,k)=Pres1(i,j-1,k)
            U2(i,j,k)=U2(i,j-1,k)
            V2(i,j,k)=V2(i,j-1,k)
            W2(i,j,k)=W2(i,j-1,k)
            Pres2(i,j,k)=Pres2(i,j-1,k)
            Titaa1(i,j,k)=Titaa1(i,j-1,k)
            Qvap1(i,j,k)=Qvap1(i,j-1,k)
            Qgot1(i,j,k)=Qgot1(i,j-1,k)
            Qllu1(i,j,k)=Qllu1(i,j-1,k)
            Qcri1(i,j,k)=Qcri1(i,j-1,k)
            Aer1(i,j,k)=Aer1(i,j-1,k)
            Fcalo(i,j,k)=0.
 1531     continue
          j=1
          i=0
            U1(i,j,k)=(U1(i,j+1,k)+U1(i+1,j,k))/2.
            V1(i,j,k)=(V1(i,j+1,k)+V1(i+1,j,k))/2.
            W1(i,j,k)=(W1(i,j+1,k)+W1(i+1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j+1,k)+Pres1(i+1,j,k))/2.
            U2(i,j,k)=(U2(i,j+1,k)+U2(i+1,j,k))/2.
            V2(i,j,k)=(V2(i,j+1,k)+V2(i+1,j,k))/2.
            W2(i,j,k)=(W2(i,j+1,k)+W2(i+1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j+1,k)+Pres2(i+1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j+1,k)+Titaa1(i+1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j+1,k)+Qvap1(i+1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j+1,k)+Qgot1(i+1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j+1,k)+Qllu1(i+1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j+1,k)+Qcri1(i+1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j+1,k)+Aer1(i+1,j,k))/2.
            Fcalo(i,j,k)=0.
          i=nx1+1
            U1(i,j,k)=(U1(i,j+1,k)+U1(i-1,j,k))/2.
            V1(i,j,k)=(V1(i,j+1,k)+V1(i-1,j,k))/2.
            W1(i,j,k)=(W1(i,j+1,k)+W1(i-1,j,k))/2.
            Pres1(i,j,k)=(Pres1(i,j+1,k)+Pres1(i-1,j,k))/2.
            U2(i,j,k)=(U2(i,j+1,k)+U2(i-1,j,k))/2.
            V2(i,j,k)=(V2(i,j+1,k)+V2(i-1,j,k))/2.
            W2(i,j,k)=(W2(i,j+1,k)+W2(i-1,j,k))/2.
            Pres2(i,j,k)=(Pres2(i,j+1,k)+Pres2(i-1,j,k))/2.
            Titaa1(i,j,k)=(Titaa1(i,j+1,k)+Titaa1(i-1,j,k))/2.
            Qvap1(i,j,k)=(Qvap1(i,j+1,k)+Qvap1(i-1,j,k))/2.
            Qgot1(i,j,k)=(Qgot1(i,j+1,k)+Qgot1(i-1,j,k))/2.
            Qllu1(i,j,k)=(Qllu1(i,j+1,k)+Qllu1(i-1,j,k))/2.
            Qcri1(i,j,k)=(Qcri1(i,j+1,k)+Qcri1(i-1,j,k))/2.
            Aer1(i,j,k)=(Aer1(i,j+1,k)+Aer1(i-1,j,k))/2.
            Fcalo(i,j,k)=0.
 1530   continue

      endif

      posxx=posx(tte)*dx1+Xnub(tte)
      posyy=posy(tte)*dx1+Ynub(tte)
