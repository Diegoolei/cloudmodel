      open(unit=60,file='f77_data/nube'//tie//bre(2*t1-1:2*t1)//'.sal'
     &      ,status='unknown',form='unformatted')

      write(60) U1,V1,W1,Titaa1,Pres1,Qvap1,Qgot1,Qllu1,Qcri1,Qnie1
     &        ,Qgra1,aer1
      close(60)