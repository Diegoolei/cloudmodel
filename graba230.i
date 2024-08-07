      close(3)
      nombre='Tita'//tie//bre(2*t1-1:2*t1)//'.m'
      open(3,file=nombre)
      do 370 i=1,nx1
       write(3,2000) (Titaa1(i,j,0),j=1,nx1)
  370 continue   
      close(3)
      nombre='Qvap'//tie//bre(2*t1-1:2*t1)//'.m'
      open(3,file=nombre)
      do 380 i=1,nx1
       write(3,2000) (Qvap1(i,j,0)+Qvap0(k),j=1,nx1)
  380 continue   
      close(3)
      nombre='Qllu'//tie//bre(2*t1-1:2*t1)//'.m'
      open(3,file=nombre)
      do 395 i=1,nx1
       write(3,2000) (Qllu1(i,j,0),j=1,nx1)
  395 continue   
      close(3)
      nombre='Aero'//tie//bre(2*t1-1:2*t1)//'.m'
      open(3,file=nombre)
      do 385 i=1,nx1
       write(3,2000) (aer1(i,j,0)+aer0(k),j=1,nx1)
  385 continue   
      close(3)

      nombre='Qgra'//tie//bre(2*t1-1:2*t1)//'.m'
      open(3,file=nombre)
      do 525 i=1,nx1
       write(3,2000) (Qgra1(i,j,0),j=1,nx1)
  525 continue   
      close(3)


 2000 format(50E11.3)
