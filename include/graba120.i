      open(unit=40,file='f77_data/inis.da')
      open(unit=41,file='f77_data/velos.da',status='unknown',
     & form='unformatted')
      rewind 41
      open(unit=42,file='f77_data/varconz.da',status='unknown',  
     &   form='unformatted')
      rewind 42

        write(40,*) Den0,Temp0,Tita0,Pres00,Qvap0
     &	               ,cc2,aer0,UU,VV
!$$
        write(41) U1,U2,V1,V2,W1,W2,Titaa1,Titaa2,Pres1,Pres2,
     &            Qvap1,Qvap2,Qgot1,Qgot2,Qllu1,Qllu2,
     &            Qcri1,Qcri2,Qnie1,Qnie2,Qgra1,Qgra2,
     &            aer1,aer2,Fcalo
        write(42)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,
     &            Qvaprel,aerrel,Eautcn,Eacrcn
	

      close(42)
      close(41)
      close(40)
