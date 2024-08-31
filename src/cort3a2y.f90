!     Este programa lee la salida del la nube en 3-D y realiza cortes
!     planos para que puedan ser levantados por el Matlab
!     Los cortes son en y en la posicion del maximo de precipitacion y
!     gotitas
module get_cut
   implicit none

contains
   subroutine get_cutted()
      use dimensions
      use constants
      use initial_z_state
      use microphysics_perturbation
      use dinamic_var_perturbation
      use memory_managment
      use model_initialization, only: initialize_model
      implicit none
      character*13 :: directory
      integer ii,i,j,k,l,m,n,ymax
      integer pasoaux(25),lvar
      integer, parameter :: nvar=13
      real*8 aux,maypre
      character*100 paso
      character*2 nombre
      character*25 archaux
      character*2 var(13)
      integer*4 variable(nx1, mod_nz1, nvar)
      integer unit_number
      integer unit_number2

      !**   datos generales
      var(1)='UU'
      var(2)='VV'
      var(3)='WW'
      var(4)='TT'
      var(5)='GO'
      var(6)='LL'
      var(7)='CR'
      var(8)='NI'
      var(9)='GR'
      var(10)='VA'
      var(11)='AE'
      var(12)='V1'
      var(13)='A1'


      paso='01020304050607080910111213141516171819202122232425'&
          //'26272829303132333435363738394041424344454647484950'
      pasoaux(1)=0
      pasoaux(2)=0
      pasoaux(3)=0
      pasoaux(4)=0
      pasoaux(5)=0
      pasoaux(6)=0
      pasoaux(7)=0
      pasoaux(8)=0
      pasoaux(9)=0
      pasoaux(10)=0
      pasoaux(11)=0
      pasoaux(12)=0
      pasoaux(13)=0
      pasoaux(14)=0
      pasoaux(15)=1
      pasoaux(16)=1
      pasoaux(17)=1
      pasoaux(18)=1
      pasoaux(19)=0
      pasoaux(20)=0
      pasoaux(21)=0
      pasoaux(22)=0
      pasoaux(23)=0
      pasoaux(24)=0
      pasoaux(25)=0
      !*    lectura de los datos generales y de base
      directory='Data/py_data/'

      open (newunit=unit_number, file=directory//"inis.da")
      read(unit_number,*) air_density_z_initial, temperature_z_initial, theta_z_initial, &
         Pres00, vapor_z_initial, cc2, aerosol_z_initial, u_z_initial, v_z_initial

      close (unit_number)

      open(newunit=unit_number,file=directory//'varconz.da',status='unknown',form='unformatted')
      read(unit_number)  Tvis,Tlvl,Tlsl,Tlvs,Telvs,Tesvs,Av,Vtnie,Vtgra0,&
         vapor_z_relative,aerosol_z_relative,Eautcn,Eacrcn
      close(unit_number)

      open(newunit=unit_number, file=directory//'cortes/possy')
      do 10 ii=1,25
         if(pasoaux(ii).eq.1) then
            open(newunit=unit_number2,file=directory//'nube'&
               //paso(ii*2-1:ii*2)//'.sal',status='unknown',form='unformatted')
            read(unit_number2) u_perturbed_base,v_perturbed_base,w_perturbed_base,theta_base,pressure_base,vapor_base,drop_base,&
               rain_base,crystal_base,snow_base,hail_base,aerosol_base
            close(unit_number2)

            !*        determinacion del maximo Y de la precipitacion
            ymax=0.
            maypre=0.
            do 300 i=1,nx1
               do 300 j=1,nx1
                  do 300 k=1,mod_nz1
                     if (maypre.lt.hail_base(i,j,k)+rain_base(i,j,k)+drop_base(i,j,k)) then
                        maypre=hail_base(i,j,k)+rain_base(i,j,k)+drop_base(i,j,k)
                        ymax=j
                     endif
300         continue

            write(unit_number,*) ii,ymax,maypre

            !*        cortes (no incluyen el piso)
            j=ymax
            do 50 i=1,nx1
               do 50 k=1,mod_nz1
                  !           calculo de la temperatua
                  aux=Pres00(k)+pressure_base(i,j,k)
                  variable(i,k,4)=(theta_z_initial(k)+theta_base(i,j,k))*aux*1e2

                  !           velocidades
                  variable(i,k,1)=(u_z_initial(k)+u_perturbed_base(i,j,k))*1e3
                  variable(i,k,2)=(v_z_initial(k)+v_perturbed_base(i,j,k))*1e3
                  variable(i,k,3)=w_perturbed_base(i,j,k)*1e3

                  !           filtrado de gotitas
                  aux=0.
                  do 200 l=i-1,i+1
                     do 200 m=j-1,j+1
                        do 200 n=k-1,k+1
                           aux=aux+drop_base(l,m,n)*1e6
200               continue
                  variable(i,k,5)=aux/27.
                  !           filtrado de lluvia
                  aux=0.
                  do 210 l=i-1,i+1
                     do 210 m=j-1,j+1
                        do 210 n=k-1,k+1
                           aux=aux+rain_base(l,m,n)*1e6
210               continue
                  variable(i,k,6)=aux/27.*.5+rain_base(i,j,k)*1e6*.5
                  !           filtrado de cristales
                  aux=0.
                  do 220 l=i-1,i+1
                     do 220 m=j-1,j+1
                        do 220 n=k-1,k+1
                           aux=aux+crystal_base(l,m,n)*1e6
220               continue
                  variable(i,k,7)=aux/27.
                  !           filtrado de nieve
                  aux=0.
                  do 230 l=i-1,i+1
                     do 230 m=j-1,j+1
                        do 230 n=k-1,k+1
                           aux=aux+snow_base(l,m,n)*1e6
230               continue
                  variable(i,k,8)=aux/27.
                  !           filtrado de granizos
                  aux=0.
                  do 240 l=i-1,i+1
                     do 240 m=j-1,j+1
                        do 240 n=k-1,k+1
                           aux=aux+hail_base(l,m,n)*1e6
240               continue
                  variable(i,k,9)=aux/27.*.5+hail_base(i,j,k)*1e6*.5

                  variable(i,k,10)=(vapor_z_initial(k)+vapor_base(i,j,k))*1e6
                  variable(i,k,11)=(aerosol_z_initial(k)+aerosol_base(i,j,k))*1e1
                  variable(i,k,12)=vapor_base(i,j,k)*1e6
                  variable(i,k,13)=aerosol_base(i,j,k)*1e1
50          continue

            do 80 lvar=1,nvar
               archaux=var(lvar)&
               //'y'//paso(ii*2-1:ii*2)
               ! write(*,*) archaux
               open(10+lvar,file=directory//'cortes/'//archaux)
               do 90 k=1,mod_nz1
                  write(10+lvar,1000) (variable(j,k,lvar),j=1,nx1)
90             continue
               close(10+lvar)
80          continue
         endif
10    continue
      close(unit_number)
1000  format(50i7)
      return
   end subroutine get_cutted
end module get_cut

