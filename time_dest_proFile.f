      program time_dest_profile
      implicit real*8 (a-h,l-z)
      implicit integer (i-k)
      
      common /a/ chilometro,solarmass,ggg,au,year, pi
      common /b/ distancesun, parsec

      character*8 namecluster
      character*12 buffer12
      character*4  buffer4
 
 
c    this program calculates the time for destabilizing an orbit of a 
c    candidate life host planet in a gc 

       chilometro = 1.d5    ! cm
       solarmass  = 1.99d33 ! gr
       parsec     = 3.09d18 ! cm
       ggg        = 6.67d-8 ! dyne cm**2/g
       au         = 1.5d13  ! cm
       year       = 3.2d7   ! s
       pi=4*atan(1.0)
       
1     format (a8,a12,f6.2,a12,f4.2,a12,f3.2,a12,f4.2,a12,f5.2) 
2     format (f6.2,a4,f4.2) 
       
       
c     read from file name cluster and other parametrers proper of the 
c     cluster itself

    
      open(unit=0, file='Dati_GCs.dat', action='read', status='old', 
     & iostat= kode0)
     
      if (kode0 .ne. 0) then 
      write (*,*) 'Dati_GCs.dat cannot be opened'
      end if
      
      read(0,*) ! skip the first line
      read(0,*) ! skip the second line
     
      read(0, fmt=1) namecluster,buffer12,masstotsm,buffer12,
     & radiuscorepc,buffer12,massmeansm,buffer12,distancesun,
     & buffer12,delta_sis_dist       
      write(*,*) "namecluster ",namecluster
      write(*,*) "total mass in solar masses =",  masstotsm
      write(*,*) "core radius in pc = ", radiuscorepc
      write(*,*) "mean mass in solar masses =",massmeansm
      write(*,*) "distance from sun in kpc =",distancesun
      write(*,*) "delta_dist_from_sun in pc =",delta_sis_dist
      close(0)
      
      radiuscore = radiuscorepc*parsec
      masstot    = masstotsm *solarmass
      massmean   = massmeansm *solarmass
c     calculus of the capture radius of gliese 74
c     capture radius is 0.04 AU
      massstartestsm = 0.1 ! in solar masses
c     this parameter is going to change as in figure
      massstartest = massstartestsm*solarmass
      orbradius    = 0.04*au
      dest = orbradius * (2*massmeansm/massstartestsm)**(1./3.)
      
c calcolo il tempo di di destabilizzazione per il raggio del core
      
      sigmaquad    = masstot*ggg/(6*radiuscore)
      sigma  = sqrt(sigmaquad)
      densitycore = 3*masstotsm/(8*massmeansm*pi*radiuscore**3)
      focusingg  = (dest**2 + (ggg*massstartest*dest)/
     & sigmaquad)
      time = 1/(4*sqrt(pi)*sigma*densitycore*focusingg) 
      timeyear = time/year 
      
      write(*,*) "massstartest in g =", massstartestsm*solarmass
      write(*,*) "orbital radius in cm =", orbradius 
      write(*,*) "destabilization radius in cm =", dest 
      write(*,*) "mean velocity in cm/s = ", sigma
      write(*,*) "stellar density mean in (cm)^-3 =", densitycore 
      write(*,*) "time in yr =", timeyear 
      
c calculus of the time of destabilization of an orbital planet traiectory
c from the data presents on the online database 
c https://people.smp.uq.edu.au/HolgerBaumgardt/globular/veldis.html
c put in the files: Dati_GCs.dat and NGC_0000_r.dat and NGC_0000_v.dati
c and so on for the error bands
      
3     format (f6.2) 
      
      
      open(unit=1, file='NGC_6397_r.dat', action='read', status='old',
     & iostat= kode1 )
      open(unit=2, file='NGC_6397_v.dat', action='read', status='old',
     & iostat= kode2 )
      open(unit=3, file='NGC_6397_error_upper_v.dat', action='read',
     & status='old',iostat= kode3 )
      open(unit=4, file='NGC_6397_error_below_v.dat', action='read', 
     & status='old',iostat= kode4 )
      read(1,*) ! skip the first two lines in every file
      read(1,*) 
      read(2,*) 
      read(2,*) 
      read(3,*) 
      read(3,*) 
      read(4,*) 
      read(4,*) 
      open(unit=5, file='NGC_6397_lifetime_year.dat', action='write', 
     & status='replace',iostat= kode5 )
     
11    do       
      
      read (1, fmt=3) r
      read (2, fmt=3) v
      read (3, fmt =3) delta_v_up
      read (4, fmt =3) delta_v_low
      write(*,*) "r (arcsec)= ", r
      write(*,*) "v (km/s)= ", v
      if(r == 0 .or. v == 0 .or. delta_v_up == 0 .or. delta_v_low == 0)
     & go to 12
      v_cm = chilometro*v
      delta_v_up_cm = chilometro*delta_v_up
      delta_v_low_cm = chilometro*delta_v_low
      r_cm = r_arcsec_2_cm(r)
      write(*,*) "r (cm)= ", r_cm
      write(*,*) "v (cm/s)= ", v_cm
      massa = (r_cm)*6*(v_cm**2)/ggg
      
      call  calculus_time_r(r_cm, v_cm, massa, dest, massmean, 
     & massstartest)
      
      delta_massa_up = 2*massa*delta_v_up_cm/(v_cm)
      delta_massa_low = 2*massa*delta_v_low_cm/(v_cm)
      delta_massa_sis = massa*delta_sis_dist/(distancesun*1000)
      write(*,*) "massa_r (g)= ", massa
     
       
      if(kode0 .lt. 1 .and. kode1 .lt. 1 .and. kode2 .lt. 1 
     & .and. kode3 .lt. 1) go to 11
      end do

12    continue

      
      write(*,*) "massa_tot (solarmasses)= ", massa/solarmass
      write(*,*) "delta_m_up(solarmasses)= ", delta_massa_up/solarmass
      write(*,*) "delta_m_low(solarmasses)= ", delta_massa_low/solarmass
      write(*,*) "delta_m_sis(solarmasses)= ", delta_massa_sis/solarmass
      write(*,*) "massa_tot_prec (solarmasses)= ", masstotsm
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      
      end program time_dest_profile
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     functions and subroutines


c     function to convert radius from arcsec to parsec and then to cm
 
      real*8 function r_arcsec_2_cm(r_arc)
      real*8 r_arc, distancesun, parsec
      common /b/ distancesun, parsec
      r_arcsec_2_cm = r_arc*distancesun*1000*parsec/206265 
      return 
      end 
      
c     subroutine for the calculus of the lifetime destabilization at arbitrary r
      
      subroutine calculus_time_r (r_cm, v_cm, massa, dest, massmean, 
     & massstartest)
     
      implicit real*8 (a-h,l-z)
      implicit integer (i-k)
      
      common /a/ chilometro,solarmass,ggg,au,year, pi
      common /b/ distancesun, parsec
4     format (f6.2,a4,f18.2) 
      
      
      r = r_cm
      v_r = v_cm
      m_r = massa
      sigmaquad_r = v_r**2
      density_r = 3*m_r/(8*massmean*pi*r**3)
      focusingg_r  = (dest**2 + (ggg*massstartest*dest)/
     & sigmaquad_r)
      
      time_r = 1/(4*sqrt(pi)*v_r*density_r*focusingg_r) 
      timeyear_r = time_r/year 
      
      write (5,fmt =4) r/parsec,'    ', timeyear_r 
      
      end 
      
       



