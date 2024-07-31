      PROGRAM time_dest_proFile
      COMMON /ALFA(8)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      PARAMETER (Nmax=10000)
      REAL*8 massTot,massMean,radiusCore,pi
      REAL*8 rStep
      CHARACTER*16 nameCluster
 
 
C    THIS PROGRAM CALCULATES THE TIME FOR DESTABILIZING AN ORBIT OF A 
C    CANDIDATE LIFE HOST PLANET IN A GC 
       Chilometro = 1.d5    ! cm
       solarMass  = 1.99d33 ! gr
       parsec     = 3.09d18 ! cm
       GGG        = 6.67d-8 ! dyne cm**2/g
       SolarLum   = 3.9d33  ! erg/s
  
       pi=4*atan(1.0)
C     Legge da file il nome del cluster e tutti gli altri parametri
      OPEN(UNIT=0, FILE='Dati_GCs.dat', STATUS=OLD)
      
      
      
