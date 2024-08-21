      PROGRAM TIME_DEST_PROFILE
c      IMPLICIT REAL*8 (a-h,o-z)
c      IMPLICIT INTEGER (i-n)
c      COMMON /ALFA(8)/
      COMMON CHILOMETRO,SOLARMASS,PARSEC,GGG,SOLARLUM
c      PARAMETER (Nmax=10000)
      REAL*8 massTot,massMean,radiusCore,PI
c      REAL*8 rStep
      CHARACTER*8 nameCluster
      CHARACTER*12 buffer12
 
 
C    THIS PROGRAM CALCULATES THE TIME FOR DESTABILIZING AN ORBIT OF A 
C    CANDIDATE LIFE HOST PLANET IN A GC 
       CHILOMETRO = 1.d5    ! cm
       SOLARMASS  = 1.99d33 ! gr
       PARSEC     = 3.09d18 ! cm
       GGG        = 6.67d-8 ! dyne cm**2/g
       SOLARLUM   = 3.9d33  ! erg/s
1 	   FORMAT (A8,A12,f6.2,A12,f4.2,A12,f4.2) 

       PI=4*atan(1.0)
C     Legge da file il nome del cluster e tutti gli altri parametri
      OPEN(UNIT=0, FILE='Dati_GCs.dat', ACTION='READ', STATUS='OLD', 
     & IOSTAT= KODE)
     
      IF (KODE .NE. 0) THEN 
      WRITE (*,*) 'Dati_GCs.dat cannot be opened'
      END IF
      
      READ(0,*) ! SKIP THE FIRST LINE
      READ(0,*) ! SKIP THE SECOND LINE
     

      
      READ(0, FMT=1) nameCluster,buffer,massTot,buffer,radiusCore,
     & buffer,massMean          
      WRITE(*,*) nameCluster,  massTot, radiusCore, massMean
      CLOSE(0)
      
      

        
      
      
      
      
      END PROGRAM TIME_DEST_PROFILE
