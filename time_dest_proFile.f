      PROGRAM TIME_DEST_PROFILE
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-k)
      
      
c      COMMON /ALFA(8)/
      COMMON CHILOMETRO,SOLARMASS,PARSEC,GGG,AU,SOLARLUM
c      PARAMETER (Nmax=10000)
     
      REAL*8 massTotSM, massTot,massMeanSM, massMean,PI,massStarTestSM,
     & massStarTest
c      REAL*8 rStep
      CHARACTER*8 nameCluster
      CHARACTER*12 buffer12
      CHARACTER*4  buffer4
 
 
C    THIS PROGRAM CALCULATES THE TIME FOR DESTABILIZING AN ORBIT OF A 
C    CANDIDATE LIFE HOST PLANET IN A GC 
       CHILOMETRO = 1.d5    ! cm
       SOLARMASS  = 1.99d33 ! gr
       PARSEC     = 3.09d18 ! cm
       GGG        = 6.67d-8 ! dyne cm**2/g
       SOLARLUM   = 3.9d33  ! erg/s
       AU         = 1.5d13  ! cm
       YEAR       = 3.2d7   ! s
1     FORMAT (A8,A12,f6.2,A12,f4.2,A12,f4.2) 
2     FORMAT (F6.2,A4,f4.2) 
       PI=4*atan(1.0)
C     Legge da file il nome del cluster e tutti gli altri parametri
      OPEN(UNIT=0, FILE='Dati_GCs.dat', ACTION='READ', STATUS='OLD', 
     & IOSTAT= KODE)
     
      IF (KODE .NE. 0) THEN 
      WRITE (*,*) 'Dati_GCs.dat cannot be opened'
      END IF
      
      READ(0,*) ! SKIP THE FIRST LINE
      READ(0,*) ! SKIP THE SECOND LINE
     
      READ(0, FMT=1) nameCluster,buffer12,massTotSM,buffer12,
     & radiusCorePC,buffer12,massMeanSM          
      WRITE(*,*) nameCluster,  massTotSM, radiusCorePC, massMeanSM
      CLOSE(0)
      
      radiusCore = radiusCorePC*PARSEC
      massTot    = massTotSM *SOLARMASS
      massMean   = massMeanSM *SOLARMASS
c     calcolo il raggio di cattura da Gliese 74
c     raggio orbitale =  0.04*AU
      massStarTestSM = 0.1 ! in solar masses: per ogni valore se 
c     questo parametro cambia deve cambiare anche orbRadius come 
c     in figura
      massStarTest = massStarTestSM*SOLARMASS
      orbRadius    = 0.04*AU
      dest = orbRadius * (2*massMeanSM/massStarTestSM)**(1./3.)
      
c calcolo il tempo di di destabilizzazione per il raggio del core
      
      sigmaQuad    = massTot*GGG/(2*radiusCore)
      sigma  = SQRT(sigmaQuad)
      densityCore = 3*massTotSM/(8*massMeanSM*PI*radiusCore**3)
      focusingG  = (dest**2 + (GGG*massStarTest*dest)/
     & sigmaQuad)
      time = 1/(4* SQRT(PI) * sigma *densityCore*focusingG) 
      timeYear = time/Year 
      
      WRITE(*,*) "massStarTestSM in g =", massStarTestSM*SOLARMASS
      WRITE(*,*) "orbital Radius in cm =", orbRadius 
      WRITE(*,*) "destabilization radius in cm =", dest 
      WRITE(*,*) "mean velocity in cm/s = ", sigma
      WRITE(*,*) "stellar density mean =", densityCore 
      WRITE(*,*) "time in yr =", timeYear 
      
c calcolo un profilo di tempo per i tempi di destabilizzazione 
c con i dati forniti dal DataBase online
c https://people.smp.uq.edu.au/HolgerBaumgardt/globular/veldis.html
c messi in un apposito file
      
3     FORMAT (F6.2) 
      
      
      OPEN(UNIT=1, FILE='NGC_6397_r.dat', ACTION='READ', STATUS='OLD',
     & IOSTAT= KODE0 )
      OPEN(UNIT=2, FILE='NGC_6397_v.dat', ACTION='READ', STATUS='OLD',
     & IOSTAT= KODE1 )
      READ(1,*) ! SKIP THE FIRST TWO LINES IN BOTH FILES
      READ(1,*) 
      READ(2,*) 
      READ(2,*) 
      READ (1, FMT =3) RSTEP0
      READ (2, FMT =3) VSTEP0
      WRITE(*,*) RSTEP0, VSTEP0, KODE0, KODE1
   
       
11    DO       
      
      READ (1, FMT=3) RSTEP
      READ (2, FMT=3) VSTEP
      IF(RSTEP == 0 .OR. VSTEP == 0) GO TO 12
      WRITE(*,*) RSTEP, VSTEP, KODE0, KODE1
      RSTEP0 = RSTEP
      VSTEP0 = VSTEP
      IF(KODE0 .LT. 1 .AND. KODE1 .LT. 1) GO TO 11
      END DO

12    CONTINUE
      
      CLOSE(1)
      CLOSE(2)

      
      END PROGRAM TIME_DEST_PROFILE
