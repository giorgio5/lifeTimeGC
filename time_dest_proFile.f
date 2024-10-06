      PROGRAM TIME_DEST_PROFILE
      IMPLICIT REAL*8 (A-H,L-Z)
      IMPLICIT INTEGER (I-K)
      
      COMMON /A/ CHILOMETRO,SOLARMASS,GGG,AU,YEAR
      COMMON /B/ DISTANCESUN,PARSEC

      CHARACTER*8 NAMECLUSTER
      CHARACTER*12 BUFFER12
      CHARACTER*4  BUFFER4
 
 
C    THIS PROGRAM CALCULATES THE TIME FOR DESTABILIZING AN ORBIT OF A 
C    CANDIDATE LIFE HOST PLANET IN A GC 
       CHILOMETRO = 1.d5    ! cm
       SOLARMASS  = 1.99d33 ! gr
       PARSEC     = 3.09d18 ! cm
       GGG        = 6.67d-8 ! dyne cm**2/g
       AU         = 1.5d13  ! cm
       YEAR       = 3.2d7   ! s
1     FORMAT (A8,A12,F6.2,A12,F4.2,A12,F3.2,A12,F4.2,A12,F5.2) 
2     FORMAT (F6.2,A4,F4.2) 
       PI=4*ATAN(1.0)
C     LEGGE DA FILE IL NOME DEL CLUSTER E TUTTI GLI ALTRI PARAMETRI

    
      OPEN(UNIT=0, FILE='Dati_GCs.dat', ACTION='READ', STATUS='OLD', 
     & IOSTAT= KODE0)
     
      IF (KODE0 .NE. 0) THEN 
      WRITE (*,*) 'Dati_GCs.dat CANNOT BE OPENED'
      END IF
      
      READ(0,*) ! SKIP THE FIRST LINE
      READ(0,*) ! SKIP THE SECOND LINE
     
      READ(0, FMT=1) NAMECLUSTER,BUFFER12,MASSTOTSM,BUFFER12,
     & RADIUSCOREPC,BUFFER12,MASSMEANSM,BUFFER12,DISTANCESUN,
     & BUFFER12,DELTA_SIS_DIST       
      WRITE(*,*) "namecluster ",NAMECLUSTER
      WRITE(*,*) "total mass in solar masses =",  MASSTOTSM
      WRITE(*,*) "Core radius in pc = ", RADIUSCOREPC
      WRITE(*,*) "mean mass in solar masses =",MASSMEANSM
      WRITE(*,*) "distance from sun in kpc =",DISTANCESUN
      WRITE(*,*) "delta_dist_from_Sun in pc =",DELTA_SIS_DIST
      CLOSE(0)
      
      RADIUSCORE = RADIUSCOREPC*PARSEC
      MASSTOT    = MASSTOTSM *SOLARMASS
      MASSMEAN   = MASSMEANSM *SOLARMASS
C     CALCOLO IL RAGGIO DI CATTURA DA GLIESE 74
C     RAGGIO ORBITALE =  0.04*AU
      MASSSTARTESTSM = 0.1 ! IN SOLAR MASSES: PER OGNI VALORE SE 
C     QUESTO PARAMETRO CAMBIA DEVE CAMBIARE ANCHE ORBRADIUS COME 
C     IN FIGURA
      MASSSTARTEST = MASSSTARTESTSM*SOLARMASS
      ORBRADIUS    = 0.04*AU
      DEST = ORBRADIUS * (2*MASSMEANSM/MASSSTARTESTSM)**(1./3.)
      
C CALCOLO IL TEMPO DI DI DESTABILIZZAZIONE PER IL RAGGIO DEL CORE
      
      SIGMAQUAD    = MASSTOT*GGG/(6*RADIUSCORE)
      SIGMA  = SQRT(SIGMAQUAD)
      DENSITYCORE = 3*MASSTOTSM/(8*MASSMEANSM*PI*RADIUSCORE**3)
      FOCUSINGG  = (DEST**2 + (GGG*MASSSTARTEST*DEST)/
     & SIGMAQUAD)
      TIME = 1/(4*SQRT(PI)*SIGMA*DENSITYCORE*FOCUSINGG) 
      TIMEYEAR = TIME/YEAR 
      
      WRITE(*,*) "massstartestsm in g =", MASSSTARTESTSM*SOLARMASS
      WRITE(*,*) "orbital radius in cm =", ORBRADIUS 
      WRITE(*,*) "destabilization radius in cm =", DEST 
      WRITE(*,*) "mean velocity in cm/s = ", SIGMA
      WRITE(*,*) "stellar density mean in (cm)^-3 =", DENSITYCORE 
      WRITE(*,*) "time in yr =", TIMEYEAR 
      
C CALCOLO UN PROFILO DI TEMPO PER I TEMPI DI DESTABILIZZAZIONE 
C CON I DATI FORNITI DAL DATABASE ONLINE
C https://people.smp.uq.edu.au/HolgerBaumgardt/globular/veldis.html
C MESSI IN UN APPOSITO FILE
      
3     FORMAT (F6.2) 
      
      
      OPEN(UNIT=1, FILE='NGC_6397_r.dat', ACTION='READ', STATUS='OLD',
     & IOSTAT= KODE1 )
      OPEN(UNIT=2, FILE='NGC_6397_v.dat', ACTION='READ', STATUS='OLD',
     & IOSTAT= KODE2 )
      OPEN(UNIT=3, FILE='NGC_6397_error_upper_v.dat', ACTION='READ',
     & STATUS='OLD',IOSTAT= KODE3 )
      OPEN(UNIT=4, FILE='NGC_6397_error_below_v.dat', ACTION='READ', 
     & STATUS='OLD',IOSTAT= KODE4 )
      READ(1,*) ! SKIP THE FIRST TWO LINES IN EVERY FILE
      READ(1,*) 
      READ(2,*) 
      READ(2,*) 
      READ(3,*) 
      READ(3,*) 
      READ(4,*) 
      READ(4,*) 
      
11    DO       
      
      READ (1, FMT=3) R
      READ (2, FMT=3) V
      READ (3, FMT =3) DELTA_V_UP
      READ (4, FMT =3) DELTA_V_LOW
      WRITE(*,*) "r (arcsec)= ", R
      WRITE(*,*) "v (km/s)= ", V
      IF(R == 0 .OR. V == 0 .OR. DELTA_V_UP == 0 .OR. DELTA_V_LOW == 0)
     & GO TO 12
      V_CM = CHILOMETRO*V
      DELTA_V_UP_CM = CHILOMETRO*DELTA_V_UP
      DELTA_V_LOW_CM = CHILOMETRO*DELTA_V_LOW
      R_CM = R_ARCSEC_2_CM(R)
      WRITE(*,*) "r (cm)= ", R_CM
      WRITE(*,*) "v (cm/s)= ", V_CM
      MASSA = (R_CM)*6*(V_CM**2)/GGG
      DELTA_MASSA_UP = 2*MASSA*DELTA_V_UP_CM/(V_CM)
      DELTA_MASSA_LOW = 2*MASSA*DELTA_V_LOW_CM/(V_CM)
      DELTA_MASSA_SIS = MASSA*DELTA_SIS_DIST/(DISTANCESUN*1000)
      WRITE(*,*) "massa_r (g)= ", MASSA
     
      
      IF(KODE0 .LT. 1 .AND. KODE1 .LT. 1 .AND. KODE2 .LT. 1 
     & .AND. KODE3 .LT. 1) GO TO 11
      END DO

12    CONTINUE

      
      WRITE(*,*) "massa_tot (solarmasses)= ", MASSA/SOLARMASS
      WRITE(*,*) "delta_m_up(solarmasses)= ", DELTA_MASSA_UP/SOLARMASS
      WRITE(*,*) "delta_m_low(solarmasses)= ", DELTA_MASSA_LOW/SOLARMASS
      WRITE(*,*) "delta_m_sis(solarmasses)= ", DELTA_MASSA_SIS/SOLARMASS
      WRITE(*,*) "massa_tot_prec (solarmasses)= ", MASSTOTSM
      
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)

      
      END PROGRAM TIME_DEST_PROFILE
    
C     CONVERTO IL RAGGIO DA ARCOSECONDI A PARSEC E POI IN CENTIMETRI
C     CON UNA FUNZIONE
 
      REAL*8 FUNCTION R_ARCSEC_2_CM(R_ARC)
      REAL*8 R_ARC, DISTANCESUN, PARSEC
      COMMON /B/ DISTANCESUN, PARSEC
      R_ARCSEC_2_CM = R_ARC*DISTANCESUN*1000*PARSEC/206265 
      RETURN 
      END 
      
      
      
      
       



