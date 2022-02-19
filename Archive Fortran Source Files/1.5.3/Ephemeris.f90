module Ephemeris
    
    use MathAndProgConstants
    
    private :: FSIZER1, FSIZER2, FSIZER3, PLEPH, INTERP, SPLIT, STATE, CONST_EPH
    public :: PlanetState
    
    contains
    
        ! About function, opened for F2PY wrapping
        function ephemeris_about() result(i)
            integer :: i
            i = 0
            write(*,*) 'Ephemeris: Version 1.0.'
            write(*,*) 'Author: Maksim Shirobokov.'
            write(*,*) 'Date: 19.01.2022.'
        end function ephemeris_about
    
        ! Opened for F2PY wrapping
        function PlanetState(JulianDateDays,TargetPlanet,CenterPlanet) result(x)
        
            ! DE430 in JPLEPH
        
            implicit none
        
            real(dp) :: JulianDateDays
            character(*) :: TargetPlanet, CenterPlanet
            real(dp), dimension(6) :: x
            integer :: TargetPlanetIndex, CenterPlanetIndex
            
            if (TargetPlanet .eq. 'Earth') then
                TargetPlanetIndex = 3
            elseif (TargetPlanet .eq. 'Moon') then
                TargetPlanetIndex = 10
            elseif (TargetPlanet .eq. 'Venus') then
                TargetPlanetIndex = 2
            elseif (TargetPlanet .eq. 'Jupiter') then  
                TargetPlanetIndex = 5
            elseif (TargetPlanet .eq. 'Saturn') then 
                TargetPlanetIndex = 6
            elseif (TargetPlanet .eq. 'Mercury') then 
                TargetPlanetIndex = 1
            elseif (TargetPlanet .eq. 'Mars') then
                TargetPlanetIndex = 4
            elseif (TargetPlanet .eq. 'Uranus') then
                TargetPlanetIndex = 7
            elseif (TargetPlanet .eq. 'Neptune') then
                TargetPlanetIndex = 8
            elseif (TargetPlanet .eq. 'Pluto') then
                TargetPlanetIndex = 9
            elseif (TargetPlanet .eq. 'Sun') then
                TargetPlanetIndex = 11
            elseif (TargetPlanet .eq. 'SolarSystemBarycenter') then
                TargetPlanetIndex = 12
            elseif (TargetPlanet .eq. 'EarthMoonBarycenter') then
                TargetPlanetIndex = 13
            else
                write(*,*) 'Unknown TargetPlanet'
                stop
            endif
            
            if (CenterPlanet .eq. 'Earth') then
                CenterPlanetIndex = 3
            elseif (CenterPlanet .eq. 'Moon') then
                CenterPlanetIndex = 10
            elseif (CenterPlanet .eq. 'Venus') then
                CenterPlanetIndex = 2
            elseif (CenterPlanet .eq. 'Jupiter') then  
                CenterPlanetIndex = 5
            elseif (CenterPlanet .eq. 'Saturn') then 
                CenterPlanetIndex = 6
            elseif (CenterPlanet .eq. 'Mercury') then 
                CenterPlanetIndex = 1
            elseif (CenterPlanet .eq. 'Mars') then
                CenterPlanetIndex = 4
            elseif (CenterPlanet .eq. 'Uranus') then
                CenterPlanetIndex = 7
            elseif (CenterPlanet .eq. 'Neptune') then
                CenterPlanetIndex = 8
            elseif (CenterPlanet .eq. 'Pluto') then
                CenterPlanetIndex = 9
            elseif (CenterPlanet .eq. 'Sun') then
                CenterPlanetIndex = 11
            elseif (CenterPlanet .eq. 'SolarSystemBarycenter') then
                CenterPlanetIndex = 12
            elseif (CenterPlanet .eq. 'EarthMoonBarycenter') then
                CenterPlanetIndex = 13
            else
                write(*,*) 'Unknown CenterPlanet'
                stop
            endif
            
            call pleph(JulianDateDays,TargetPlanetIndex,CenterPlanetIndex,x)
            x(1:3) = x(1:3)*149597870.700D0 ! *AU
            x(4:6) = x(4:6)*1.731456836805555D+03 ! *AUpDAY
            
        
        end function PlanetState
    
        ! Opened for F2PY wrapping
        function MoonLibration(JulianDateDays) result(angles)
        
            ! DE430 in JPLEPH
        
            implicit none
            
            real(dp) :: JulianDateDays
            real(dp) :: angles(3)
            real(dp) :: pleh_output(6)
            
            call pleph(JulianDateDays,15,0,pleh_output)
            angles = pleh_output(1:3)
        
        end function MoonLibration
        
        ! Opened for F2PY wrapping
        function JulianDate(Year,Month,Day,Hour,Min,Sec) result(JD)
        
            implicit none
        
            integer :: Year, Month, Day, Hour, Min, Sec
            integer :: YearCorr, MonthCorr
            real(dp) :: PrecDay, JD, dayFraction
        
            YearCorr = Year
            MonthCorr = Month
            if (Month <= 2) then
                YearCorr = Year - 1
                MonthCorr = Month + 12
            endif
            
            dayFraction = (Hour + Min/60.0D0 + Sec/3600.0D0)/24.0D0
            
            PrecDay = floor(365.25*(YearCorr+4716.0)) + &
                      floor(30.6001*(MonthCorr+1.0)) + 2.0 - &
                      floor(YearCorr/100.0) + &
                      floor(floor(YearCorr/100.0)/4.0) + Day - 1524.5
            
            JD = dayFraction + PrecDay
            
        end function JulianDate
        
        SUBROUTINE FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)

            !   Version 1.0 uses the INQUIRE statement to find out the the record length 
            !   of the direct access file before opening it.  This procedure is non-standard, 
            !   but seems to work for VAX machines. 

            !  THE SUBROUTINE ALSO SETS THE VALUES OF  NRECL, NRFILE, AND NAMFIL.

            !  *****************************************************************
            !  *****************************************************************
            !
            !  THE PARAMETERS NAMFIL, NRECL, AND NRFILE ARE TO BE SET BY THE USER
            !
            !  *****************************************************************

            !  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

            CHARACTER*80  NAMFIL

            NAMFIL='JPLEPH'

            !  *****************************************************************

            !  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
            !  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
            !  (for a VAX, it is probably 1)
            !
            NRECL= 4

            !  *****************************************************************

            !  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE

            NRFILE=12

            !  *****************************************************************

            !   FIND THE RECORD SIZE USING THE INQUIRE STATEMENT  


            !      IRECSZ=0

            INQUIRE(FILE=NAMFIL,RECL=IRECSZ)

            ! IF 'INQUIRE' DOES NOT WORK, USUALLY IRECSZ WILL BE LEFT AT 0 

            IF(IRECSZ .LE. 0) THEN
                write(*,*) ' INQUIRE STATEMENT PROBABLY DID NOT WORK'
            ENDIF

            KSIZE = IRECSZ/NRECL

            RETURN

        END SUBROUTINE

        SUBROUTINE FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)

            !  THIS SUBROUTINE OPENS THE FILE, 'NAMFIL', WITH A PHONY RECORD LENGTH, READS 
            !  THE FIRST RECORD, AND USES THE INFO TO COMPUTE KSIZE, THE NUMBER OF SINGLE 
            !  PRECISION WORDS IN A RECORD.  

            !  THE SUBROUTINE ALSO SETS THE VALUES OF  NRECL, NRFILE, AND NAMFIL.

            IMPLICIT DOUBLE PRECISION(A-H,O-Z)

            SAVE

            INTEGER OLDMAX
            PARAMETER (OLDMAX = 400)
            INTEGER NMAX
            PARAMETER (NMAX = 1000)

            CHARACTER*6 TTL(14,3),CNAM(NMAX)
            CHARACTER*80 NAMFIL

            DIMENSION SS(3)

            INTEGER IPT(3,13)

            !  *****************************************************************
            !  *****************************************************************

            !  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER
            !
            !  *****************************************************************

            !  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
            !  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
            !  (for UNIX, it is probably 4)

            NRECL=4

            !  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE

            NRFILE=12

            !  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

            NAMFIL='JPLEPH' 

            !  *****************************************************************
            !  *****************************************************************

            !  **  OPEN THE DIRECT-ACCESS FILE AND GET THE POINTERS IN ORDER TO 
            !  **  DETERMINE THE SIZE OF THE EPHEMERIS RECORD

            MRECL = NRECL*1000

            OPEN(NRFILE,FILE=NAMFIL,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=MRECL,STATUS='OLD')

            READ(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT,&
                ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)

            CLOSE(NRFILE)

            !  FIND THE NUMBER OF EPHEMERIS COEFFICIENTS FROM THE POINTERS

            KMX = 0
            KHI = 0

            DO I = 1,13
                IF (IPT(1,I) .GT. KMX) THEN
                    KMX = IPT(1,I)
                    KHI = I
                ENDIF
            ENDDO

            ND = 3
            IF (KHI .EQ. 12) ND=2

            KSIZE = 2*(IPT(1,KHI)+ND*IPT(2,KHI)*IPT(3,KHI)-1)

            RETURN

        END SUBROUTINE

        SUBROUTINE FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

            !  THE SUBROUTINE SETS THE VALUES OF  NRECL, KSIZE, NRFILE, AND NAMFIL.

            SAVE

            CHARACTER*80 NAMFIL

            !  *****************************************************************
            !  *****************************************************************
            !
            !  THE PARAMETERS NRECL, NRFILE, AND NAMFIL ARE TO BE SET BY THE USER

            !  *****************************************************************

            !  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
            !  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
      
            NRECL=4

            !  *****************************************************************

            !  NRFILE IS THE INTERNAL UNIT NUMBER USED FOR THE EPHEMERIS FILE (DEFAULT: 12)

            NRFILE=12

            !  *****************************************************************

            !  NAMFIL IS THE EXTERNAL NAME OF THE BINARY EPHEMERIS FILE

            NAMFIL='JPLEPH'

            !  *****************************************************************

            !  KSIZE must be set by the user according to the ephemeris to be read

            !  For  de200, set KSIZE to 1652
            !  For  de405, set KSIZE to 2036
            !  For  de406, set KSIZE to 1456
            !  For  de414, set KSIZE to 2036
            !  For  de418, set KSIZE to 2036
            !  For  de421, set KSIZE to 2036
            !  For  de422, set KSIZE to 2036
            !  For  de423, set KSIZE to 2036
            !  For  de424, set KSIZE to 2036
            !  For  de430, set KSIZE to 2036

            KSIZE = 2036

            !  *******************************************************************

            RETURN

        END SUBROUTINE

        SUBROUTINE PLEPH ( ET, NTARG, NCENT, RRD )

            !  NOTE : Over the years, different versions of PLEPH have had a fifth argument:
            !  sometimes, an error return statement number; sometimes, a logical denoting
            !  whether or not the requested date is covered by the ephemeris.  We apologize
            !  for this inconsistency; in this present version, we use only the four necessary 
            !  arguments and do the testing outside of the subroutine.

            !     THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS
            !     AND GIVES THE POSITION AND VELOCITY OF THE POINT 'NTARG'
            !     WITH RESPECT TO 'NCENT'.
            !
            !     CALLING SEQUENCE PARAMETERS:
            !
            !       ET = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION
            !            IS WANTED.
            !
            !       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **
            !          THE REASON FOR THIS OPTION IS DISCUSSED IN THE 
            !          SUBROUTINE STATE
            !
            !     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
            !
            !     NCENT = INTEGER NUMBER OF CENTER POINT.
            !
            !            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
            !
            !                1 = MERCURY           8 = NEPTUNE
            !                2 = VENUS             9 = PLUTO
            !                3 = EARTH            10 = MOON
            !                4 = MARS             11 = SUN
            !                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
            !                6 = SATURN           13 = EARTH-MOON BARYCENTER
            !                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
            !                            15 = LIBRATIONS, IF ON EPH FILE
            !
            !             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
            !              SET NTARG = 15. SET NCENT=0.)
            !
            !      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
            !            OF POINT 'NTARG' RELATIVE TO 'NCENT'. THE UNITS ARE AU AND
            !            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
            !            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
            !            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
            !            RADIANS AND RADIANS/DAY.
            !
            !            The option is available to have the units in km and km/sec.
            !            For this, set km=.true. in the STCOMX common block.
            !

            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

            INTEGER NMAX
            PARAMETER (NMAX = 1000)

            DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
            DIMENSION PVST(6,11),PNUT(4)
            DIMENSION SS(3),CVAL(NMAX),PVSUN(6),ZIPS(2)
            DATA ZIPS/2*0.d0/

            LOGICAL BSAVE,KM,BARY
            LOGICAL FIRST
            DATA FIRST/.TRUE./

            INTEGER LIST(12),IPT(39),DENUM

            COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT

            COMMON/STCOMX/KM,BARY,PVSUN

            !     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT

            ET2(1)=ET
            ET2(2)=0.D0
            GO TO 11

            !     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT 
            !          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

            ENTRY DPLEPH(ET2Z,NTARG,NCENT,RRD)

            ET2(1)=ET2Z(1)
            ET2(2)=ET2Z(2)

  11        IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT)
            FIRST=.FALSE.

  96        IF(NTARG .EQ. NCENT) RETURN

            DO I=1,12
                LIST(I)=0
            ENDDO

            !     CHECK FOR NUTATION CALL

            IF(NTARG.NE.14) GO TO 97
            IF(IPT(35).GT.0) THEN
                LIST(11)=2
                CALL STATE(ET2,LIST,PVST,PNUT)
                DO I=1,4
                    RRD(I)=PNUT(I)
                ENDDO
                RRD(5) = 0.d0
                RRD(6) = 0.d0
                RETURN
            ELSE
                DO I=1,4
                    RRD(I)=0.d0
                ENDDO
                WRITE(6,297)
  297           FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
                STOP
            ENDIF

            !     CHECK FOR LIBRATIONS

  97        CONTINUE
            DO I=1,6
                RRD(I)=0.d0
            ENDDO

            IF(NTARG.NE.15) GO TO 98
            IF(IPT(38).GT.0) THEN
                LIST(12)=2
                CALL STATE(ET2,LIST,PVST,PNUT)
                DO I=1,6
                    RRD(I)=PVST(I,11)
                ENDDO
                RETURN
            ELSE
                WRITE(6,298)
  298           FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
                STOP
            ENDIF

            !       FORCE BARYCENTRIC OUTPUT BY 'STATE'

  98        BSAVE=BARY
            BARY=.TRUE.

            !       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL

            DO I=1,2
                K=NTARG
                IF(I .EQ. 2) K=NCENT
                IF(K .LE. 10) LIST(K)=2
                IF(K .EQ. 10) LIST(3)=2
                IF(K .EQ. 3)  LIST(10)=2
                IF(K .EQ. 13) LIST(3)=2
            ENDDO

            !       MAKE CALL TO STATE

            CALL STATE(ET2,LIST,PVST,PNUT)

            DO I=1,10
                DO J = 1,6
                    PV(J,I) = PVST(J,I)
                ENDDO
            ENDDO

            IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
                DO I=1,6
                    PV(I,11)=PVSUN(I)
                ENDDO
            ENDIF

            IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
                DO I=1,6
                    PV(I,12)=0.D0
                ENDDO
            ENDIF

            IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
                DO I=1,6
                    PV(I,13) = PVST(I,3)
                ENDDO
            ENDIF

            IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
                DO I=1,6
                    PV(I,3)=0.D0
                ENDDO
                GO TO 99
            ENDIF

            IF(LIST(3) .EQ. 2) THEN
                DO I=1,6
                    PV(I,3)=PVST(I,3)-PVST(I,10)/(1.D0+EMRAT)
                ENDDO
            ENDIF

            IF(LIST(10) .EQ. 2) THEN
                DO I=1,6
                    PV(I,10) = PV(I,3)+PVST(I,10)
                ENDDO
            ENDIF

  99        DO I=1,6
                RRD(I)=PV(I,NTARG)-PV(I,NCENT)
            ENDDO

            BARY=BSAVE

            RETURN
            
        END SUBROUTINE

        SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)

            !     THIS SUBROUTINE DIFFERENTIATES AND INTERPOLATES A
            !     SET OF CHEBYSHEV COEFFICIENTS TO GIVE POSITION AND VELOCITY
            !
            !     CALLING SEQUENCE PARAMETERS:
            !
            !       INPUT:
            !
            !         BUF   1ST LOCATION OF ARRAY OF D.P. CHEBYSHEV COEFFICIENTS OF POSITION
            !
            !           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
            !               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
            !               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
            !               INTERVAL IN INPUT TIME UNITS.
            !
            !         NCF   # OF COEFFICIENTS PER COMPONENT
            !
            !         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
            !
            !          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
            !               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
            !
            !          IFL  INTEGER FLAG: =1 FOR POSITIONS ONLY
            !                             =2 FOR POS AND VEL
            !
            !
            !       OUTPUT:
            !
            !         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
            !               EXPECTED IS PV(NCM,IFL), DP.
        
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

            SAVE
            
            DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)

            DATA NP/2/
            DATA NV/3/
            DATA TWOT/0.D0/
            DATA PC(1),PC(2)/1.D0,0.D0/
            DATA VC(2)/1.D0/

            !       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
            !       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
            !       WITHIN THAT SUBINTERVAL.

            DNA=DBLE(NA)
            DT1=DINT(T(1))
            TEMP=DNA*T(1)
            L=IDINT(TEMP-DT1)+1

            !         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)

            TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0

            !       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
            !       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
            !       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
            !       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)

            IF(TC.NE.PC(2)) THEN
                NP=2
                NV=3
                PC(2)=TC
                TWOT=TC+TC
            ENDIF

            !       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
            !       AND ARE STORED IN THE ARRAY 'PC'.

            IF(NP.LT.NCF) THEN
                DO 1 I=NP+1,NCF
                PC(I)=TWOT*PC(I-1)-PC(I-2)
    1           CONTINUE
                NP=NCF
            ENDIF
            !
            !       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
            !
            DO 2 I=1,NCM
                PV(I,1)=0.D0
                DO 3 J=NCF,1,-1
                    PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
    3               CONTINUE
    2               CONTINUE
            IF(IFL.LE.1) RETURN
            !
            !       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
            !       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
            !
            VFAC=(DNA+DNA)/T(2)
            VC(3)=TWOT+TWOT
            IF(NV.LT.NCF) THEN
                DO 4 I=NV+1,NCF
                    VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
    4               CONTINUE
                NV=NCF
            ENDIF
            !
            !       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
            !
            DO 5 I=1,NCM
                PV(I,2)=0.D0
                DO 6 J=NCF,2,-1
                    PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
    6               CONTINUE
                PV(I,2)=PV(I,2)*VFAC
    5           CONTINUE
!
            RETURN
!
        END SUBROUTINE

        SUBROUTINE SPLIT(TT,FR)

            !     THIS SUBROUTINE BREAKS A D.P. NUMBER INTO A D.P. INTEGER
            !     AND A D.P. FRACTIONAL PART.
            !
            !     CALLING SEQUENCE PARAMETERS:
            !
            !       TT = D.P. INPUT NUMBER
            !
            !       FR = D.P. 2-WORD OUTPUT ARRAY.
            !            FR(1) CONTAINS INTEGER PART
            !            FR(2) CONTAINS FRACTIONAL PART
            !
            !            FOR NEGATIVE INPUT NUMBERS, FR(1) CONTAINS THE NEXT
            !            MORE NEGATIVE INTEGER; FR(2) CONTAINS A POSITIVE FRACTION.
            !
            !       CALLING SEQUENCE DECLARATIONS
            !
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

            DIMENSION FR(2)

            !       MAIN ENTRY -- GET INTEGER AND FRACTIONAL PARTS

            FR(1)=DINT(TT)
            FR(2)=TT-FR(1)

            IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN

            !       MAKE ADJUSTMENTS FOR NEGATIVE INPUT NUMBER

            FR(1)=FR(1)-1.D0
            FR(2)=FR(2)+1.D0

            RETURN

      END SUBROUTINE

        SUBROUTINE STATE(ET2,LIST,PV,PNUT)

            ! THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE
            !
            !     CALLING SEQUENCE PARAMETERS:
            !
            !     INPUT:
            !
            !         ET2   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
            !               IS WANTED.  ANY COMBINATION OF ET2(1)+ET2(2) WHICH FALLS
            !               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
            !
            !                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
            !                   ENTIRE EPOCH IN ET2(1) AND SET ET2(2)=0.
            !
            !                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET ET2(1) =
            !                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
            !                   EPOCH AND SET ET2(2) = FRACTIONAL PART OF A DAY
            !                   ELAPSED BETWEEN ET2(1) AND EPOCH.
            !
            !                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
            !                   ET2(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
            !                   AND ET2(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
            !
            !        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
            !               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
            !
            !                         LIST(I)=0, NO INTERPOLATION FOR BODY I
            !                                =1, POSITION ONLY
            !                                =2, POSITION AND VELOCITY
            !
            !               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
            !
            !                         I = 1: MERCURY
            !                           = 2: VENUS
            !                           = 3: EARTH-MOON BARYCENTER
            !                           = 4: MARS
            !                           = 5: JUPITER
            !                           = 6: SATURN
            !                           = 7: URANUS
            !                           = 8: NEPTUNE
            !                           = 9: PLUTO
            !                           =10: GEOCENTRIC MOON
            !                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
            !                           =12: LUNAR LIBRATIONS (IF ON FILE)
            !
            !     OUTPUT:
            !
            !          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
            !               QUANTITIES (OTHER THAN NUTATION, STOERD IN PNUT).  
            !               THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
            !               STATE IN THE ARRAY STARTING AT PV(1,I).  
            !               (ON ANY GIVEN CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE 
            !                AFFECTED BY THE  FIRST 10 'LIST' ENTRIES, AND BY LIST(12)
            !                IF LIBRATIONS ARE ON THE FILE, ARE SET.  
            !                THE REST OF THE 'PV' ARRAYIS UNTOUCHED.)  
            !               THE ORDER OF COMPONENTS STARTING IN PV(1,I) IS: X,Y,Z,DX,DY,DZ.
            !
            !               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
            !               EQUATOR AND EQUINOX OF J2000 IF THE DE NUMBER IS 200 OR
            !               GREATER; OF B1950 IF THE DE NUMBER IS LESS THAN 200. 
            !
            !               THE MOON STATE IS ALWAYS GEOCENTRIC; THE OTHER NINE STATES 
            !               ARE EITHER HELIOCENTRIC OR SOLAR-SYSTEM BARYCENTRIC, 
            !               DEPENDING ON THE SETTING OF COMMON FLAGS (SEE BELOW).
            !
            !               LUNAR LIBRATIONS, IF ON FILE, ARE PUT INTO PV(K,11) IF
            !               LIST(12) IS 1 OR 2.
            !
            !         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
            !               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
            !               QUANTITIES IN NUT IS:
            !
            !                        D PSI  (NUTATION IN LONGITUDE)
            !                        D EPSILON (NUTATION IN OBLIQUITY)
            !                        D PSI DOT
            !                        D EPSILON DOT
            !
            !           *   STATEMENT # FOR ERROR RETURN, IN CASE OF EPOCH OUT OF
            !               RANGE OR I/O ERRORS.
            !
            !     COMMON AREA STCOMX:
            !
            !          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
            !               STATES. KM = .TRUE., KM AND KM/SEC
            !                          = .FALSE., AU AND AU/DAY
            !               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
            !               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
            !
            !        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
            !               ONLY THE 9 PLANETS ARE AFFECTED.
            !                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
            !                             = .FALSE. =\ CENTER IS SUN
            !               DEFAULT VALUE = .FALSE.
            !
            !       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
            !               VELOCITY OF THE SUN.
            !
            !
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

            SAVE

            INTEGER OLDMAX
            PARAMETER ( OLDMAX = 400)
            INTEGER NMAX
            PARAMETER ( NMAX = 1000)

            DIMENSION ET2(2),PV(6,11),PNUT(4),T(2),PJD(4),BUF(1500),&
                SS(3),CVAL(NMAX),PVSUN(6)

            INTEGER LIST(12),IPT(3,13)

            LOGICAL FIRST
            DATA FIRST/.TRUE./

            CHARACTER*6 TTL(14,3),CNAM(NMAX)
            CHARACTER*80 NAMFIL

            LOGICAL KM,BARY

            COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
            COMMON/CHRHDR/CNAM,TTL
            COMMON/STCOMX/KM,BARY,PVSUN

!
!       ENTRY POINT - 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
!
            IF(FIRST) THEN
          
                FIRST=.FALSE.

                ! ************************************************************************
                ! ************************************************************************

                ! THE USER MUST SELECT ONE OF THE FOLLOWING BY DELETING THE 'C' IN COLUMN 1

                ! ************************************************************************

                !        CALL FSIZER1(NRECL,KSIZE,NRFILE,NAMFIL)
                !        CALL FSIZER2(NRECL,KSIZE,NRFILE,NAMFIL)
                CALL FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)

                IF(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'

                ! ************************************************************************
                ! ************************************************************************

                IRECSZ=NRECL*KSIZE
                NCOEFFS=KSIZE/2

                OPEN(NRFILE,FILE=NAMFIL,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=IRECSZ,STATUS='OLD')

                READ(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
                    ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3), &
                    (CNAM(L),L=OLDMAX+1,NCON)

                IF(NCON .LE. OLDMAX)THEN
                    READ(NRFILE,REC=2)(CVAL(I),I=1,OLDMAX)
                ELSE
                    READ(NRFILE,REC=2)(CVAL(I),I=1,NCON)
                ENDIF

                NRL=0

            ENDIF

            !       ********** MAIN ENTRY POINT **********

            IF(ET2(1) .EQ. 0.D0) RETURN

            S=ET2(1)-.5D0
            CALL SPLIT(S,PJD(1))
            CALL SPLIT(ET2(2),PJD(3))
            PJD(1)=PJD(1)+PJD(3)+.5D0
            PJD(2)=PJD(2)+PJD(4)
            CALL SPLIT(PJD(2),PJD(3))
            PJD(1)=PJD(1)+PJD(3)

            !       ERROR RETURN FOR EPOCH OUT OF RANGE

            IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98

            !       CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL

            NR=IDINT((PJD(1)-SS(1))/SS(3))+3
            IF(PJD(1).EQ.SS(2)) NR=NR-1

            tmp1 = DBLE(NR-3)*SS(3) + SS(1)
            tmp2 = PJD(1) - tmp1
            T(1) = (tmp2 + PJD(4))/SS(3)

            !       READ CORRECT RECORD IF NOT IN CORE

            IF(NR.NE.NRL) THEN
                NRL=NR
                READ(NRFILE,REC=NR,ERR=99)(BUF(K),K=1,NCOEFFS)
            ENDIF

            IF(KM) THEN
                T(2)=SS(3)*86400.D0
                AUFAC=1.D0
            ELSE
                T(2)=SS(3)
                AUFAC=1.D0/AU
            ENDIF

            !   INTERPOLATE SSBARY SUN

            CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)

            DO I=1,6
                PVSUN(I)=PVSUN(I)*AUFAC
            ENDDO

            !   CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED

            DO 4 I=1,10
                IF(LIST(I).EQ.0) GO TO 4

                CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I),&
                    LIST(I),PV(1,I))

                DO J=1,6
                    IF(I.LE.9 .AND. .NOT.BARY) THEN
                        PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
                    ELSE
                        PV(J,I)=PV(J,I)*AUFAC
                    ENDIF
                ENDDO

   4        CONTINUE

            !       DO NUTATIONS IF REQUESTED (AND IF ON FILE)

            IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0) THEN
                CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12),LIST(11),PNUT)
            ENDIF

            !       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)

            IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0) THEN
                CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13),LIST(12),PV(1,11))
            ENDIF
                
            RETURN

  98        WRITE(*,198)ET2(1)+ET2(2),SS(1),SS(2)
 198        FORMAT(' ***  Requested JED,',f12.2,' not within ephemeris limits,',2f12.2,'  ***')

            STOP

   99       WRITE(*,'(2F12.2,A80)')ET2,'ERROR RETURN IN STATE'

            STOP

        END SUBROUTINE

        SUBROUTINE CONST_EPH(NAM,VAL,SSS,N)

            !
            !     THIS ENTRY OBTAINS THE CONSTANTS FROM THE EPHEMERIS FILE
            !
            !     CALLING SEQEUNCE PARAMETERS (ALL OUTPUT):
            !
            !       NAM = CHARACTER*6 ARRAY OF CONSTANT NAMES
            !
            !       VAL = D.P. ARRAY OF VALUES OF CONSTANTS
            !
            !       SSS = D.P. JD START, JD STOP, STEP OF EPHEMERIS
            !
            !         N = INTEGER NUMBER OF ENTRIES IN 'NAM' AND 'VAL' ARRAYS
            !
            IMPLICIT DOUBLE PRECISION (A-H,O-Z)

            SAVE

            INTEGER NMAX
            PARAMETER (NMAX = 1000)

            CHARACTER*6 NAM(*),TTL(14,3),CNAM(NMAX)

            DOUBLE PRECISION VAL(*),SSS(3),SS(3),CVAL(NMAX),ZIPS(2)
            DOUBLE PRECISION PVST(6,11),PNUT(4)
            DATA ZIPS/2*0.d0/

            INTEGER IPT(3,13),DENUM,LIST(12)
            logical first
            data first/.true./

            COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
            COMMON/CHRHDR/CNAM,TTL

            !  CALL STATE TO INITIALIZE THE EPHEMERIS AND READ IN THE CONSTANTS

            IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT)
            first=.false.

            N=NCON

            DO I=1,3
                SSS(I)=SS(I)
            ENDDO

            DO I=1,N
                NAM(I)=CNAM(I)
                VAL(I)=CVAL(I)
            ENDDO

            RETURN

        END SUBROUTINE


end module