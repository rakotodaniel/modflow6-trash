      module GwfResSubs
        
        use GwfResModule, only: SGWF2RES7PNT, SGWF2RES7PSV
        
        private
        public :: GWF2RES7AR, GWF2RES7RP
        
      contains

      SUBROUTINE GWF2RES7AR(IN,IGRID)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR RESERVOIRS, AND READ RESERVOIR
C     LOCATIONS, LAYER, CONDUCTANCE, BOTTOM ELEVATION, AND BED THICKNESS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,DELR,DELC,
     1                       PERLEN,NSTP,TSMULT
      USE GWFBASMODULE, ONLY:DELT
      USE GWFRESMODULE, ONLY:NRES,IRESCB,NRESOP,IRESPT,NPTS,
     1                       IRES,IRESL,BRES,CRES,BBRES,HRES,HRESSE
      use SimModule, only: ustop
!      use InputOutputModule, only: U2DINT, U2DREL
      use ArrayReadersMF5Module, only: ReadArray
C
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'      RESERVOIR LOCATION'/
      DATA ANAME(2) /'   RESERVOIR LAYER INDEX'/
      DATA ANAME(3) /'RESERVOIR LAND SURF ELEV'/
      DATA ANAME(4) /'  RES. BED VERT HYD COND'/
      DATA ANAME(5) /' RESERVOIR BED THICKNESS'/
C     ------------------------------------------------------------------
      ALLOCATE (NRES,IRESCB,NRESOP,IRESPT,NPTS)
C
C1------IDENTIFY PACKAGE AND INITIALIZE
      WRITE(IOUT,1)IN
    1 FORMAT(/,'RES7 -- RESERVOIR PACKAGE, VERSION 7, 2/28/2006',
     1' INPUT READ FROM UNIT',I3)
C
C2------READ & PRINT NUMBER OF RESERVOIRS AND FLAGS FOR
C2------RESERVOIR OPTIONS
      READ(IN,2) NRES,IRESCB,NRESOP,IRESPT,NPTS
    2 FORMAT(5I10)
C
C2A-----CHECK TO SEE THAT NUMBER OF RESERVOIRS IS AT LEAST 1,
C2A-----PRINT VALUE
      IF(NRES.GT.0) THEN
       WRITE(IOUT,6) NRES
    6 FORMAT(1X,'TOTAL NUMBER OF RESERVOIRS: ',I3)
      ELSE
       WRITE (IOUT,7)
    7 FORMAT(1X,'ABORTING, NUMBER OF RESERVOIRS LESS THAN 1...')
      CALL USTOP(' ')
      ENDIF 
C
C2B-----CHECK FLAG FOR CELL-BY-CELL OUTPUT, PRINT VALUE
      IF(IRESCB.GT.0) WRITE(IOUT,10) IRESCB
 10   FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE RECORDED ON UNIT',I3)
C2C-----CHECK TO SEE THAT RESERVOIR LAYER OPTION FLAG IS LEGAL,
C2C-----PRINT VALUE
      IF(NRESOP.GE.1.AND.NRESOP.LE.3)GO TO 200
C
C2C1----IF ILLEGAL PRINT A MESSAGE AND ABORT SIMULATION
      WRITE(IOUT,8)
    8 FORMAT(1X,'ILLEGAL OPTION CODE. SIMULATION ABORTING')
      CALL USTOP(' ')
C
C2C2----IF OPTION IS LEGAL PRINT OPTION CODE.
 200  CONTINUE
      IF(NRESOP.EQ.1) WRITE(IOUT,201)
  201 FORMAT(1X,'OPTION 1 -- RESERVOIR CONNECTED TO TOP LAYER')
      IF(NRESOP.EQ.2) WRITE(IOUT,202)
  202 FORMAT(1X,'OPTION 2 -- RESERVOIR CONNECTED TO ONE SPECIFIED',
     1        ' NODE IN EACH VERTICAL COLUMN')
      IF(NRESOP.EQ.3) WRITE(IOUT,203)
  203 FORMAT(1X,'OPTION 3 -- RESERVOIR CONNECTED TO HIGHEST',
     1        ' ACTIVE NODE IN EACH VERTICAL COLUMN')
C
C2D-----PRINT VALUE FOR RESERVIOR PRINT OPTION FLAG
      IF(IRESPT.GT.0) WRITE(IOUT,14) 
 14   FORMAT(1X,'RESERVOIR HEADS, AREAS, AND VOLUMES ',
     1 'WILL BE PRINTED EACH TIME STEP')
C2E-----PRINT NUMBER OF POINTS TO BE USED IN CALCULATING TABLE
C2E-----OF RESERVOIR STAGE VS. AREA AND VOLUME
      IF(NPTS.LT.1) THEN
       WRITE(IOUT,*) ' Table of reservoir areas and volumes ',
     1 'will not be calculated.'
      ELSE
       WRITE(IOUT,9) NPTS
 9     FORMAT(I5,' points will be used in constructing table of ',
     1  'reservoir areas and volumes.')
      ENDIF
C
C3------ALLOCATE SPACE FOR ARRAYS.
      ALLOCATE (IRES(NCOL,NROW))
      ALLOCATE (IRESL(NCOL,NROW))
      ALLOCATE (BRES(NCOL,NROW))
      ALLOCATE (CRES(NCOL,NROW))
      ALLOCATE (BBRES(NCOL,NROW))
      ALLOCATE (HRES(NRES))
      ALLOCATE (HRESSE(2,NRES))
C
C4------READ INDICATOR ARRAY SHOWING LOCATIONS OF RESERVOIRS (IRES)
      KK=1
      CALL ReadArray(IRES,ANAME(1),NROW,NCOL,KK,IN,IOUT)
C5------VERIFY LOCATIONS EXIST FOR ALL RESERVOIRS
      DO 36 N=1,NRES
      NCELL=0
      DO 30 I=1,NROW
      DO 20 J=1,NCOL
      IF(IBOUND(J,I,1).LE.0) IRES(J,I)=0
      IF(IRES(J,I).EQ.N) NCELL=NCELL+1
   20 CONTINUE
   30 CONTINUE
      IF(NCELL.GT.0) THEN
       WRITE(IOUT,32) N,NCELL
   32  FORMAT(1X,'NUMBER OF CELLS IN RESERVOIR ',I2,':',I6)
      ELSE
       WRITE(IOUT,34)
   34 FORMAT(1X,'NO ACTIVE CELLS FOUND FOR RESERVOIR ',I2,'.',
     1 '  ABORTING...')
      ENDIF
   36 CONTINUE
C
C6------IF NRESOP=2 THEN A LAYER INDICATOR ARRAY IS NEEDED.
      IF (NRESOP.NE.2)GO TO 37
      CALL ReadArray(IRESL,ANAME(2),NROW,NCOL,0,IN,IOUT)
C7------READ IN BOTTOM ELEVATION, BED CONDUCTIVITY, AND BED THICKNESS
   37 CALL ReadArray(BRES,ANAME(3),NROW,NCOL,KK,IN,IOUT)
      CALL ReadArray(CRES,ANAME(4),NROW,NCOL,KK,IN,IOUT)
      CALL ReadArray(BBRES,ANAME(5),NROW,NCOL,KK,IN,IOUT)
C8------CONVERT RESERVOIR BED HYDRAULIC CONDUCTIVITY TO CONDUCTANCE
C8------BED THICKNESS TO ELEVATION OF BOTTOM OF RESERVOIR BED  
      DO 40 I=1,NROW
      DO 38 J=1,NCOL
      IF(IRES(J,I).LE.0) GO TO 38
      IF(IRES(J,I).GT.NRES) GO TO 38
       CRES(J,I)=CRES(J,I)*DELC(I)*DELR(J)/BBRES(J,I)
       BBRES(J,I)=BRES(J,I)-BBRES(J,I)
   38 CONTINUE
   40 CONTINUE
C9------MAKE STAGE-VOLUME TABLE FOR EACH RESERVOIR
      DO 60 N=1,NRES
C9A-----FIND MAX AND MIN BOTTOM ELEVATION
      ELMIN=9.99E10
      ELMAX=-9.99E10
      DO 44 I=1,NROW
      DO 42 J=1,NCOL
      IF(IRES(J,I).NE.N) GO TO 42
      IF(BRES(J,I).LT.ELMIN) ELMIN=BRES(J,I)
      IF(BRES(J,I).GT.ELMAX) ELMAX=BRES(J,I)
   42 CONTINUE
   44 CONTINUE
C9B-----CONSTRUCT TABLE
      WRITE(IOUT,46) N,ELMIN
   46 FORMAT(1X,'STAGE-VOLUME TABLE FOR RESERVOIR',I2,/,6X,
     1 'STAGE       VOLUME         AREA',/,
     2 3X,36('-'),/,1X,G12.5,2(10X,'0.0'))
      IF(NPTS.LT.1) GO TO 60
      DEL=(ELMAX-ELMIN)/FLOAT(NPTS)
      STAGE=ELMIN
      DO 56 NP=1,NPTS
      STAGE=STAGE+DEL
      VOL=0.0
      TAREA=0.0
      DO 50 I=1,NROW
      DO 48 J=1,NCOL
      IF(IRES(J,I).NE.N) GO TO 48
      IF(STAGE.GT.BRES(J,I))THEN
       AREA=DELR(J)*DELC(I)
       TAREA=TAREA+AREA
       VOL=VOL+AREA*(STAGE-BRES(J,I))
      ENDIF
   48 CONTINUE
   50 CONTINUE
      WRITE(IOUT,54) STAGE,VOL,TAREA
   54 FORMAT(1X,G12.5,2G13.5)
   56 CONTINUE
      WRITE(IOUT,58)
   58 FORMAT(1X,' ')
   60 CONTINUE
C
C10-----RETURN
      CALL SGWF2RES7PSV(IGRID)
      RETURN
      END SUBROUTINE GWF2RES7AR
      
      SUBROUTINE GWF2RES7RP(IN,IGRID)
C     ******************************************************************
C     READ START AND END HEADS FOR EACH RESERVOIR FOR CURRENT
C     STRESS PERIOD
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFRESMODULE, ONLY:NRES,HRESSE
C     ------------------------------------------------------------------
      CALL SGWF2RES7PNT(IGRID)
C
      DO 80 N=1,NRES
      READ(IN,64) HRESSE(1,N),HRESSE(2,N)
   64 FORMAT(2F10.0)
   80 CONTINUE
C
C8------RETURN
      RETURN
      END SUBROUTINE GWF2RES7RP
      
      end module GwfResSubs
      
