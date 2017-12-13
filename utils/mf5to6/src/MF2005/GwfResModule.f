      MODULE GwfResModule
        
        public
        
        INTEGER, SAVE,POINTER   ::NRES,IRESCB,NRESOP,IRESPT,NPTS
        INTEGER, SAVE, DIMENSION(:,:), POINTER ::IRES
        INTEGER, SAVE, DIMENSION(:,:), POINTER ::IRESL
        double precision,    SAVE, DIMENSION(:,:), POINTER ::BRES
        double precision,    SAVE, DIMENSION(:,:), POINTER ::CRES
        double precision,    SAVE, DIMENSION(:,:), POINTER ::BBRES
        double precision,    SAVE, DIMENSION(:),   POINTER ::HRES
        double precision,    SAVE, DIMENSION(:,:), POINTER ::HRESSE
      TYPE GWFRESTYPE
        INTEGER,POINTER   ::NRES,IRESCB,NRESOP,IRESPT,NPTS
        INTEGER,  DIMENSION(:,:), POINTER ::IRES
        INTEGER,  DIMENSION(:,:), POINTER ::IRESL
        double precision,     DIMENSION(:,:), POINTER ::BRES
        double precision,     DIMENSION(:,:), POINTER ::CRES
        double precision,     DIMENSION(:,:), POINTER ::BBRES
        double precision,     DIMENSION(:),   POINTER ::HRES
        double precision,     DIMENSION(:,:), POINTER ::HRESSE
      END TYPE
      TYPE(GWFRESTYPE), SAVE :: GWFRESDAT(10)
      
      contains
      
      SUBROUTINE GWF2RES7DA(IGRID)
C  Deallocate RES MEMORY
C
        CALL SGWF2RES7PNT(IGRID)
        DEALLOCATE(NRES)
        DEALLOCATE(IRESCB)
        DEALLOCATE(NRESOP)
        DEALLOCATE(IRESPT)
        DEALLOCATE(NPTS)
        DEALLOCATE(IRES)
        DEALLOCATE(IRESL)
        DEALLOCATE(BRES)
        DEALLOCATE(CRES)
        DEALLOCATE(BBRES)
        DEALLOCATE(HRES)
        DEALLOCATE(HRESSE)
C
      RETURN
      END SUBROUTINE GWF2RES7DA
      
      SUBROUTINE SGWF2RES7PNT(IGRID)
C  Change RES data to a different grid.
C
        NRES=>GWFRESDAT(IGRID)%NRES
        IRESCB=>GWFRESDAT(IGRID)%IRESCB
        NRESOP=>GWFRESDAT(IGRID)%NRESOP
        IRESPT=>GWFRESDAT(IGRID)%IRESPT
        NPTS=>GWFRESDAT(IGRID)%NPTS
        IRES=>GWFRESDAT(IGRID)%IRES
        IRESL=>GWFRESDAT(IGRID)%IRESL
        BRES=>GWFRESDAT(IGRID)%BRES
        CRES=>GWFRESDAT(IGRID)%CRES
        BBRES=>GWFRESDAT(IGRID)%BBRES
        HRES=>GWFRESDAT(IGRID)%HRES
        HRESSE=>GWFRESDAT(IGRID)%HRESSE
C
      RETURN
      END SUBROUTINE SGWF2RES7PNT
      
      SUBROUTINE SGWF2RES7PSV(IGRID)
C  Save RES data for a grid.
C
        GWFRESDAT(IGRID)%NRES=>NRES
        GWFRESDAT(IGRID)%IRESCB=>IRESCB
        GWFRESDAT(IGRID)%NRESOP=>NRESOP
        GWFRESDAT(IGRID)%IRESPT=>IRESPT
        GWFRESDAT(IGRID)%NPTS=>NPTS
        GWFRESDAT(IGRID)%IRES=>IRES
        GWFRESDAT(IGRID)%IRESL=>IRESL
        GWFRESDAT(IGRID)%BRES=>BRES
        GWFRESDAT(IGRID)%CRES=>CRES
        GWFRESDAT(IGRID)%BBRES=>BBRES
        GWFRESDAT(IGRID)%HRES=>HRES
        GWFRESDAT(IGRID)%HRESSE=>HRESSE
C
      RETURN
      END SUBROUTINE SGWF2RES7PSV
      
      END MODULE GwfResModule

