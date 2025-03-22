      SUBROUTINE get_gas_opacity_corrk_wrapper
          include 'rcommons.h'
          WRITE(*,*) "In get_gas_opacity_corrk_wrapper"
          call get_gas_opacity_corrk(METALLICITY,0.0, FBASEFLUX, GASCON, with_TiO_and_VO)
      END SUBROUTINE get_gas_opacity_corrk_wrapper
      
      SUBROUTINE get_gas_opacity_corrk(METALLICITY,C_TO_O, FBASEFLUX, GASCON, with_TiO_and_VO)
          implicit none
          COMMON/CORRKGAS/OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC,
     &                    INT_SPEC, OPAC_CIA
          COMMON/PLANCK_INT/PLANCK_INTS, PLANCK_TS
          INTEGER :: I, J, K, L
          REAL, INTENT(IN) :: METALLICITY, C_TO_O, FBASEFLUX
          REAL :: TINT
          INTEGER :: TINT_INDEX
          REAL :: GASCON, MMW
          REAL :: PLANCK_INTS(11, 3925), PLANCK_TS(3925)
          REAL :: with_TiO_and_VO

          CHARACTER(LEN=32) :: dummyMETALLICITY_str, dummyC_TO_O_str, tiovo_str

          CHARACTER(len=80) :: FNM, CIAFNM
          REAL :: TEMP_OPAC_CORRK(8,73*20*11) ! 8 gauss pts, 73 T, 20 P, 11 wavelengths
          REAL :: TEMP_OPAC_CIA(11,73*20) ! 73 T, 20 P, 11 wavelengths
          ! COMMON BLOCK SHAPES HERE:
          REAL :: OPAC_CORRK(73, 20, 11, 8)
          REAL :: TS_CORRK(73), PS_CORRK(20), TS_LOG_CORRK(73), PS_LOG_CORRK(20), WGTS_CORRK(8) 
          REAL :: WNO_EDGES(12), WNO_CTRS(11), STEL_SPEC(11), INT_SPEC(11)
          REAL :: OPAC_CIA(73, 20, 11)

          
          MMW = 8314.462618 / GASCON ! Mean molecular weight in g/mol (or H masses per molecule)

          ! Convert METALLICITY and C_TO_O to strings
          if (METALLICITY .eq. 0.0) then
            dummyMETALLICITY_str = '+000'
          else
            write(*,*) "Thomas hasn't coded non-solar metallicity k-tables yet"
            stop
          end if
          if (C_TO_O .eq. 0.0) then
            dummyC_TO_O_str = '100'
          else
            write(*,*) "Thomas hasn't coded non-solar C/O k-tables yet"
            stop
          end if
          ! write(*,*) 'with_TiO_and_VO: ', with_TiO_and_VO
          if (with_TiO_and_VO .eq. 1.) then
            tiovo_str = '_witiovo'
          else if (with_TiO_and_VO .eq. 2.) then
            tiovo_str = '_notiovo'
          else
            write(*,*) "tiovo flag not recognized, setting with_TiO_and_VO = 1 for the first loop"
            tiovo_str = '_witiovo'
          end if
          ! write(*,*) 'METALLICITY: ', dummyMETALLICITY_str
          ! write(*,*) 'C_TO_O: ', dummyC_TO_O_str

          FNM= "../kcoeffs/feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)//trim(tiovo_str)//".txt"
          CIAFNM= "../CIA/cia_feh"//trim(dummyMETALLICITY_str)//"_co_"//trim(dummyC_TO_O_str)//trim(tiovo_str)//".txt"
          ! Read in k-coefficients from file
          write(*,*) 'Reading in k-coefficients from file: ', FNM
          write(*,*) 'Reading in CIA coefficients from file: ', CIAFNM
          OPEN(UNIT=1, FILE=FNM)
          READ(1,*) TEMP_OPAC_CORRK
          CLOSE(1)
          OPEN(UNIT=1, FILE=CIAFNM)
          READ(1,*) TEMP_OPAC_CIA
          CLOSE(1)

          ! Reshape into 4D array
          ! Index order of new array: T, P, wavelength, gauss pt
          ! write(*,*) 'test gauss: ', TEMP_OPAC_CORRK(1,2)
          DO I = 1, 73
            DO J = 1, 20
              DO K = 1, 11
                DO L = 1, 8
                  ! After exponential, we're in cm^2/molecule units which we convert to m^2/molecule since RT uses SI units
                  ! avagadro + MMW/1000 conversion factor get us to m^2/kg
                  OPAC_CORRK(I,J,K,L) = LOG10(EXP(TEMP_OPAC_CORRK(L,(I-1)*20*11+(J-1)*11+K)) * 1e-4 * 6.022e23 / (MMW/1000))
                END DO
              END DO
            END DO
          END DO

          DO I = 1, 73
            DO J = 1, 20
              DO K = 1, 11
                ! after this line we have log10(m^2/kg)
                ! file is in units of m^2/molecule
                OPAC_CIA(I,J,K) = LOG10(TEMP_OPAC_CIA(K,(I-1)*20+J) * 6.022e23 / (MMW/1000))
              END DO
            END DO
          END DO
          write(*,*) 'OPAC_CIA: ', OPAC_CIA
          write(*,*) 'loc(OPAC_CIA:) ', loc(OPAC_CIA)
          write(*,*) loc(OPAC_CORRK)
          write(*,*) ''
          write(*,*) ''
          write(*,*) ''
          write(*,*) ''
          ! write(*,*) 'test gauss: ', OPAC_CORRK(35,16,7,8)

          ! Set up the gauss weights of the 8 gauss points
          WGTS_CORRK = (/0.16523105, 0.30976895, 0.30976895, 0.16523105, 0.00869637,
     &              0.01630363, 0.01630363, 0.00869637/)

          ! Read in the temperature and pressure points of our grid
          OPEN(UNIT=1, FILE="../kcoeffs/Tpts.txt")
          READ(1,*) TS_CORRK
          CLOSE(1)

          OPEN(UNIT=1, FILE="../kcoeffs/Ppts.txt")
          READ(1,*) PS_CORRK
          CLOSE(1)
          PS_CORRK = PS_CORRK * 1.0E5 ! Convert bars to Pascals, which the RT uses

          ! Invert T and log P for interpolation
          TS_LOG_CORRK = LOG10(TS_CORRK)
          PS_LOG_CORRK = LOG10(PS_CORRK)

          ! Read in the wavenumber points of our grid
          OPEN(UNIT=1, FILE="../kcoeffs/wno_edges.txt")
          READ(1,*) WNO_EDGES
          CLOSE(1)

          OPEN(UNIT=1, FILE="../kcoeffs/wno_ctrs.txt")
          READ(1,*) WNO_CTRS
          CLOSE(1)

          ! Read in the stellar spectrum
          CALL BIN_STELLAR_SPECTRUM("../kcoeffs/stellar_spectrum.txt", STEL_SPEC, WNO_EDGES)
          ! WRITE(*,*) 'STELLAR SPECTRUM: ', STEL_SPEC
          ! WRITE(*,*) 'total stellar flux: ', SUM(STEL_SPEC)

          ! Grab Planck function integrals for each tmeperature and bin:
          OPEN(UNIT=1, FILE="../kcoeffs/planck_integrals.txt")
          READ(1,*) PLANCK_INTS
          CLOSE(1)
          OPEN(UNIT=1, FILE="../kcoeffs/planck_temps.txt")
          READ(1,*) PLANCK_TS
          CLOSE(1)
          ! window, temp
          ! write(*,*) 'Planck integrals: ', PLANCK_INTS(1,1), PLANCK_INTS(1,2), PLANCK_INTS(2,1)
          ! write(*,*) 'Planck temps: ', PLANCK_TS(1), PLANCK_TS(2), PLANCK_TS(3)

          ! Now split up internal flux by channel:
          TINT = (FBASEFLUX / 5.670367E-8) ** 0.25 ! Convert flux to temperature
          TINT_INDEX = MINLOC(ABS(PLANCK_TS - TINT), 1)
          INT_SPEC = PLANCK_INTS(:, TINT_INDEX)
          INT_SPEC = INT_SPEC / SUM(INT_SPEC)
          write(*,*) 'TINT: ', TINT, 'TINT_INDEX: ', TINT_INDEX
          write(*,*) 'INT_SPEC: ', INT_SPEC
      END SUBROUTINE get_gas_opacity_corrk

      SUBROUTINE BIN_STELLAR_SPECTRUM(file_name, STEL_SPEC, WNO_EDGES)
          implicit none
          CHARACTER(len=80) :: file_name
          REAL :: STEL_SPEC(11), STEL_SPEC_IN(2, 1000)
          REAL :: WNO_EDGES(12)
          INTEGER :: I, J
          REAL :: WNO(1000), FLUX(1000)
          ! Read in the stellar spectrum
          OPEN(UNIT=1, FILE=file_name)
          READ(1,*) STEL_SPEC_IN
          CLOSE(1)
          DO I = 1, 1000
            WNO(I) = STEL_SPEC_IN(1,I)
            FLUX(I) = STEL_SPEC_IN(2,I)
          END DO
          ! Trapezoidal integration in each spectral band
          DO J = 1, 11
            STEL_SPEC(J) = 0.0
            DO I = 1, 999
              IF (WNO(I) >= WNO_EDGES(J) .AND. WNO(I) < WNO_EDGES(J+1)) THEN
                STEL_SPEC(J) = STEL_SPEC(J) + ((FLUX(I) + FLUX(I+1))/2 *((WNO(I+1) - WNO(I)) * 3e10))
              END IF
            END DO
          END DO

          ! Normalize the spectrum so the fort.7's SOLC_IN will define the total flux:
          STEL_SPEC = STEL_SPEC / SUM(STEL_SPEC)

      END SUBROUTINE BIN_STELLAR_SPECTRUM