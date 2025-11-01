
      PROGRAM RTMAIN
      call RT()
      END
      SUBROUTINE RT()
        ! IMPLICIT NONE
      include 'params.i'
      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL
     +,NWW=1+(MM-1)/MOCT)
      PARAMETER(IGA=NWJ2*NHEM,IGB=IDB*NHEM,IGC=MGPP*NHEM,IGD=IDD*NHEM
     +,IGG=IDG*NHEM,IGL=IDL*NHEM,IGM=IDM*NHEM,IGN=IDN*NHEM
     +,IGO=IGA+IGA,IGP=IGB+IGB,NFTWG=(5+NTRAC)*NL+3
     +,NFTGW=(6+3*NTRAC)*NL+2,NFTGD=(3+NTRAC)*NL,NLTR=NL*NTRAC)
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)
     +              ,BEGDOY,DOY
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)
     +              ,TKP(NL),C(NL2),SQH(NNP)
     +              ,MF,MFP,JZF,NF
     +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT
     +              ,RD,RV,CPD,CLATNT
     +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH
     +              ,LFLUX
     +              ,LBALAN,LRESTIJ
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT
     +              ,LNOISE,NFP
      COMPLEX EZ,AIOCT
    !   LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ
    !  +       ,LFLUX,LNOISE
    !  +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT
C
C
C     Legendre polynomials and information about gaussian latitudes
C
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC
      COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN, TAULIMIT
      LOGICAL LPLOTMAP
      COMMON/RESTIJ/TTRES(IGB)
     + ,DTNS,DTEP,DTTRP,FAC(NL),DDAMP(NL),TFRC(NL),YRLEN,TRS(NL)
     +  ,RESTTT(NL),REDTEP(NL)
     +  ,ALR,ZTROP,TGR
      COMPLEX TTRES
      COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM(3), AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3)
    !   LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
    !  + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS
      ! CHARACTER(30) :: AEROSOLMODEL
      ! CHARACTER(30) :: AEROSOLCOMP
      ! REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(nl+1),TCON(nl+1)
      ! REAL MOLEF(13)
      ! LOGICAL DELTASCALE,HAZES
      COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,MOLEF,AERLAYERS
      NAMELIST/COMMENT/THECOMMENT
      CHARACTER(70) :: THECOMMENT
      NAMELIST/INPPL/GA,GASCON,RADEA,AKAP,WW,P0,RV,CLATNT
      REAL :: GA,GASCON,RADEA,AKAP,WW,P0,RV,CLATNT
      NAMELIST/INPRN/KRUN,BEGDAY,TSPD,KITS,PNU,TDISS,LFLUX,BEGDOY,NDEL,T0,LRSTRT,LSTRETCH,LSHORT,LTVEC,LBALAN,LRESTIJ,
     &  LCLIM,LPERPET,L22L,LOROG,LCSFCT,KOLOUR,LNOISE,LMASCOR,LMASOLD,LMASPRT
      REAL :: BEGDAY, TSPD, PNU, TDISS, BEGDOY, T0
      INTEGER :: KRUN, KITS, NDEL
      LOGICAL :: LFLUX, LRSTRT, LSTRETCH, LSHORT, LTVEC, LBALAN, LRESTIJ
      LOGICAL :: LCLIM, LPERPET, L22L, LOROG, LCSFCT
      LOGICAL :: LNOISE, LMASCOR, LMASOLD, LMASPRT
      ! INTEGER :: KOLOUR(NTRAC)
      
      NAMELIST/INVARPARAM/OOM_IN, NTSTEP_IN, FBASEFLUX, OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, PORB, OBLIQ, ECCEN
      REAL :: OOM_IN, FBASEFLUX, OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, PORB, OBLIQ, ECCEN
      NAMELIST/INPRSIJ/RESTIM,RESTTT,TGR ! newtonian forcing timescale and initial condition as functions of layer and approx deep temp
      REAL :: RESTIM(NL), RESTTT, TGR

      NAMELIST/INSIMPRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,JSKIPLON,JSKIPLAT,DOSWRAD,DOLWRAD,LWSCAT,FLXLIMDIF,SURFEMIS,RAYSCAT,
     & RAYSCATLAM,AEROSOLS,ABSSW,ABSLW,ALBSW,NEWTB,NEWTE,with_TiO_and_VO,picket_fence_optical_depths
      LOGICAL :: LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD,LWSCAT,FLXLIMDIF,RAYSCAT, AEROSOLS, picket_fence_optical_depths
      REAL :: RAYSCATLAM,ABSSW,ABSLW,ALBSW,NEWTB,NEWTE
      INTEGER :: WITH_TIO_AND_VO

      NAMELIST/INCLOUDY/AEROSOLMODEL,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,AERTOTTAU,
     & CLOUDBASE,CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,ASYMLW,DELTASCALE,SIG_AREA,PHI_LON
      CHARACTER(30) :: AEROSOLMODEL, AEROSOLCOMP
      REAL :: MTLX, METALLICITY, MOLEF(13)
      INTEGER :: AERLAYERS
      LOGICAL HAZES, PICKET_FENCE_CLOUDS, DELTASCALE
      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),MAXTAU,TCON(NL+1)
      NAMELIST/ONED/MUSTEL,ISF
      REAL :: MUSTEL, ISF
      REAL :: fracchange, maxfracchange
      integer :: TSTEPINTERVAL, TSTEPINTERVALRECORD, NEXTRESTART
      CHARACTER(len=20) :: myCharacterString


    !   COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)
    !  :              ,TTVD(IGC,NL),QTVD(IGC,NL),TG(IGC,NL)
    !  :              ,TRAG(IGC,NL,NTRAC)
    !  :              ,PLG(IGC),TYBL(IGC),TXBL(IGC)
    !  :              ,SPG(IGC),VPG(IGC),TTRD(IGC,NL)
    !  :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),UNLG(IGC,NL)
    !  :              ,VNLG(IGC,NL),TTLR(IGC,NL),UTRAG(IGC,NL,NTRAC)
    !  :              ,TTCR(IGC,NL),VTRAG(IGC,NL,NTRAC)
    !  :              ,UTVD(IGC,NL),VTVD(IGC,NL)
    !  :         ,ASSBL(IGC),ASHBL(IGC),ASLBL(IGC),ARRCR(IGC),ARRLR(IGC)
    !  :         ,arflux(igc,6),asfld(igc,6),acld(igc,4)
    !  :         ,SSBL(IGC),SHBL(IGC),SLBL(IGC),RRCR(IGC),RRLR(IGC)
    !  :         ,rflux(igc,6),sfld(igc,6),cld(igc,4)
    !   REAL :: TTRD
      COMMON/ONEDRT/TTKPD(NL)
      REAL TTKPD

      ! Here we declare the variables we care the most about
      REAL :: T(NL+1) ! temp at layer centers + bottom edge (K)
      REAL :: PR(NL+1) ! pressure at layer centers + bottom edge (Pa)
      INTEGER :: I
      INTEGER :: TSTEPSTOT
      CALL INISET
      CALL INITAL


      ! ! Let's get what we need from the fort.7 first
      REWIND(7)
      READ(7,INPRN)
      WRITE(2,INPRN)
      REWIND(7)
      CALL INISET


      ITSPD = TSPD
      EZ=1.0/SQRT(.375)                                                   
      CV=RADEA*WW                                                         
      CG=CV*CV                                                            
      CT=CG/GASCON                                                        
      CQ=1000.0                                                           
      PFAC=0.5*CV*CV*1.0E5/GA                                             
      SQR2=SQRT(2.0)                                                      
      RSQR2=1.0/SQR2                                                      
      EAM1=SQR2/3.                                                        
      EAM2=SQRT(2./45.) 

      ! CALL get_cloud_scattering_properties_wrapper
      ! IC setup:
      ! just took this from a 7 OOM 50 lay run with the regular GCM, should eventually make this more general though
      write(*,*) 'SIGMA:', SIGMA
      
      DO I = 1, NL
        PR(I) = SIGMA(I)*P0
        T(I) = RESTTT(I)
      END DO
      PR(NL+1) = P0
      T(NL+1) = TGR
      ! T(:) = RESTTT
      READ(7,ONED)
      WRITE(2,ONED)
      write(*,*) 'MUSTEL, ISF: ', MUSTEL, ISF
      CALL RADIATION(0, 1, T, PR, MUSTEL, ISF)
      write(*,*) ""
      WRITE(*,*) TTKPD ! is correct heating rate output in units of K/day (I have assumed real day and not planet day)
      ! now our time-stepping loop, the general idea is to call RADIATION at each time step and get heating rates
      ! then add those on to the T array appropriately, maybe with RK4 tstepping?
      TSTEPSTOT = 0
      KOUNTP = 100 * TSPD
      write(*,*) 'KRUN=', KRUN
      TSTEPINTERVAL = 1
      TSTEPINTERVALRECORD = TSTEPINTERVAL
      NEXTRESTART = 100 * TSPD
      DO WHILE (TSTEPSTOT < KRUN)
        TTKPD = 0.0
        CALL RADIATION(0, 1, T, PR, MUSTEL, ISF)
        maxfracchange = ABS(TTKPD(1)) / TSPD / 86400.0 * PI2 / WW * TSTEPINTERVAL / T(1)
        DO I = 1, NL
          fracchange = ABS(TTKPD(I)) / TSPD / 86400.0 * PI2 / WW * TSTEPINTERVAL / T(I)
          IF (fracchange .GT. maxfracchange) THEN
            maxfracchange = fracchange
          END IF
          T(I) = T(I) + TTKPD(I)/86400.0 / (TSPD) * PI2 / WW * TSTEPINTERVAL ! converting from K/day to K/s and multiplying by time step
        END DO

        TSTEPSTOT = TSTEPSTOT + TSTEPINTERVAL
        IF (TSTEPSTOT .EQ. NEXTRESTART) THEN
          WRITE(*,*) 'After time step ', TSTEPSTOT
          WRITE(*,*) 'day', TSTEPSTOT/TSPD, '/', KRUN/TSPD
          WRITE(*,*) T(1)
          OPEN(UNIT=8, FILE='fort.8', STATUS='REPLACE')
          WRITE(8,*) 'Day ', TSTEPSTOT/TSPD, ': P(Pa) T(K) dTdt(K/day)'
          DO I=1, NL
            WRITE(8,*) PR(I), T(I), TTKPD(I)
          END DO
          CLOSE(8)
          write(myCharacterString, '(F8.0)') TSTEPSTOT/TSPD
          CALL RENAME('fort.8', TRIM(ADJUSTL(myCharacterString))//'day_fort.8')

          TSTEPINTERVAL = TSTEPINTERVALRECORD
          NEXTRESTART = NEXTRESTART + KOUNTP
        END IF
        
        
        if (maxfracchange .gt. 0.01) then
          TSTEPINTERVAL = MAX(1, (TSTEPINTERVAL / 2))
        else if (maxfracchange .lt. 0.0005) then
          ! formally, we'd really like to cap this at the highest power of 2 smaller than KOUNTP
          TSTEPINTERVAL = MIN(TSTEPINTERVAL * 2, 64)
        end if

        if ((TSTEPSTOT + TSTEPINTERVAL) > NEXTRESTART) then
          write(*,*) 'Adjusting next time step interval from ', TSTEPINTERVAL, ' to ', NEXTRESTART - TSTEPSTOT
          TSTEPINTERVALRECORD = TSTEPINTERVAL
          TSTEPINTERVAL = NEXTRESTART - TSTEPSTOT
        end if

        if (ISNAN(T(1))) then
          write(*,*) 'Temperature blew up to NaN at time step ', TSTEPSTOT
          OPEN(UNIT=8, FILE='fort.8', STATUS='REPLACE')
          WRITE(8,*) 'Day ', TSTEPSTOT/TSPD, ': P(Pa) T(K) dTdt(K/day)'
          DO I=1, NL
            WRITE(8,*) PR(I), T(I), TTKPD(I)
          END DO
          CLOSE(8)
          write(myCharacterString, '(F8.0)') KRUN/TSPD
          CALL RENAME('fort.8', TRIM(ADJUSTL(myCharacterString))//'day_fort.8')
          EXIT
        end if
      END DO
      WRITE(*,*) 'Final temperatures after ', TSTEPSTOT, ' time steps:'
      WRITE(*,*) T
      ! now write to a file called fort.8
      
      END SUBROUTINE RT
