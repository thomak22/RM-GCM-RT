      SUBROUTINE TWOSTR(TAUL,solar_calculation_indexer,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
     &  SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &  GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &  WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &  uG0,uTAUL,W0,uW0,uopd,U1S,
     &  U1I,TOON_AK,B1,B2,EE1,EM1,
     &  EM2,EL1,EL2,GAMI,AF,
     &  BF,EF,SFCS,B3,CK1,CK2,
     &  CP,CPB,CM,CMB,DIRECT,EE3,
     &  EL3,FNET,TMI,AS,DF,
     &  DS,XK,DIREC,DIRECTU,DINTENT,
     &  UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &  tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, num_layers)
!
!    ******************************************************************
!    *  Purpose             :  Defines matrix properties and sets up  *
!    *                         matrix coefficients that do not depend *
!    *                         on zenith angle or temperature.        *
!    *  Subroutines Called  :  None                                   *
!    *  Input               :  W0, G0                                 *
!    * ****************************************************************
!
      use corrkmodule, only : MINWNOSTEL
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NTOTAL,2*NL+2), WOL(NTOTAL,2*NL+2), WAVE(NTOTAL+1), TT(NL+1), Y3(NTOTAL,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NTOTAL), PTEMPT(NTOTAL), G0(NTOTAL,2*NL+2), OPD(NTOTAL,2*NL+2), PTEMP(NTOTAL,2*NL+2)
      REAL uG0(NTOTAL,2*NL+2), uTAUL(NTOTAL,2*NL+2), W0(NTOTAL,2*NL+2), uW0(NTOTAL,2*NL+2), uopd(NTOTAL,2*NL+2),  U1S(NTOTAL)
      REAL U1I(NTOTAL), TOON_AK(NTOTAL,2*NL+2), B1(NTOTAL,2*NL+2), B2(  NTOTAL,2*NL+2), EE1(NTOTAL,2*NL+2), EM1(NTOTAL,2*NL+2)
      REAL EM2(NTOTAL,2*NL+2), EL1(NTOTAL,2*NL+2), EL2(NTOTAL,2*NL+2), GAMI(NTOTAL,2*NL+2), AF(NTOTAL,4*NL+4)
      REAL BF(NTOTAL,4*NL+4), EF(NTOTAL,4*NL+4), SFCS(NTOTAL), B3(NTOTAL,2*NL+2), CK1(NTOTAL,2*NL+2), CK2(NTOTAL,2*NL+2)
      REAL CP(NTOTAL,2*NL+2), CPB(NTOTAL,2*NL+2), CM(NTOTAL,2*NL+2), CMB(NTOTAL,2*NL+2), DIRECT(NTOTAL,2*NL+2), EE3(NTOTAL,2*NL+2)
      REAL EL3(NTOTAL,2*NL+2), FNET(NTOTAL,2*NL+2), TMI(NTOTAL,2*NL+2), AS(NTOTAL,4*NL+4), DF(NTOTAL,4*NL+4)
      REAL DS(NTOTAL,4*NL+4), XK(NTOTAL,4*NL+4), DIREC(NTOTAL,2*NL+2), DIRECTU(NTOTAL,2*NL+2), DINTENT(NTOTAL,3,2*NL+2)
      REAL UINTENT(NTOTAL,3,2*NL+2), TMID(NTOTAL,2*NL+2), TMIU(NTOTAL,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      real, dimension(NTOTAL,2*NL+2) :: TAUL
      integer solar_calculation_indexer

      U1I(:) = 0.
      EF(:,:) = 0.0
      !  WRITE(*,*), 'solar_calculation_indexer', solar_calculation_indexer, 'LLA', LLA
       DO 10 L    =  MAX(solar_calculation_indexer,MINWNOSTEL*8),LLA
          if( L .LE. NSOL )then
            U1I(L) = SQ3
          else
            U1I(L) = 2.0
          endif
  10      U1S(L)  =  TPI/U1I(L) ! TPI = 2pi


!      HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
!      OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
!      NEEDED FOR MATRIX.


!
       DO 14 J          =  1,NLAYER
          DO 14 L       =  MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)    =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)    =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             TOON_AK(L,J)    = SQRT(ABS(B1(L,J)*B1(L,J) - B2(L,J)*B2(L,J)))
             GAMI(L,J)  =  B2(L,J)/(B1(L,J) + TOON_AK(L,J))
             EE1(L,J)   =  EXP(-TOON_AK(L,J)*TAUL(L,J))
             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)  !e1
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J) !e2
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)       !e3
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)       !e4
  14  CONTINUE

      DO 15 J = 1,NDBL
          DO 15 L = NSOL+1,NTOTAL
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)      =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)      =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             TOON_AK(L,J) = SQRT(ABS(B1(L,J)*B1(L,J) - B2(L,J)*B2(L,J)))
             GAMI(L,J)    =  B2(L,J)/(B1(L,J) + TOON_AK(L,J))
             EE1(L,J)     =  EXP(-TOON_AK(L,J)*TAUL(L,J))

             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)  !e1
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J) !e2
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)       !e3
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)       !e4

  15  CONTINUE

!
!     WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
!     L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
!     AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
!     FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
!     ANGLE OR TEMPARATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD.
!
      J                 =  0
      DO 18 JD          =  2,JN,2
         J              =  J + 1
         DO 18 L        =  MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
!          HERE ARE THE EVEN MATRIX ELEMENTS
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)  = EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  18  CONTINUE
      J                 =  0
      DO 19 JD          =  2,JN2,2
         J              =  J + 1
         DO 19 L        =  NSOL+1,NTOTAL
!          HERE ARE THE EVEN MATRIX ELEMENTS
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)  = EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  19  CONTINUE

!
!     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
!     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
!     NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.



      DO 20 L        = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLE) = EL1(L,NLAYER)-RSFX(L)*EL2(L,NLAYER)
         BF(L,JDBLE) = EM1(L,NLAYER)-RSFX(L)*EM2(L,NLAYER)
         EF(L,JDBLE) = 0.0
  20  CONTINUE
      DO 21 L        = NSOL+1,NTOTAL
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLEDBLE) = EL1(L,NDBL)-RSFX(L)*EL2(L,NDBL)
         BF(L,JDBLEDBLE) = EM1(L,NDBL)-RSFX(L)*EM2(L,NDBL)
         EF(L,JDBLEDBLE) = 0.0
  21  CONTINUE


      RETURN
      END
