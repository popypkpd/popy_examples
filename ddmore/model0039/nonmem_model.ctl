$PROBLEM    MODEL FOR INSULIN,GLUCOSE and HOT GLUCOSE after IVGTT
$INPUT      OID TIME AMT ODV INSU TOTG CMT BW EVID AXID DAY DROP DROP
            DROP FLAG TRAC=DROP DROP DROP DROP RATE DROP DROP DROP DV
            SRT3 ED2 EILG=DROP ID TYPE FRST=DROP
;  DV contains insulin, glucose as well as hot glucose
$DATA ../Simulated_IVGTT_IGI_commas.csv 
    IGNORE=# 
    ;IGNORE=(ID.NE.15)
    ;IGNORE=(TYPE.GT.3)
    ;IGNORE=(TIME.GT.24)
    IGNORE=(CMT.EQ.0)
    ;IGNORE=(CMT.EQ.1)
    ;IGNORE=(CMT.EQ.2)
    ;IGNORE=(CMT.EQ.3)
    ;IGNORE=(CMT.EQ.4)
    ;IGNORE=(CMT.EQ.5)
    ;IGNORE=(CMT.EQ.6)
    ;IGNORE=(CMT.EQ.7)
    ;IGNORE=(CMT.EQ.8)
    ;IGNORE=(CMT.EQ.9)
    
    ; IGNORE(ID.EQ.9)IGNORE(ID.EQ.12)

$ABBREVIATED DERIV2=NO
$SUBROUTINE ADVAN6 TRANS1 TOL=5
$MODEL      COMP=TOT_C ;1
            COMP=INS_C ;2
            COMP=HOT_C ;3
            COMP=TOT_E ;4 effect compartment for glucose
            COMP=PER_C ;5
            COMP=PER_I ;6
            COMP=PER_H ;7
            COMP=INS_E ;8 effect compartment for insulin
            COMP=INS_S ;9 secretion compartment for insulin


$THETA  1 FIX ; 1 Not in use
$THETA  1 FIX ; 2 Not in use
$THETA  (0,9.16612318909858) ; 3 VG  (L)
$THETA  (0,0.0856544825911882) ; 4 KG  (/min)
$THETA  (0,1.2179912083266) ; 5 KI  (/min)
$THETA  (0,6.71024822654209) ; 6 KGI (/min/pM)
$THETA  (0,0.0590522068938238) ; 7 KEOG(/min) (0.08 for HV)
$THETA  (0,661.784806539811) ; 8 IFST(pM/[mg/dl]) note:(0,3) for HV
$THETA  (0,0.472640233864405) ; 9 Q
$THETA  (0,8.55073592707769) ; 10 VP
$THETA  (0,0.0284945031485903) ; 11 K26
$THETA  (0,.000799) FIX ; 12 K62
$THETA  (0,0.021637630446828) ; 13 KEOI
$THETA  (0,1.43134552100246,5) ; 14 IPRG
$THETA  (3,6.04643669479359) ; 15 VI
$THETA  1 FIX ; 16 BASG_FAC; D2 - duration of ins. infusion
$THETA  1 FIX ; 17 BASI FACILAG INSULIN DOSE
$THETA  0 FIX ; 18 GPRI
$THETA  (0,0.0435584975348915) ; 19 TOTG PROP
$THETA  (0,0.253989394960252) ; 20 INSU PROP
$THETA  (0,0.105393878307098) ; 21 HOTG PROP
$THETA  (0,3.29549314960977) ; 22 EARLY ERROR
$THETA  (-5,-2.88150187087608,0) ; 23 GPRG
$THETA  (0,0.3751721008462) ; 24 KIS
$THETA  0 FIX ; 25 GELI
;--------patients-------
$THETA  (0,0.0283941818021748) ; 26 KG  (/min)
$THETA  (0.5,0.892501172846716,2) ; 27 BASG_FACpat
$THETA  (0,2.9744817285232) ; 28 KGI (/min/pM)
$THETA  0 FIX ; 29 IFST(pM/[mg/dl]) note:(0,3) for HV
$THETA  (.5,0.930501518914744,2) ; 30 BASI_FACpat
$THETA  (0,0.0510225189046737) ; 31 RESHpat
$THETA  (0,15.1573029861128) ; 32 CLGIhvins
$THETA  (0,5.01845628714497) ; 33 CLGIclamp
$THETA  0 FIX ; ; 34 IPRG
$THETA  0 FIX ; 35 GPRI
$THETA  0 FIX ; 36 GPRG
$THETA  0 FIX ; 37 GELI
;---------
$OMEGA  1.00  FIX  ;      1 GSS
$OMEGA  1.00  FIX  ;      2 ISS
$OMEGA  BLOCK(3)
 0.0910824403788188  ;       3 VG
 -0.198169492436692 0.752054276962733  ;        4 Q
 0.0807784065062901 -0.0939603544874258 0.173060574312854  ;       5 VI
$OMEGA  0.0488457040826792  ; 6 CLGI hv+clamp+hvins
$OMEGA  0  FIX  ; 7 CLGI hvins
$OMEGA  0.424782381451217  ;     8 IFST
$OMEGA  0.362564306724732  ;       9 KG
$OMEGA  0.0885639396517315  ;     10 VP_
$OMEGA  0.608902226443451  ;    11 K26_
$OMEGA  0.228265906658762  ; 12 CLGI pat
$OMEGA  0.350734462920572  ;   13 KEOI_
$OMEGA  0.128100473310518  ;   14 IPRG_
$OMEGA  0.0829671896846808  ;      15 KI
$OMEGA  0  FIX  ;      16 D2
$OMEGA  0  FIX  ;    17 ILAG
$OMEGA  0  FIX  ;   18 GPRI_
$OMEGA  0  FIX  ; 19 REST can't be used with HYBRID
$OMEGA  0  FIX  ; 20 RESI can't be used with HYBRID
$OMEGA  0  FIX  ; 21 RESH can't be used with HYBRID
$OMEGA  0  FIX  ; 22 RESE can't be used with HYBRID
$OMEGA  0  FIX  ;   23 GPRG_
$OMEGA  0  FIX  ;    24 KIS_

$SIGMA  1  FIX

$PK   
    ;"FIRST
    ;"      COMMON/PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
    ;"INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
    ;"IMAX=5000000
    
    ;---------------------Handling baseline values--------------
    IF(NEWIND.NE.2) BASG=TOTG
    IF(NEWIND.NE.2) BASI=INSU
    
    GSS = BASG*THETA(16) ;;*EXP(ETA(1)*THETA(19))
    ISS = BASI*THETA(17) ;;*EXP(ETA(2)*THETA(20))

    IF(TYPE.EQ.2) GSS = BASG*THETA(27) ;;*EXP(ETA(1)*THETA(19))    ; theta(27) and theta(30) are scaling factors as in patients the baseline was not coming back to the initial value
    IF(TYPE.EQ.2) ISS = BASI*THETA(30) ;;*EXP(ETA(2)*THETA(20))
    
    IF(TYPE.EQ.3) ISS = INSU  
    
    ;----------------------Initializing other parameters--------
    ;
    VG  = THETA(3)*BW/70 ;;*EXP(ETA(3))     ;Central volume for glucose
    CLG  = THETA(4) ;;*EXP(ETA(9))            ;k10 for glucose disappear
    CLI  = THETA(5) ;;*EXP(ETA(15))            ;k10 for insulin disappear
    CLGI = THETA(6)/1000 ;;*EXP(ETA(6))       ;2nd order rate constant for insulin-dep. glucose elimination

    IF(TYPE.EQ.3) CLGI = THETA(33)*EXP(ETA(6))/1000
    ;IF(TYPE.EQ.4) CLGI = THETA(32)*EXP(ETA(6))/1000

    KEOG= THETA(7);*EXP(ETA(7))            ;keo for glucose 
    IFST= THETA(8) ;;*EXP(ETA(8))            ;primary insulin response (lacking for patients)
    Q   = THETA(9) ;;*EXP(ETA(4))            ;k12 for glucose
    VP  = THETA(10) ;;*EXP(ETA(10))          ;k21 for glucose
    K15 = Q/VG
    K51 = Q/VP
    K26 = THETA(11) ;;*EXP(ETA(11))          ;k12 for insulin
    K62 = THETA(12) ;*EXP(ETA(12))          ;k21 for insulin
    KEOI= THETA(13) ;;*EXP(ETA(13))          ;keo for insulin
    IPRG= THETA(14) ;;*EXP(ETA(14))          ;influence of glucose on insulin production
    VI  = THETA(15)*BW/70 ;;*EXP(ETA(5))    ;central volume for insulin
    ;ILAG= THETA(17)*EXP(ETA(17))          ;lag-time for insulin infusion
    GPRI= THETA(18) ;;*EXP(ETA(18))          ;influence of insulin on glucose production
    REST= THETA(19) ;;*EXP(ETA(19))          ;CV for total glucose residual error
    RESI= THETA(20) ;;*EXP(ETA(20))          ;CV for       insulin residual error
    RESH= THETA(21) ;;*EXP(ETA(21))          ;CV for hot   glucose residual error
    RESE= THETA(22) ;;*EXP(ETA(22))          ;correction factor for early residual error
    GPRG= THETA(23) ;;*EXP(ETA(23))          ;influence of glucose on glucose production
    KIS = THETA(24) ;;*EXP(ETA(24))          ;k92 "ka" for insulin
    GELI= THETA(25)                       ;

    IF (TYPE.EQ.2) THEN
    ;----------------------Initializing other parameters--------
        CLG  = THETA(26) ;;*EXP(ETA(9))            ;k10 for glucose disappear
        ;  VG  = THETA(27)*EXP(ETA(3))            ;VG for patients
        CLGI = THETA(28)/1000 ;;*EXP(ETA(12))       ;2nd order rate constant for insulin-dep. glucose elimination
        IFST= THETA(29);;*EXP(ETA(8))            ;primary insulin response (lacking for patients)
        RESH = THETA(31)                       ;RESH for glucose
        ; K51 = THETA(31);;*EXP(ETA(10))          ;k21 for glucose
        ; K26 = THETA(32);;*EXP(ETA(11))          ;k12 for insulin
        ; K62 = THETA(33);;*EXP(ETA(12))          ;k21 for insulin
        ; IPRG= THETA(34);;*EXP(ETA(14))          ;influence of glucose on insulin production
        GPRI= THETA(35);;*EXP(ETA(18))          ;influence of insulin on glucose production
        GPRG= THETA(36);;*EXP(ETA(23))          ;influence of glucose on glucose production
        GELI= THETA(37)                       ;
    ENDIF
    
    ;----------------production and initializing compartments------------
    KG  = CLG/VG                          ;k10 for glucose disappear
    KGI = CLGI/VG                         ;2nd order rate constant for insulin-dep. glucose elimination
    KI  = CLI/VI                          ;duration of insulin infusion
    GPRO= GSS*(KG+KGI*ISS)*VG             ;baseline glucose production
    IPRO= KI*ISS*VI                       ;baseline insulin production

    ;IF (A_0FLG.EQ.1) THEN 
    ;    IF (TYPE.EQ.2) A_0(9)=IFST
    ;ENDIF
        
    F1  = BW/10                            ;correction for glucose in mg/dl
    IF(AMT.EQ.1) F1=GSS*VG
    
    F2  = 0
    IF(SRT3.EQ.0) F2=ISS*VI
    
    F3  = BW/10        
    F4  = GSS
    F5  = K15*GSS*VG/K51
    F6  = GSS       ;K26*ISS*VI/K62
    F7  = 0
    F8  = ISS
    F9  = IFST

    IF(TYPE.EQ.2) F1  = 1/10                            ;correction for glucose in mg/dl
    IF(TYPE.EQ.2.AND.AMT.EQ.1) F1=GSS*VG                         
    IF(TYPE.EQ.2)   F2  = ISS*VI
    IF(TYPE.EQ.2.AND.TIME.GT.0) F2=1
    IF(TYPE.EQ.2)  F3  = 0.26/10
    
    IINF = 0
    IF(TYPE.EQ.3)  F1  = 1/1000
    IF(TYPE.EQ.3.AND.AMT.EQ.1)  F1  = GSS*VG
    IF(TYPE.EQ.3)  F3  = 1/1000
    IF(TYPE.EQ.3)  F9  = 0
    IF(TYPE.EQ.3.AND.TIME.GT.150) IINF = 1*BW

    IF(TYPE.EQ.4)  F1  = BW/10
    IF(TYPE.EQ.4.AND.AMT.EQ.1)  F1  = GSS*VG
    IF(TYPE.EQ.4)  F3  = BW
    IF(TYPE.EQ.4.AND.TIME.GT.0) F2=1

$DES  
    PTOT   = K15*A(1)-K51*A(5)             ;distribution of tot glucose
    ;PINS   = K26*A(2)-K62*A(6)             ;distribution of     insulin
    PHOT   = K15*A(3)-K51*A(7)             ;distribution of hot glucose
    IGPR   = (A(6)/GSS+.0001)**IPRG              ;glucose on insulin production
    GIPR   = (A(8)/ISS+.0001)**GPRI              ;insulin on glucose production
    GGPR   = (A(4)/GSS+.0001)**GPRG              ;glucose on glucose production
    GPRDG   = GPRO*GIPR*GGPR                ;glucose production
    IPRDG   = IPRO*IGPR                     ;insulin production
    GIEL   = (A(8)/ISS+.0001)**GELI       ;insulin on glucose elimination
    ;TELI   = (KG + KGI*A(8)*GIEL)*A(1)          ;tot glucose elimination
    ;IELI   =  KI*A(2)                           ;    insulin elimination
    HELI   = (KG + KGI*A(8)*GIEL)*A(3)          ;hot glucose elimination
    
    IF(TYPE.EQ.3) IPRD=IPRO
    
    DADT(1) = GPRDG - PTOT - KG*A(1) - (KGI*A(8)*GIEL)*A(1)
    DADT(5) = PTOT

    DADT(2) = IINF + IPRDG + KIS*A(9) - KI*A(2)
    DADT(9) =              - KIS*A(9)
    DADT(3) = - HELI - PHOT
    DADT(7) =          PHOT

    DADT(4) = KEOG*(A(1)/VG - A(4))        ;tot glucose effect compartment
    DADT(6) = K26*(A(1)/VG - A(6))        ;tot glucose effect compartment
    DADT(8) = KEOI*(A(2)/VI - A(8))        ;insulin effect-compartment

$ERROR  
    CTOT    =A(1)/VG
    CINS    =A(2)/VI
    CHOT    =A(3)/VG
    
    ;TEFF    =A(4)
    ;TPER    =A(5)
    ;IPER    =A(6)
    ;HPER    =A(7)
    ;IEFF    =A(8)
    ;ISEC    =A(9)

    ;----------additional calculations----------- 
    ;GRAT = TEFF/GSS                            ;glucose eff comp to GSS ratio (concentrations) 
    ;IRAT = IEFF/ISS                            ;insulin eff comp to ISS ratio (concentrations)
    ;PRAT = IPER/(K26*ISS*VI/K62)               ;insulin per comp to initial conc ratio (amounts)
    ;GLEF = (GRAT+.0001)**IPRG                  ;glucose eff on insulin prod 
    ;INPR = (KI*ISS*VI)*(A(4)/GSS+.0001)**IPRG  ;total insulin production 
    ;PRO3 = (INPR+K62*IPER)                ;"actual" production
    ;PRO4 = INPR/(INPR+K62*IPER)                ;"actual" production
    ;ELIM = KI*CINS                             ;insulin elimination
    ;IOUT = K26*(CINS-ISS)                       ;delta flow out of comp 2
    ;IIN  = K62*(IPER-K26*ISS*VI/K62)            ;delta flow out of comp 6
    ;PNET = IOUT-IIN                            ;net flow in comp 2

    IF (CMT.EQ.1) IPRED = CTOT+.00001
    IF (CMT.EQ.2) IPRED = CINS+.00001
    IF (CMT.EQ.3) IPRED = CHOT+.00001
    IF (CMT.GE.4) IPRED = 1

    IF (CMT.EQ.1) W = REST
    IF (CMT.EQ.2) W = RESI
    IF (CMT.EQ.3) W = RESH
    IF (CMT.GE.4) W = 1
    
    IF(TIME.LE.2) W = W*RESE

    IRES  = DV-LOG(IPRED+1e-6)
    IWRES = IRES/W
    Y     = LOG(IPRED+1e-6) + ERR(1)*W
    LOGIPRED = LOG(IPRED+1e-6)
    OBS   = ODV
    
    OBJVt = IWRES**2 + 2*LOG(W)

;$SIM (87568767)

$ESTIMATION MAXEVAL=0 PRINT=1 POSTHOC MSFO=msfb592 METH=1

;$COVARIANCE PRINT=R

;$TABLE ID TIME TYPE CMT GRAT IRAT PRAT GLEF INPR PRO3 

$TABLE 
    ID TIME TYPE EVID CMT 
    F1 F2 F3 F4 F5 F6 F7 F8 F9
    GPRO GPRI GPRG 
    GSS ISS VG VI
    KG KGI KI KEOI GELI KIS
    CTOT W IPRED IRES IWRES OBJVt
    NOPRINT file=pttab.csv

;   PRO4 PNET                                       NOPRINT ONEHEAD FILE=mytab590
;$TABLE      ID TIME TYPE GPRDG IPRDG CTOT CINS CHOT TEFF IEFF TPER IPER HPER ISEC CMT NOPRINT ONEHEAD FILE=mytab592
;$TABLE      ID TIME DV IPRED PRD OBS CTOT CINS CHOT NOAPP NOPRINT ONEHEAD FILE=mytab592
;$TABLE      ID TIME AMT RATE CMT DV BW TYPE IPRED IWRES NOPRINT ONEHEAD FILE=sdtab592
;$TABLE      ID ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA8 ETA9 ETA10 ETA11 ETA12 ETA13 ETA14 ETA15 NOPRINT ONEHEAD FILE=patab592
;$TABLE      ID ETA3 ETA4 ETA5 ETA6 ETA8 ETA9 ETA10 ETA11 ETA12 ETA13 ETA14 ETA15 FIRST NOAPP NOPRINT ONEHEAD FILE=patab592

