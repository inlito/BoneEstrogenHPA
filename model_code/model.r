model <- function(t, y, parms, stress, EST_group, renal_state, PTH_state) {
  PTH <- y[1]
  S <- y[2]
  PTmax <- y[3]
  B <- y[4]
  SC <- y[5]
  A  <- y[6]
  P <- y[7]
  ECCPhos <- y[8]
  T <- y[9]
  R <- y[10]
  HAp <- y[11]
  OBfast <- y[12]
  OBslow <- y[13]
  PhosGut <- y[14]
  IntraPO <- y[15]
  OC <- y[16]
  ROB1 <- y[17]
  TGFB <- y[18]
  TGFBact <- y[19]
  L <- y[20]
  RNK <- y[21]
  M <- y[22]
  N <- y[23]
  O <- y[24]
  Q <- y[25]
  Qbone <- y[26]
  RX2 <- y[27]
  CREB <- y[28]
  BCL2 <- y[29]
  TERISC <- y[30]
  BMDfn <- y[31]


# hpa compartments --------------------------------------------------------

  CRH <- y[32]
  ACTH <- y[33]
  cortisol <- y[34]
  corticotroph_mass <- y[35]
  adrenal_mass <- y[36]
  GR_resistance <- y[37]
  CR <- y[38]
  external_CRH <- y[39]
  
# parameter functions --------------------------------------------------------
  
  GR <- function(cortisol, kGR) {
    return(1 / ((cortisol/ kGR) ** 3 + 1))
  }
  
  MR <- function(cortisol) 1/cortisol
  
  external_crh <- function(t) return(0)
  # external_crh <- function(t) return(ifelse(t < 30, 20, 0))
  
  stress_test <- function(t, u0, stress=NULL) {
    if (is.null(stress)) {
      return(u0)
    }
    return(ifelse(stress$start < t && t < stress$end, 4, 1))
  }
 
  # Glomerular filtration rate basec on renal state
  get_GFR <- function(t, GFR0, renal_state=NULL) {
    if (is.null(renal_state))  {
      return(GFR0)
    } else if (renal_state == "progressive")  {
      GFR0 <- GFR0 * exp(-t * 1.862775 * 1e-5)   # fig 6a) - progressive renal insufficiency
    } else if (renal_state == "immediate")  {
      GFR0 <- GFR0/6   # fig 6b) - immediate renal insufficiency
    }
    return(GFR0)
  }
  
  # PTH based on disease state
  get_PTH <- function(t, PTH_state=NULL) {
    if (is.null(PTH_state))  {
      return(0)
      # fig 5 - primary hyperparathyroidism through a longitudinal increase in endogenous PTH production
    } else if (PTH_state == "primary_hyper_long")  {
      PTH_conc <- 5 * tanh(3*t/(12*30*24))
    } else if (PTH_state == "primary_hyper_teri") {
      PTH_conc <- 0
    } else if (PTH_state == "primary_hypo") {
      # PTH_conc <- -0.5 * tanh(3*t/(12*30*24))    ("divide EPTH by 2"?)
      PTH_conc <- -0.5
    }
    return(PTH_conc)
  }
  
  
  
  # get EST concentration in pmol to match PTH unit in Riggs.
  get_EST <- function(t, EST_group){
    cycle_len <- 28 # cycle length
    day <- t/24 # time is in hours, convert to day
    if (EST_group == "baseline")  {
      EST_conc <- 363
    } else if (EST_group == "postmeno_E2")  {
      EST_conc <- 26 # menopause: <26 pmol/L https://pubmed.ncbi.nlm.nih.gov/30981845/
    } else if (EST_group == "OC_EE")  {
      EST_conc <- 407 # OC with ethinyl estradiol - cumulative exposure is 11.3 nmol/L for a cycle = 11300 pmol/L --> distributed over 28 days  --> 407 pmol/L
      # E2_conc <- 120 # OC without ethinyl estradiol: 110-183 pmol/L https://link.springer.com/article/10.1007/s00198-019-05103-6
      # E2_conc <- 120 # "All OC users used combined ethinylestradiol and progesterone contraception," OC with EE https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8606688/
    } else if (EST_group == "premeno_E2") {
      EST_conc <- 112 + 
        100 * exp(-((day%%cycle_len-(cycle_len-17))**2)/cycle_len) + 
        700 * exp(-((day%%cycle_len-(cycle_len-15))**2)/4) + 
        500 * exp(-((day%%cycle_len-(cycle_len-5))**2)/12)
    }
    return(EST_conc * parms$V1) # volume "Vascularate" compartment from Riggs-model
  }

  
  
  Hplus <- function(S, T, n)
    return(S**n/(T**n + S**n))
  
  Hminus <- function(S, T, n) 
    return(T**n/(T**n + S**n))


  yn <- c()
  with(parms, {
    T13 = (CaDay / 24) / Q0
    
    T15 = CaDay / (2.35 * 14 * 24)
    
    
    
    Osteoclast =  OC
    
    J14OC50 = exp(log((J14OCmax * OC0 ^ J14OCgam / T13) - OC0 ^ J14OCgam) / J14OCgam)
    
    OCeqn = (J14OCmax * Osteoclast ^ J14OCgam) / (Osteoclast ^ J14OCgam + J14OC50 ^ J14OCgam)
    
    kinRNK = (koutRNK * RNK0 + k3 * RNK0 * RANKL0 - k4 * M0) / TGFBact0 ^ kinRNKgam
    
    
    MOCratio = M / Osteoclast
    
    MOCratio0 = M0 / OC0
    
    MOCratioEff = (MOCratio / MOCratio0) ^ MOCratioGam
    
    
    J14OCdepend = OCeqn * Q0 * FracJ14 * MOCratioEff
    
    J14 = T13 * Q0 * (1 - FracJ14) + J14OCdepend
    
    J41 = 0.464 * J14
    
    
    PicOCkin = Pic0
    
    bigDb = kb * OB0 * Pic0 / ROB0
    
    
    kinTGF = koutTGF0 * TGFB0
    
    koutTGF = koutTGF0
    
    koutTGFact = koutTGF0 * 1000
    
    
    E0PicROB = FracPicROB * Pic0
    
    EC50PicROBparen = (EmaxPicROB * TGFBact0 ^ PicROBgam / (Pic0 - E0PicROB)) - TGFBact0 ^
      PicROBgam
    
    EC50PicROB = exp(log(EC50PicROBparen) / PicROBgam)
    
    
    Dr = kb * OB0 / Pic0
    
    PicROB = E0PicROB + EmaxPicROB * TGFBact ^ PicROBgam / (TGFBact ^ PicROBgam + EC50PicROB ^
                                                              PicROBgam)
    
    ROBin2 = Dr * PicROB
    
    ROBin = ROBin2
    
    
    E0PicOB = FracPicOB * Pic0
    
    EC50PicOBparen = (EmaxPicOB * TGFBact0 ^ PicOBgam / (Pic0 - E0PicOB)) - TGFBact0 ^ PicOBgam
    
    EC50PicOB = exp(log(EC50PicOBparen) / PicOBgam)
    
    
    PicOB = E0PicOB + EmaxPicOB * TGFBact ^ PicOBgam / (TGFBact ^ PicOBgam + EC50PicOB ^ PicOBgam)
    
    KPT = 1 * (bigDb / PicOB)
    
    
    D = ROB1
    
    EC50MeffOC = exp(log(M0 ^ kinOCgam * EmaxMeffOC / (1 - E0Meff) - M0 ^ kinOCgam) / kinOCgam)
    
    MeffOC = E0Meff + (EmaxMeffOC * M ^ kinOCgam / (M ^ kinOCgam + EC50MeffOC ^ kinOCgam))
    
    kinOC2 = Da * PicOCkin * MeffOC * OC0
    
    
    E0PicOC = FracPicOC * Pic0
    
    EC50PicOCparen = (EmaxPicOC * TGFBact0 ^ PicOCgam / (Pic0 - E0PicOC)) - TGFBact0 ^ PicOCgam
    
    EC50PicOC = exp(log(EC50PicOCparen) / PicOCgam)
    
    PicOC = E0PicOC + ((EmaxPicOC * TGFBact ^ PicOCgam) / (TGFBact ^ PicOCgam + EC50PicOC ^ PicOCgam))
    
    
    PiL0 = (k3 / k4) * RANKL0
    
    PiL = M / 10
    
    
    EC50survInPar = (E0RANKL - EmaxL) * (PiL0 ^ LsurvOCgam / (E0RANKL - 1)) - PiL0 ^ LsurvOCgam
    
    EC50surv = exp(log(EC50survInPar) / LsurvOCgam)
    
    LsurvOC = E0RANKL - (E0RANKL - EmaxL) * (PiL ^ LsurvOCgam / (PiL ^ LsurvOCgam + EC50surv ^ LsurvOCgam))
    
    KLSoc = Da * PicOC * LsurvOC #*3
    
    
    C4 = PTH / V1
    
    T66 = (T67 ^ AlphOHgam + 3.85 ^ AlphOHgam) / 3.85 ^ AlphOHgam
    
    k15a = k14a * QboneInit / Q0
    
    J14a = k14a * Qbone
    
    J15a = k15a * Q
    
    
    kLShap = 1 / HApMRT
    
    kHApIn = kLShap / OB0
    
    
    J15 = T15 * P * (1 - FracJ15) + T15 * P * FracJ15 * HAp
    
    J42 = 0.464 * J15
    
    
    OBfast0 = OB0 * FracOBfast
    
    Osteoblast = OBfast + OBslow
    
    OsteoEffect = (Osteoblast / OB0) ^ OsteoEffectGam
    
    PTH50 = EmaxLpth * 3.85 - 3.85
    
    
    PTHconc = C4
    
    LpthEff = EmaxLpth * (PTHconc) / ((PTH50 * OsteoEffect ^ TESTPOWER) + (PTHconc))
    
    
    kinLbase = koutL * RANKL0
    
    kinL = kinLbase * (OsteoEffect) * LpthEff
    
    pObase = kO * OPG0
    
    
    pO = pObase * (D / ROB0) * ((PTHconc + (opgPTH50 * (D / ROB0))) / (2 * PTHconc)) + IO
    
    RX2Kin = RX2Kout0 * RX20
    
    
    EC50PTHRX2x = ((EmaxPTHRX2x * 3.85) / (RX2Kout0 - E0rx2Kout)) - 3.85
    
    RX2Kout = E0rx2Kout + EmaxPTHRX2x * PTHconc / (PTHconc + EC50PTHRX2x)
    
    
    EC50PTHcreb = ((EmaxPTHcreb * 3.85) / (1 - E0crebKin)) -  3.85
    
    crebKin0 = crebKout * CREB0
    
    
    crebKin = crebKin0 * (E0crebKin + EmaxPTHcreb * PTHconc / (PTHconc + EC50PTHcreb))
    
    bcl2Kin = RX2 * CREB * 0.693
    
    
    CaConc = P / 14
    
    C2 = ECCPhos / V1
    
    
    PO4inhPTH = (C2 / 1.2) ^ PO4inhPTHgam
    
    
    PhosEffTop = (PhosEff0 - 1) * (1.2 ^ PhosEffGam + PhosEff50 ^ PhosEffGam)
    
    PhosEffBot = PhosEff0 * 1.2 ^ PhosEffGam
    
    PhosEffMax =  PhosEffTop / PhosEffBot
    
    PhosEff = PhosEff0 - (PhosEffMax * PhosEff0 * C2 ^ PhosEffGam / (C2 ^ PhosEffGam  + PhosEff50 ^ PhosEffGam))
    
    PhosEffect <- ifelse(C2 > 1.2, PhosEff, 1)
    
    
    T68 = T66 * C4 ^ AlphOHgam / (T67 ^ AlphOHgam * PO4inhPTH + C4 ^ AlphOHgam)
    
    SE = T65 * T68 * PhosEffect
    
    C8 = B / V1
    
    C1 = P / V1
    
    T36 = T33 + (T34 - T33) * (C8 ^ CaPOgam / (T35 ^ CaPOgam + C8 ^ CaPOgam))
    
    T37 = T34 - (T34 - T33) * (C8 ^ CaPOgam / (T35 ^ CaPOgam + C8 ^ CaPOgam))
    
    
    # # RENAL CALCIUM HANDLING
    # Calcium filtration rate in the kidney
    # Assume 50/50 PTH-independent PTH-dependent reabsorption
    # 
    
    CaFilt = 0.6 * 0.5 * get_GFR(t, GFR0, renal_state) * C1
    
    
    # Maximum calcium reabsorption in the kidney - PTH sensitiv
    ReabsMax = (0.3 * get_GFR(t, GFR0, renal_state) * 2.35 - 0.149997) * (Reabs50 + 2.35) / 2.35
    
    # Effect of PTH on calcium reabsorption
    T17 = 3.85 * T16 - 3.85
    
    ReabsPTHeff = (T16 * C4) / (C4 + T17)
    
    
    # PTH-sensitive calcium reabsorption in kidney
    # Reabs50 = 1.573 = H(4-u)-delta 
    
    CaReabsActive = (ReabsMax * C1 / (Reabs50 + C1)) * ReabsPTHeff
    
    T20 = CaFilt - CaReabsActive
    
    T10 = T7 * C8 / (C8 + T9)
    
    # Temporary calcium excretion rate. J27 will be the flux of calcium out of the plasma via the kidney
    
    J27 = max((2 - T10) * T20, 0)
    
    
    ScaEff =  (2.35 / CaConc) ^ ScaEffGam
    
    T72 = 90 * ScaEff
    
    T73 = T71 * (C8 - T72)
    
    T74 = (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))
    
    T75 = T70 * (0.85 * (1 + T74) + 0.15)
    
    T76 = T70 * (0.85 * (1 - T74) + 0.15)
    
    
    # phosphate renal excretion
    T47 = T46 * 0.88 * get_GFR(t, GFR0, renal_state)
    
    J48 = max(0.88 * get_GFR(t, GFR0, renal_state) * C2 - T47, 0)
    
    # phosphate oral absorption
    J53 = T52 * PhosGut
    
    J54 = T49 * C2
    
    J56 = T55 * IntraPO
    
    #  Parameters describing TGF-beta effects on Osteoblast and osteoclast differentiation and apoptosis
    E0PicOBkb = MultPicOBkb * Pic0
    
    EmaxPicOBkb = FracPic0kb * Pic0
    
    EC50PicOBparenKb = ((E0PicOBkb - EmaxPicOBkb) * TGFBact0 ^ PicOBgamkb) / (E0PicOBkb - Pic0)  - TGFBact0 ^ PicOBgamkb
    
    EC50PicOBkb = exp(log(EC50PicOBparenKb) / PicOBgamkb)
    
    PicOBkb = E0PicOBkb - (E0PicOBkb  - EmaxPicOBkb) * TGFBact ^ PicOBgamkb / (TGFBact ^ PicOBgamkb + EC50PicOBkb ^ PicOBgamkb)
    
    # Estrogen effect that propogates through to OB apoptosis
    PicOBkbEff = PicOBkb / Pic0  # *(1/(pow(EST,ESTscalePicB1)))
    
    # Parameters describing osteoblast apoptosis as affected by PTH (continuous vs intermitent)
    E0RUNX2kbEff = E0RUNX2kbEffFACT * kb
    
    RUNX2 <- ifelse(BCL2 > 105, BCL2 - 90, 10)
    
    RUNkbMax = E0RUNX2kbEff * RUNkbMaxFact
    
    INparen = (RUNkbMax * RUNX20 ^ RUNkbGAM) / (E0RUNX2kbEff - kb) - RUNX20 ^ RUNkbGAM
    
    RUNkb50 = exp(log(INparen) / RUNkbGAM)
    
    RUNX2kbPrimeEff = RUNkbMax * RUNX2 ^ RUNkbGAM / (RUNX2 ^ RUNkbGAM + RUNkb50 ^ RUNkbGAM)
    
    kbprime = E0RUNX2kbEff * PicOBkbEff - RUNX2kbPrimeEff
    
    kbslow = kbprime * Frackb
    
    kbfast = (kb * OB0 + kbslow * OBfast0 - kbslow * OB0) / OBfast0
    
    Frackb2 = kbfast / kbprime
    
    
    # Equations relating to calcium movement to/from the gut
    
    T29 = (T28 * T0 - 0.17533 * T0) / 0.17533
    
    T31 = T28 * T / (T + T29)
    
    # R is calcitriol-dependent gut Ca2+ absorption
    T83 = R / 0.5
    
    # J40 = calcium flux from gut to plasma 
    J40 = T31 * T * T83 / (T + T81) + T87 * T
    
    # T85 relates to extent of absorption of orally-administered dose
    T85Rpart = R ^ T80 / (R ^ T80 + T81 ^ T80)
    
    T85 = T77 * T85Rpart
    
    F11 = T85
    
    # Calcitriol equations
    
    INparenCtriol = ((CtriolMax - CtriolMin) * C8 ^ CtriolPTgam) / (CtriolMax - 1) - C8 ^ CtriolPTgam
    
    Ctriol50 = exp(log(INparenCtriol) / CtriolPTgam)
    
    CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * C8 ^ CtriolPTgam / (C8 ^ CtriolPTgam + Ctriol50 ^ CtriolPTgam)
    
    PTin = PTout * CtriolPTeff
    
    # S is the PTH gland pool
    FCTD = (S / 0.5) * PTmax
    
    
    INparenCa = (T58 - T61) * 2.35 ^ T59 / (T58 - 385) - 2.35 ^ T59
    
    T60 = exp(log(INparenCa) / T59)
    
    T63 =  T58 - (T58 - T61) * (CaConc) ^ T59 / ((CaConc) ^ T59 + T60 ^ T59)
    
    
    EPTH = T63 * FCTD * (1 + get_PTH(t, PTH_state))
    
    IPTH = 0.693 * SC +  IPTHinf
    
    
    SPTH = EPTH + IPTH
    
    
    kout = T57 / 14
    
    TERIPK = TERISC * TERICL / TERIVC
    
    
    
    gamOCfn = gamOCfnBAS   
    
    koutBMDfn = koutBMDfnBAS 
    
    kinBMDfn =  koutBMDfn*BMDfn_0
    

    # hpa calculations --------------------------------------------------------
    
    baseline_EST <- get_EST(t, EST_group="baseline")
    
    koutTGFeqn = koutTGF * TGFB * ((Osteoclast / OC0) ^ OCtgfGAM) * (0.5 + 1.5*Hplus(get_EST(t, EST_group), baseline_EST, 2))
    
    
    
# DIFFERENTIAL EQUATIONS ------------------------------------------------------


# ca pth differential equations -----------------------------------------------

        
    ## Parathyroid (PTH)
    yn[1] = SPTH - kout * PTH + TERIPK
    
    ##  S - PT gland pool
    yn[2] = (1 - S) * T76 - (S * T75)
    
    ## PTmax - PT maximum capacity
    yn[3] = PTin - PTout * PTmax
    
    ## B  -  Plasma calcitriol
    yn[4] = A - T69 * B
    
    ## SC - SubQ
    yn[5] = IPTHint - 0.693 * SC
    
    ## A - 1-alpha hydroxylase
    yn[6] = SE - T64 * A
    
    ## P - plasma calcium
    yn[7] = J14 - J15 - J27 + J40
    
    ## ECCphos - phosphate
    yn[8] = J41 - J42 - J48 + J53 - J54 + J56
    
    ## T - oral calcium gut
    yn[9] = OralCa * F11 - J40
    
    ## R - intestinal calcium
    yn[10] = T36 * (1 - R) - T37 * R
    
    ## HAp - Hydroxyapatite
    yn[11] = kHApIn * Osteoblast - kLShap * HAp
    
    ## OBfast - Fast osteoblast
    yn[12] = (bigDb / PicOB) * D * FracOBfast * Frackb2 - kbfast * OBfast * (0.5 + 1.5*Hminus(get_EST(t, EST_group), baseline_EST, 2))
    
    ## OBslow - Slow osteoblast
    yn[13] = (bigDb / PicOB) * D * (1 - FracOBfast) * Frackb - kbslow * OBslow * (0.5 + 1.5*Hminus(get_EST(t, EST_group), baseline_EST, 2))
    
    ## PhosGut - oral phosphate
    yn[14] = OralPhos * F12 - J53
    
    ## IntraPO - intraceulular phosphate
    yn[15] = J54 - J56
    
    ## OC - Osteoclast
    yn[16] = kinOC2 - KLSoc * OC
    
    ## ROB1 - responding osteoblast
    yn[17] = ROBin *(0.5 + 1.5*Hminus(get_EST(t, EST_group), baseline_EST, 2)) - KPT * ROB1
    
    ## TGFB - Latent TGFbeta
    yn[18] = kinTGF * ((Osteoblast / OB0) ^ OBtgfGAM) * (0.5 + 1.5*Hminus(get_EST(t, EST_group), baseline_EST, 2)) - koutTGFeqn
    
    ## TGFBact - active TGFbeta
    yn[19] = koutTGFeqn - koutTGFact * TGFBact
    
    ## L - RANKL
    yn[20] = kinL - koutL * L - k1 * O * L + k2 * N - k3 * RNK * L + k4 * M
    
    ## RNK - RANK
    yn[21] = kinRNK * TGFBact ^ kinRNKgam - koutRNK * RNK - k3 * RNK * L  + k4 * M
    
    ## M - RANK-RANKL
    yn[22] = k3 * RNK * L - k4 * M
    
    ## N - RANKL-OPG
    
    yn[23] = k1 * O * L - k2 * N
    
    ## O - OPG
    yn[24] = pO - k1 * O * L + k2 * N - kO * O
    
    ## Q - Bone Ca (IC)
    yn[25] = J15 - J14 + J14a - J15a
    
    ## Qbone - Bone Ca (non-IC)
    yn[26] = J15a - J14a
    
    ## RX2 - RunX2
    yn[27] = RX2Kin - RX2Kout * RX2
    
    ## CREB
    yn[28] = crebKin - crebKout * CREB
    
    ## BCL2 - Bcl-2
    yn[29] = bcl2Kin - bcl2Kout * BCL2
    
    ## Teriparatide SQ dosing compartment
    yn[30] = -TERISC * TERICL / TERIVC
    
    ## Bone mineral density, femur neck (fetched from OpenBoneMin cpp model)
    yn[31] = kinBMDfn * ((Osteoblast/OB0)**gamOB) - koutBMDfn * ((OC/OC0)**gamOCfn) * BMDfn
    
  
# hpa differential equations --------------------------------------------------
    
    # stim_PTH <- (1 + 0.1*(PTH**2/(PTH**2 + 53**2)))
    # stim_PTH <- 1
    # stim_PTH <- 1 + 0.1 * (PTH - 55)
    
    # CRH with added stimulation
    stim_crh <- CRH + external_CRH
    
    # Stimularory effect of cortisol on PTH? 2021: https://www.frontiersin.org/articles/10.3389/fendo.2021.692722/full
    # "Serum levels of OCN and cortisol, rather than PTH and calcium, are associated with the development of anxiety and depression symptoms in PHPT patients. " 
    # Stimulatory effect of PTH on cortisol levels 2001:  https://journals.physiology.org/doi/full/10.1152/ajpendo.2001.280.2.EST09
    stim_PTH <- 1 + tanh((PTH - 52)/10)
    
    # Increase in total cortisol by ethinyl estradiol
    stim_EE <- ifelse(EST_group =="OC_EE", 2, 1)

# Equations from model code by Karin et al. (added stim_PTH and inh_EST) --------
    
    # CRH
    yn[32] = gamma_CRH * (stim_EE*stress_test(t, u0, stress) * MR(cortisol) * GR(cortisol, kGR/GR_resistance) - CRH)
    
    # ACTH
    yn[33] = gamma_ACTH * (corticotroph_mass * GR(cortisol, kGR/GR_resistance) * stim_crh - ACTH)
    
    # cortisol
    yn[34] = gamma_cortisol * (adrenal_mass * ACTH*stim_PTH - cortisol/CR)
    
    # corticotroph mass
    yn[35] = gamma_corticotroph_mass * corticotroph_mass * (stim_crh - 1)
    
    # adrenal mass
    yn[36] = gamma_adrenal_mass * adrenal_mass * (ACTH - 1)
    
    # GR resistance
    yn[37] = gamma_GR_resistance * (1 - (1 + (cortisol)**2)*GR_resistance)

    # CR
    yn[38] = gamma_CR * (cortisol - CR)

    # external CRH
    yn[39] = gamma_external_CRH * (external_crh(t) - external_CRH)
    
    list(c(yn))
  })
  
}
