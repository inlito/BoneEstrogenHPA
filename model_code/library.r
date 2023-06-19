source("model_code/initial.r")
source("model_code/parameters.r")
source("model_code/compartments.r")
source("model_code/model.r")

load.model <- function() {
  ## Loads functions for
  ## generating lists of compartments, initials,
  ## parameters, and the model function
  ##
  ## Returns: hpacamod list with the following items:
  ##   $init named list of initial values
  ##   $param named list of model parameters
  ##   $cmt ordered vector of compartment names
  ## Also, the following functions are made available:
  ##   initial()
  ##   compartments()
  ##   parameters()
  
  hpacamod <-
    list(
      init = initial(),
      param = parameters(),
      cmt = compartments(),
      model = model
    )
  
  hpacamod
}


derive.init <- function(hpacamod) {
  ## Derive some initial values from other initials and parameters.
  ## Returns: updated model list after calling copy.init()
  
  hpacamod$init$A <- hpacamod$init$B / 10
  hpacamod$init$TGFB <- hpacamod$param$Pic0 * 1000
  hpacamod$init$TGFBact <- hpacamod$param$Pic0
  hpacamod$init$OBfast <- hpacamod$init$OB * hpacamod$param$FracOBfast
  hpacamod$init$OBslow <- hpacamod$init$OB * (1 - hpacamod$param$FracOBfast)
  hpacamod$init$M <- hpacamod$param$k3 * hpacamod$init$RNK * hpacamod$init$L / hpacamod$param$k4
  hpacamod$init$N <- hpacamod$param$k1 * hpacamod$init$O * hpacamod$init$L / hpacamod$param$k2
  
  copy.init(hpacamod)
}

copy.init <- function(hpacamod) {
  ## Copy some initial compartment values into the parameter list.
  ## Returns: hpacamod list with updated parameters
  
  hpacamod$param$Q0 <- hpacamod$init$Q
  hpacamod$param$OC0 <- hpacamod$init$OC
  hpacamod$param$RNK0 <- hpacamod$init$RNK
  hpacamod$param$RANKL0 <- hpacamod$init$L
  hpacamod$param$RNKL0 <- hpacamod$init$L
  hpacamod$param$RNK0 <- hpacamod$init$RNK
  hpacamod$param$OB0 <- hpacamod$init$OB
  hpacamod$param$ROB0 <- hpacamod$init$ROB1
  hpacamod$param$QboneInit <- hpacamod$init$Qbone
  hpacamod$param$OPG0 <- hpacamod$init$O
  hpacamod$param$RX20 <- hpacamod$init$RX2
  hpacamod$param$CREB0 <- hpacamod$init$CREB
  hpacamod$param$M0 <- hpacamod$init$M
  hpacamod$param$TGFBact0 <- hpacamod$init$TGFBact
  hpacamod$param$TGFB0 <- hpacamod$init$TGFB
  hpacamod$param$BMDfn_0 <- hpacamod$init$BMDfn

  hpacamod
}


responses <- function(out, hpacamod) {
  ## Derives some outcomes of interest from raw simulated output
  ##
  ## Arguments: data frame output from lsoda() call and model list
  ## Returns: updated simulation output with some markers derived
  ##          and some amounts of interest expressed as concentrations
  
  ## BSAP = bone-specific alkaline phosphatase
  ## sCTx = serum C-terminal telopeptide of type I collagen
  ## PTHconc = plasma pth concentration (pg/mL)
  ## PTHpM = plasma pth concentration (pM)
  ## CaConc = plasma calcium concentration (mM)
  ## calcitriol = plasma calcitriol (pM)
  ## V1 - vascular volume in litre
  
  
  out$PTHpM <- out$PTH / hpacamod$param$V1
  out$PTHconc <- out$PTHpM * 9.4
  out$BSAP <- with(out, OBfast + OBslow)
  out$sCTx <- out$OC
  
  out$CaConc <- out$P / hpacamod$param$V1
  out$CalcitriolConc <- out$B / hpacamod$param$V1
  out
  
}
