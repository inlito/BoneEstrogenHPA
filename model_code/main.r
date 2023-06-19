library(deSolve)

# Load model, initial values and parameters
source("model_code/library.r")
hpacamod <- load.model()
hpacamod <- derive.init(hpacamod)

# Load script containing function for creating steady states
source("model_code/steady_states.r")

## Teriparatide dosing for the teriparatide simulation
teri.func <- function() {
  ## TERIPARATIDE DOSING EVENTS (TIMES IN HOURS)
  teri.times <- seq(12, 24, 12*30*24)
  teri.dose.mcg <- 20 #mcg
  teri.dose <- teri.dose.mcg * 1E6 / 4117.8
  
  events <- data.frame(var = "TERISC", time = teri.times, value = teri.dose, method = "add")

  ## ADD DOSING EVENTS TO EVALUATION TIMES
  times <- sort(unique(c(times, events$time)))
  
  events <- list(data=events)
  return(events)
}

## Simulate using lsoda, outputs dataframe
simulation <- function(hpacamod=hpacamod, 
                       stress=stress, 
                       EST_group=EST_group, 
                       renal_state=renal_state, 
                       PTH_state=PTH_state, 
                       teri=NULL, 
                       steady_state=steady_state, 
                       duration=duration) {
  df <- as.data.frame(
    lsoda(
      y = unlist(steady_state),
      times = seq(0, duration, 1),
      func = hpacamod$model,
      parms = hpacamod$param,
      stress = stress,
      EST_group=EST_group, 
      renal_state=renal_state,
      PTH_state=PTH_state,
      rtol = 1e-6, # relative tolerance
      atol = 1E-6, # absolute tolerance
      ynames = F,
      events = teri
    )
  )
  df$OB <- df$OBfast + df$OBslow
  return(df)
}


## Simulation with three levels of estradiol
simulate_state <- function(stress = NULL,
                           renal_state = NULL, 
                           PTH_state = NULL, 
                           teri = FALSE, 
                           prefix = prefix, 
                           steady_state=NULL,
                           duration=12*30*24) {
    if (teri) {
      df <- simulation(hpacamod, 
                       EST_group = "baseline", 
                       teri = teri.func(), 
                       steady_state = get("baseline_ss")
                       )
      df$EST_group <- "baseline"
      saveRDS(df, file = "simulations/teriparatide.RDS")
    }
    else {
      for (EST_group in c("baseline", "postmeno_E2", "OC_EE")) {
        ss_name <- paste("base", EST_group, "ss", sep = "_")
        if (prefix != "base") {
          steady_state <- get(ss_name)
        }
        df <- simulation(hpacamod=hpacamod, 
                         stress=stress, 
                         renal_state=renal_state, 
                         PTH_state=PTH_state,
                         EST_group=EST_group, 
                         steady_state=steady_state, 
                         duration=duration
                         )
        df$EST_group <- EST_group
        saveRDS(df, file = paste0("simulations/", prefix, "_", EST_group, ".RDS"))
        }
    }
}



# Run simulations. Simulation duration must be given in hours. -----------------


# base
run_base <- menu(c("Yes, changes were made to the model", 
                   "No, changes were not made to the model"), 
                 title="Run a new base model?")
if (run_base == 1) {
  simulate_state(prefix = "base", steady_state = hpacamod$init[hpacamod$cmt], duration = 24*30*12*12)
}

# Load steady states into environment
load_steady_states()

# # demonstration of new base model
simulate_state(prefix = "constant")

# stress
simulate_state(stress = list(start=24*100, end=24*190), prefix = "stress")

# secondary hyperparathyroidism immediate with GFR drop:
simulate_state(renal_state = "immediate", prefix = "secondary_hyper_immediate")

# primary hyperparathyroidism with longitudinal PTH increase:
simulate_state(PTH_state = "primary_hyper_long", prefix = "primary_hyper_long")

# primary hypoparathyroidism
simulate_state(PTH_state = "primary_hypo", prefix = "primary_hypo")

# 24 hour test
# simulate_state(stress = list(start=12, end=36), prefix = "stress_day", duration = 24*2)

# # teriparatide
# simulate_state(teri = TRUE, prefix = "teriparatide")

# # stress early follicular
# simulate_state(stress = list(start=6, end=17), prefix = "stress_follicular", duration = 3*28)
# 
# # stress luteal
# simulate_state(stress = list(start=17, end=28), prefix = "stress_luteal", duration = 3*28)






