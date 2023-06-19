load_steady_states <- function() {
  tryCatch(
    {
      # Create steady state names based on EST level
      states = c()
      for (EST in c("baseline", "OC_EE", "postmeno_E2")) {
        name = paste("base", EST , sep= "_")
        do.call("<-", list(name, readRDS(paste0("simulations/", name, ".RDS"))))
        states = append(states, name)
      }
      
      for (state in states) {
        df_tail <- tail(get(state), n=1)
        # skip time, osteoblast (OB) and EST level character column
        ss <- df_tail[,2:(length(df_tail)-2)]
        name <- paste(state, "ss", sep="_")
        assign(x = name, value = ss, inherits = TRUE)
        saveRDS(name, file = paste0("steady_states/", paste(state, "ss", sep="_"), ".RDS"))
      }
      message("Steady states were successfully extracted from most recent base model")
    },
    error=function(e) {
      message(e)
      message("Steady states could not be extracted from most recent base model.")
      message("Probably because most recent base model does not exist")
    },
    warning=function(w) {
      message(w)
      message("Steady states could not be extracted from most recent base model.")
      message("Probably because most recent base model does not exist")
    }
  )
}
    