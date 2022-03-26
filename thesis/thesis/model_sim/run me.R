
load <- TRUE
if (load) {
  source("model_sim/sims.R")
  source("model_sim/measures.R")
  source("model_sim/aux_fun.R")
  source("model_sim/calculations.R")
  source("model_sim/check_fun.R")
  source("model_sim/check_all.R")
  load = FALSE
}
## 

times <- 1000

sim_big_more(times = times, size = 3, seed = 1)

sim_big_more(times = times, size = 4, seed = 1)

sim_big_more(times = times, size = 5, seed = 1)

sim_big_more(times = times, size = 6, seed = 1)


