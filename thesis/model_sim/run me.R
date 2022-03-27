#
#
#
#
# This is the file to run for producing the result in chapter 4 of my thesis.
# First load all other files below
# then the four simulations below will produce the exact result used in the thesis.
#
#
#
load <- TRUE
if (load) {
  source("model_sim/sims.R")
  source("model_sim/measures.R")
  source("model_sim/aux_fun.R")
  source("model_sim/calculations.R")
  source("model_sim/check_fun.R")
  source("model_sim/check_all.R")
  load <- FALSE
}


times <- 1000

# Simulation A. See thesis
sim_big_more(times = times, size = 3, seed = 1)

# Simulation B.
sim_big_more(times = times, size = 4, seed = 1)

# Simulation C.
sim_big_more(times = times, size = 5, seed = 1)

#Simulation D.
sim_big_more(times = times, size = 6, seed = 1)


