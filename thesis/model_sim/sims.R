#
#
#
#
# This file contains functions to simulate DAGS and SCM over those DAGS
#
#
#
#
#


# Function to simulate a DAG.
# Input: size:   number of vertices
#        p_edge: probability of creating an edge from i to j for all i > j
# Outout: a matrix with 1's and 0's indicating directed edge from column to row
sim_DAG <- function(size, p_edge = 0.7) {
  mat_dag <- matrix(0, nrow = size, ncol = size)
  for (i in 2:size) {
    for (j in 1:(i - 1)) {
      mat_dag[i,j] <- sample(c(1, 0), 1, prob = c(p_edge, 1 - p_edge))
    }
  }
  mat_dag
}

# Function to generate a linear SCM.
# Input dag:   DAG to create the model. If NULL then it is randomly generated
#       theta:   coefficient of the SCM are normal distribution with mean theta
#       sigmasq: -||- with variance sigmasq
#       size: If DAG is NULL, then size picks the number of variables.
# Output: see function.

sim_SCM <- function(dag = NULL,
                    theta = 1,
                    sigmasq = 1,
                    size = NULL) {
  if (is.null(dag)) {
    if (is.null(size)) {
      dag <- sim_DAG(6)
    } else {
      dag <- sim_DAG(size)
    }
  }
  size <- NROW(dag)
  mat_coef <- dag
  for (i in 2:size) {
    for (j in 1:(i - 1)) {
      if (mat_coef[i,j] == 1) {
        mat_coef[i,j] <- sample(c(-1, 1), 1) * rnorm(1, mean = theta, sd = sqrt(sigmasq))
      }
    }
  }
  noises <- runif(size, min = 0.5, max = 1.5)
  return(list(
    dag    = dag,
    coef   = mat_coef,
    noises = noises
  ))
}
