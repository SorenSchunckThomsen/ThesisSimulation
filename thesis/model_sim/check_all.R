#
#
#
#
# This file contains function check properties. See thesis
#
#
#
#
#

# Function to return subsets of x
# Returns a function, which for each call returns next subset of x
subset_factory <- function(x) {
  val <- rep(0, length(x))
  f <- function() {
    val <<- add_one(val)
    if (is.null(val)) {
      return(NULL)
    }
    return(x[val == 1])
  }
  return(f)
}

# Function to help subset_factory
add_one <- function(x) {
  if (all(x == 1)) {
    return(NULL)
  }
  if (x[1] == 0){
    x[1] <- 1
    return(x)
  } else {
    x <- c(0, add_one(x[-1]))
    return(x)
  }
}

# Function to simulate models and check properties
# Input: size of SCM's to be simulated
#        seed for reproducibility
# Outputs: 5x7 table with 1's and 0's. Entry ij incodes if measure i satisfied property j.
sim_big_all <- function(size = 6, seed = NULL, silent = T) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  result <- matrix(0, nrow = 5, ncol = 7)
  
  my_SCM <- sim_SCM(size = size)
  
  my_dag <- my_SCM$dag
  my_coef <- my_SCM$coef
  my_noises <- my_SCM$noises
  
  size <- length(my_noises)
  
  my_fun <- all_check_factory(my_SCM)
  
  # my_fun has 
  # Input: (property = 1, MCA = 1, y, x1, x2 = NULL, W = NULL)
  # Output: TRUE/FALSE
  
  # target = pick_Y
  pick_Y <- sample(c(-1,0) + size, 1)
  allx <- (1:size)[-pick_Y]
  
  # Checking P1
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      satis <- my_fun(property = 1, MCA = i, pick_Y, x1 = sub_set)[1]
      sub_set  <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 1] <- 1
        satis <- F
      }
    }
  }
  if (!silent) print("P1 done")
  
  # Checking P2
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      satis <- my_fun(property = 2, MCA = i, pick_Y, x1 = sub_set)[1]
      sub_set  <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 2] <- 1
        satis <- F
      }
    }
  }
  if (!silent) print("P2 done")
  
  # Checking P3
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      satis <- my_fun(property = 3, MCA = i, pick_Y, x1 = sub_set)[1]
      sub_set  <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 3] <- 1
        satis <- F
      }
    }
  }
  if (!silent) print("P3 done")
  
  # Checking P4
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      no_x1_no_y   <- (1:size)[-c(pick_Y, sub_set)] 
      pick_subset2 <- subset_factory(no_x1_no_y)
      sub_set2 <- pick_subset2()
      while(satis & !is.null(sub_set2)) {
        satis <- my_fun(property = 4, MCA = i, pick_Y, x1 = sub_set, x2 = sub_set2)[1]
        sub_set2  <- pick_subset2()
      }
      sub_set <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 4] <- 1
        satis <- F
      }
      
    }
  }
  if (!silent) print("P4 done")
  
  # Checking P5
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      pick_subset2 <- subset_factory(allx)
      sub_set2 <- pick_subset2()
      while(satis & !is.null(sub_set2)) {
        satis <- my_fun(property = 5, MCA = i, pick_Y, x1 = sub_set, x2 = sub_set2)[1]
        sub_set2  <- pick_subset2()
      }
      sub_set <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 5] <- 1
        satis <- F
      }
      
    }
  }
  if (!silent) print("P5 done")
  
  # Checking P6
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    
    while(satis) {
      no_x1_no_y   <- (1:size)[-c(pick_Y, sub_set)] 
      pick_subset2 <- subset_factory(no_x1_no_y)
      sub_set2 <- pick_subset2()
      while(satis & !is.null(sub_set2)) {
        satis <- my_fun(property = 6, MCA = i, pick_Y, x1 = sub_set, x2 = sub_set2)[1]
        sub_set2  <- pick_subset2()
 
      }
      sub_set <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 6] <- 1
        satis <- F
      }
      
    }
    
  }
  if (!silent) print("P6 done")
  
  # Checking P7
  for (i in 1:5) {
    satis <- TRUE
    
    pick_subset <- subset_factory(allx)
    sub_set  <- pick_subset()
    sub_set2 <- c()
    
    while(satis) {
      no_y_no_x <- (1:size)[-c(sub_set, pick_Y)] #her
      pick_subset2 <- subset_factory(no_y_no_x)
      set_W <- c(pick_Y, sub_set, sub_set2)
      satis <- my_fun(property = 7, MCA = i, pick_Y, x1 = sub_set, W = set_W)[1]
      sub_set2 <- pick_subset2()
      while(satis & !is.null(sub_set2)) {
        set_W <- c(pick_Y, sub_set, sub_set2)
        satis <- my_fun(property = 7, MCA = i, pick_Y, x1 = sub_set, W = set_W)[1]
        sub_set2 <- pick_subset2()
      }
      sub_set <- pick_subset()
      if (satis & is.null(sub_set)) {
        result[i, 7] <- 1
        satis <- F
      }
    }
  }
  if (!silent) print("P7 done")
  return(result)
}

# Function to use big_sim_all to simulate more than once.
# Input: times: number of SCM to simulate
#        size: number of variables in the SCM's
#        seed: set.seed for reproducibility
sim_big_more <- function(times = 100, size = NULL, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(size)) {
    size <- 6
  }
  
  s <- 0
  for (i in 1:times) {
    s <- s + sim_big_all(size = size)
  }
  return(s / times)
}
