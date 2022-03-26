###
###
# 
# tmp_dag <- matrix(c(rep(0,6),
#                     1, rep(0,5),
#                     1,1,0,0,0,0,
#                     0,1,1,0,0,0,
#                     1,1,1,1,0,0,
#                     1,1,1,1,1,0), nrow =6, byrow = T)
# tmp_dag
# tmp_scm <- sim_SCM(dag = tmp_dag)
# 
# tmp_scm$dag
# 
# W <- c(1,2,3,5,6)


is_simple_lat_proj <- function(SCM, W) {
  mat_dag <- SCM$dag
  size <- NROW(SCM$dag)
  if (length(W) == size) {
    return(TRUE)
  }
  L <- (1:size)[-W]
  
  X <- mat_dag
  X[, W] <- 0 
  X <- find_des_mat(X)
  for (i in seq_along(L)) {
    if (sum(X[W, L[i]]) > 1) {
      return(FALSE)
    }
  }
  return(TRUE)
}

check_prop1 <- function(SCM, MCA, y, x) {
  des_mat <- find_des_mat(SCM$dag)
  
  any_an <- any(des_mat[y, x] == 1)
  if (any_an) return(c(TRUE, FALSE))
  tmp <- MCA(SCM, x, y)
  return(tmp == 0)
}


# Function to check if P2 hold.
check_prop2 <- function(SCM, MCA, y, x) {
  tmp <- MCA(SCM, x, y)
  return(tmp <= 1)
}


# Function to check if P3 hold.
check_prop3 <- function(SCM, MCA, y, x) {
  tmp <- MCA(SCM, x, y)
  return(tmp >= 0)
}


# Function to check if P4 hold.
check_prop4 <- function(SCM, MCA, y, x1, x2, epsi = 0.05) {
  for (i in x1) {
    for (j in x2) {
      if (j == i) {
        #        cat("input not disjoint. \n")
        return(c(TRUE, F))
      }
    }
  }
  tmp <- MCA(SCM, x1, y) + MCA(SCM, x2, y) - MCA(SCM, c(x1, x2), y) 
  tmp2 <- abs(tmp) < epsi
  
  return(abs(tmp) < epsi)
}


# Function to check if P5 hold.
check_prop5 <- function(SCM, MCA, y, x1, x2) {
  x1x2 <- union(x1, x2) 
  tmp <- MCA(SCM, x1, y) + MCA(SCM, x2, y) - MCA(SCM, x1x2, y)
  return(tmp >= 0)
}


# Function to check if P6 hold.
check_prop6 <- function(SCM, MCA, y, x1, x2, epsi= 0.01) {
  for(i in x1) {
    for (j in x2) {
      if (i == j) {
        #       cat("Input not disjoint. \n")
        return(c(TRUE, F))
      }
    }
  }
  
  is_indep <- condi_indep(SCM, x = x2, y = y, w = x1)
  if (!is_indep) {
    return(c(TRUE, F))
  }
  tmp <- MCA(SCM, x1, y) - MCA(SCM, c(x1,x2), y)
  return(abs(tmp) < epsi)
}


# Function to check if P7 Hold.
check_prop7 <- function(SCM, MCA, y, x, W, epsi= 0.05) {
  mat_coef  <- SCM$coef
  mat_coef  <- SCM$dag
  my_noises <- SCM$noises
  size     <- length(my_noises)
  
  L <- (1:size)[-W]
  if (length(intersect(c(y, x), L) > 0)) {
    #    cat("W does not contain y and x. \n")
    return(c(TRUE, F))
  }
  
  if (!is_simple_lat_proj(SCM, W)) {
    #    cat("Not simple lat proj. \n")
    return(c(TRUE, F))
  }
  

  SCM2 <- find_lat_proj(SCM, W)
  
  a1 <- MCA(SCM, x, y)
  a2 <- MCA(SCM2, x, y)
  tmp <- a1 - a2
  return(abs(tmp) < epsi)
}

all_check_factory <- function(SCM) {
  list_of_prop <- c(check_prop1, check_prop2, check_prop3, check_prop4, 
                    check_prop5, check_prop6, check_prop7)
  list_of_MCA  <- c(measure_A_max, measure_A_rel, measure_A_con, measure_A_min, measure_A_adj)
  
  f_result <- function(property = 1, MCA = 1, y, x1, x2 = NULL, W = NULL) {
    check_fun <- list_of_prop[[property]]
    MCA       <- list_of_MCA[[MCA]]
    
    if (property == 1) {
      result <- check_fun(SCM, MCA, y, x1)
      return (result)
    }
    
    if (property == 2) {
      result <- check_fun(SCM, MCA, y, x1)
      return(result)
    }
    
    if (property == 3) {
      result <- check_fun(SCM, MCA, y, x1)
      return(result)
    }
    
    if (property == 4) {
      result <- check_fun(SCM, MCA, y, x1, x2)
      return(result)
    }
    
    if (property == 5) {
      result <- check_fun(SCM, MCA, y, x1, x2)
      return(result)
    }
    
    if (property == 6) {
      result <- check_fun(SCM, MCA, y, x1, x2)
      return(result)
    }
    
    if (property == 7) {
      result <- check_fun(SCM, MCA, y, x1, W)
      return(result)
    }
  }
  return(f_result)
}

