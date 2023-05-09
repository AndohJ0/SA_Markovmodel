## Model Function  Zimmerman

## Base case model ----
# parameters
parms <- list(
  p_N_LG = 0.001,
  p_N_D = 0.00023,
  p_LG_D = 0.00023,
  p_HG_D = 0.00023,
  p_LG_HG = 0.003,
  p_HG_LG = 0.002,
  p_HG_N = 0.001,
  p_HG_LIC = 0.002,
  p_LIC_RIC = 0.003,
  p_RIC_DIC = 0.003,
  p_LG_N = 0.003,
  p_LIC_D = 0.004,
  p_RIC_D = 0.021,
  p_DIC_D = 0.063
)

model1 <- function(.parms) {
  with(.parms, {
    #fixed
    
    n_s <- 7   # number of states
    n_t <- 840 # number of cycles in months (10 years)
    n_c <- 1e3 # cohort size
    
    st_names <- c("Normal", "LGrade", "HGrade", "LIC", "RIC", "DIC", "Death")
    
    ## Probability Transition Matrix
    # initial
    
    m_P <- matrix(0,
                  nrow = n_s,
                  ncol = n_s,
                  byrow = T,
                  dimnames = list(from = st_names,
                                  to = st_names))
    
    # prob. equations
    m_P["Normal", "LGrade"]  <- p_N_LG
    m_P["Normal", "Death"]   <- p_N_D
    m_P["LGrade", "Normal"]  <- p_LG_N
    m_P["LGrade", "HGrade"]  <- p_LG_HG
    m_P["LGrade", "Death"]   <- p_LG_D
    m_P["HGrade", "Normal"]  <- p_HG_N
    m_P["HGrade", "LGrade"]  <- p_HG_LG
    m_P["HGrade", "Death"]   <- p_HG_D
    m_P["HGrade", "LIC"]     <- p_HG_LIC
    m_P["LIC", "RIC"]        <- p_LIC_RIC
    m_P["RIC", "DIC"]        <- p_RIC_DIC
    m_P["LIC", "Death"]      <- p_LIC_D
    m_P["RIC", "Death"]      <- p_RIC_D
    m_P["DIC", "Death"]      <- p_DIC_D
    m_P["Death", "Death"]    <- 1
    
    m_P["Normal", "Normal"] <- 1 - sum(m_P["Normal",-1])
    m_P["LGrade", "LGrade"] <- 1 - sum(m_P["LGrade", -2])
    m_P["HGrade", "HGrade"] <- 1 - sum(m_P["HGrade", -3])
    m_P["LIC", "LIC"]       <- 1 - sum(m_P["LIC", -4])
    m_P["RIC", "RIC"]       <- 1 - sum(m_P["RIC", -5])
    m_P["DIC", "DIC"]       <- 1 - sum(m_P["DIC", -6]) 
    
    st_mem <- array(NA,
                    dim = c(n_t, n_s),
                    dimnames = list(cycle = 1:n_t,
                                    state = st_names))
    st_mem[1,] <- c(n_c, rep(0, n_s - 1))
    
    for (i in 2:n_t) {
      st_mem[i,] <- st_mem[i - 1,] %*% m_P
    }
    
    sum_states <- apply(st_mem, 2, sum)
    state_exp <- sum_states[1:6] / n_c
    life_exp <- sum(state_exp)
    list(prob.matrix = m_P,
         life_exp = life_exp)
  })
}


## Model with interventions ----

model2 <- function(.parms, sens,spec) { #measuring the impact of screening
  with(.parms, {
    #fixed
    
    n_s <- 7   # number of states
    n_t <- 840 # number of cycles in months (10 years)
    n_c <- 1e3 # cohort size
    
    st_names <- c("Normal", "LGrade", "HGrade", "LIC", "RIC", "DIC", "Death")
    
    # Test characteristics
    lr <- sens/(1-spec)         # likelihood ratio of test
    RR <- lr/(1+lr)
    
    ## Probability Transition Matrix
    # initial
    
    m_P <- matrix(0,
                  nrow = n_s,
                  ncol = n_s,
                  byrow = T,
                  dimnames = list(from = st_names,
                                  to = st_names))
    
    # prob. equations
    m_P["Normal", "LGrade"]  <- p_N_LG * (1-RR)
    m_P["Normal", "Death"]   <- p_N_D
    m_P["LGrade", "Normal"]  <- p_LG_N * (1+RR)
    m_P["LGrade", "HGrade"]  <- p_LG_HG * (1-RR)
    m_P["LGrade", "Death"]   <- p_LG_D 
    m_P["HGrade", "Normal"]  <- p_HG_N * (1+RR)
    m_P["HGrade", "LGrade"]  <- p_HG_LG * (1+RR)
    m_P["HGrade", "Death"]   <- p_HG_D 
    m_P["HGrade", "LIC"]     <- p_HG_LIC * (1-RR)
    m_P["LIC", "RIC"]        <- p_LIC_RIC 
    m_P["RIC", "DIC"]        <- p_RIC_DIC
    m_P["LIC", "Death"]      <- p_LIC_D 
    m_P["RIC", "Death"]      <- p_RIC_D 
    m_P["DIC", "Death"]      <- p_DIC_D 
    m_P["Death", "Death"]    <- 1
    
    m_P["Normal", "Normal"] <- 1 - sum(m_P["Normal",-1])
    m_P["LGrade", "LGrade"] <- 1 - sum(m_P["LGrade", -2])
    m_P["HGrade", "HGrade"] <- 1 - sum(m_P["HGrade", -3])
    m_P["LIC", "LIC"]       <- 1 - sum(m_P["LIC", -4])
    m_P["RIC", "RIC"]       <- 1 - sum(m_P["RIC", -5])
    m_P["DIC", "DIC"]       <- 1 - sum(m_P["DIC", -6]) 
    
    st_mem <- array(NA,
                    dim = c(n_t, n_s),
                    dimnames = list(cycle = 1:n_t,
                                    state = st_names))
    st_mem[1,] <- c(n_c, rep(0, n_s - 1))
    
    for (i in 2:n_t) {
      st_mem[i,] <- st_mem[i - 1,] %*% m_P
    }
    
    sum_states <- apply(st_mem, 2, sum)
    state_exp <- sum_states[1:6] / n_c
    life_exp <- sum(state_exp)
    list(prob.matrix = m_P,
         life_exp = life_exp)
  })
}

## Time varying probabilities ----
parms <- list(
  p_N_LG = 0.001,
  p_N_D = 0.00023,
  p_LG_D = 0.00023,
  p_HG_D = 0.00023,
  p_LG_HG = 0.003,
  p_HG_LG = 0.002,
  p_HG_N = 0.001,
  p_HG_LIC = 0.002,
  p_LIC_RIC = 0.003,
  p_RIC_DIC = 0.003,
  p_LG_N = 0.003,
  p_LIC_D = 0.004,
  p_RIC_D = 0.021,
  p_DIC_D = 0.063
)
model3 <- function(.parms) { 
  with(.parms, {
    #fixed
    
    n_s <- 7   # number of states
    n_t <- 840 # number of cycles in months (10 years)
    n_c <- 1e3 # cohort size
    
    st_names <- c("Normal", "LGrade", "HGrade", "LIC", "RIC", "DIC", "Death")
    
    ## Probability Transition "Array"
    # initial
    
    tvm_P <- array(NA,
                   dim = c(n_s, n_s, n_t),
                   dimnames = list(from  = st_names, 
                                   to    = st_names,
                                   cycle = 1:n_t)) # NA array, dim = dim of mat x cycle
    
    
    #Some initial prob. conditions
    tvm_P["Death",1:6 , 1:n_t]     <- 0
    tvm_P["Death","Death", 1:n_t]  <- 1
    
    tvm_P["DIC", 1:5, ]            <- 0
    tvm_P["DIC", "Death", 1:456]   <- p_DIC_D
    tvm_P["DIC", "Death", 457:576] <- p_DIC_D + 0.002
    tvm_P["DIC", "Death", 577:840] <- p_DIC_D + 0.005
    tvm_P["DIC", "DIC", ]          <- 1 - apply(tvm_P["DIC", , ],2, sum, na.rm = T)
    
    tvm_P["RIC", 1:4, ]            <- 0
    tvm_P["RIC", "DIC", 1:456]     <- p_RIC_DIC
    tvm_P["RIC", "DIC", 457:576]   <- p_RIC_DIC + 0.002
    tvm_P["RIC", "DIC", 577:840]   <- p_RIC_DIC + 0.005
    tvm_P["RIC", "Death", 1:456]   <- p_RIC_D
    tvm_P["RIC", "Death", 457:576] <- p_RIC_D + 0.002
    tvm_P["RIC", "Death", 577:840] <- p_RIC_D + 0.005
    tvm_P["RIC", "RIC", ]          <- 1 - apply(tvm_P["RIC", , ],2, sum, na.rm = T)
    
    tvm_P["LIC", 1:3, ]            <- 0
    tvm_P["LIC", 6, ]              <- 0
    tvm_P["LIC", "RIC", 1:456]     <- p_LIC_RIC
    tvm_P["LIC", "RIC", 457:576]   <- p_LIC_RIC + 0.002
    tvm_P["LIC", "RIC", 577:840]   <- p_LIC_RIC + 0.005
    tvm_P["LIC", "Death", 1:456]   <- p_LIC_D
    tvm_P["LIC", "Death", 457:576] <- p_LIC_D + 0.002
    tvm_P["LIC", "Death", 577:840] <- p_LIC_D + 0.005
    tvm_P["LIC", "LIC", ]          <- 1 - apply(tvm_P["LIC", , ],2, sum, na.rm = T)
    
    tvm_P["HGrade", 5:6 , ]             <- 0
    tvm_P["HGrade", "LIC", 1:456]       <- p_HG_LIC
    tvm_P["HGrade", "LIC", 457:576]     <- p_HG_LIC + 0.002
    tvm_P["HGrade", "LIC", 577:840]     <- p_HG_LIC + 0.005
    tvm_P["HGrade", "LGrade", 1:456]     <- p_HG_LG
    tvm_P["HGrade", "LGrade", 457:576]   <- p_HG_LG + 0.002
    tvm_P["HGrade", "LGrade", 577:840]   <- p_HG_LG + 0.005
    tvm_P["HGrade", "Normal", 1:456]     <- p_HG_N
    tvm_P["HGrade", "Normal", 457:576]   <- p_HG_N + 0.002
    tvm_P["HGrade", "Normal", 577:840]   <- p_HG_N + 0.005
    tvm_P["HGrade", "Death", 1:456]   <- p_HG_D
    tvm_P["HGrade", "Death", 457:576] <- p_HG_D + 0.002
    tvm_P["HGrade", "Death", 577:840] <- p_HG_D + 0.005
    tvm_P["HGrade", "HGrade", ]          <- 1 - apply(tvm_P["HGrade", , ],2, sum, na.rm = T)
    
    tvm_P["LGrade", 4:6 , ]             <- 0
    tvm_P["LGrade", "HGrade", 1:456]    <- p_LG_HG
    tvm_P["LGrade", "HGrade", 457:576]  <- p_LG_HG + 0.002
    tvm_P["LGrade", "HGrade", 577:840]  <- p_LG_HG + 0.005
    tvm_P["LGrade", "Normal", 1:456]    <- p_LG_N
    tvm_P["LGrade", "Normal", 457:576]  <- p_LG_N + 0.002
    tvm_P["LGrade", "Normal", 577:840]  <- p_LG_N + 0.005
    tvm_P["LGrade", "Death", 1:456]     <- p_LG_D
    tvm_P["LGrade", "Death", 457:576]   <- p_LG_D + 0.002
    tvm_P["LGrade", "Death", 577:840]   <- p_LG_D + 0.005
    tvm_P["LGrade", "LGrade", ]         <- 1 - apply(tvm_P["LGrade", , ],2, sum, na.rm = T)
    
    tvm_P["Normal", 3:6 , ]             <- 0
    tvm_P["Normal", "LGrade", 1:456]    <- p_N_LG
    tvm_P["Normal", "LGrade", 457:576]  <- p_N_LG + 0.002
    tvm_P["Normal", "LGrade", 577:840]  <- p_N_LG + 0.005
    tvm_P["Normal", "Death", 1:456]     <- p_N_D
    tvm_P["Normal", "Death", 457:576]   <- p_N_D + 0.002
    tvm_P["Normal", "Death", 577:840]   <- p_N_D + 0.005
    tvm_P["Normal", "Normal", ]         <- 1 - apply(tvm_P["Normal", , ],2, sum, na.rm = T)
    
    st_mem <- array(NA, dim = c(n_t, n_s), 
                    dimnames = list(cycle = 1:n_t, state = st_names)) #Empty array to initialize
    
    st_mem[1,] <- c(n_c, rep(0, n_s - 1))
    
    for (i in 2:n_t){  
      st_mem[i, ] <- st_mem[i-1, ] %*% tvm_P[ , , i-1]
    }  
    
    sum_states <- apply(st_mem, 2, sum)
    state_exp <- sum_states[1:6] / n_c
    life_exp <- sum(state_exp)
    list(life_exp = life_exp)
  })
}

