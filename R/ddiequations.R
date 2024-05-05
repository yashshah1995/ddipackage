# For the sum you would provide values in a vector

reversible_cyp_inhibition <- function(f_guti,
                                      f_gut,
                                      f_m,
                                      f_m_CYP,
                                      I_u,
                                      K_iu) {

  # Put all vectors in a list
  vectors <- list(f_guti, f_gut, f_m, f_m_CYP, I_u, K_iu)

  # Check if all vectors have the same length
  if (!all(sapply(vectors, length) == length(vectors[[1]]))) {
    stop("All vectors must have the same length.")
  }

  # Initialize the total sum for vector additions
  total_sum_firstpart <- 0
  total_sum_secondpart <- 0

  # Loop over the length of the vectors
  for (i in 1:length(f_m)) {
    # Calculate the fraction for the current set of values
    fm_fmcyp <- f_m[i] * f_m_CYP[i]
    denom <- 1 + (I_u[i]/K_iu[i])
    fraction <- fm_fmcyp / denom

    # Add the fraction to the total sum
    total_sum_firstpart <- total_sum_firstpart + fraction

    total_sum_secondpart <- total_sum_secondpart + fm_fmcyp
  }

  # Calculate the reversible CYP inhibition equation
  #(AUC_i / AUC)
  result <- (f_guti / f_guti) * (1 / (total_sum_firstpart + (1- total_sum_secondpart)))

  # Return the result
  return(result)
}

tdi <- function(f_guti,
                f_gut,
                f_m,
                f_m_CYP,
                I_u,
                K_iu,
                k_inact,
                k_deg
                ) {

  # Put all vectors in a list
  vectors <- list(f_guti, f_gut, f_m, f_m_CYP, I_u, K_iu, k_inact, k_deg)

  # Check if all vectors have the same length
  if (!all(sapply(vectors, length) == length(vectors[[1]]))) {
    stop("All vectors must have the same length.")
  }

  # Initialize the total sum for vector additions
  total_sum_firstpart <- 0
  total_sum_secondpart <- 0

  # Loop over the length of the vectors
  for (i in 1:length(f_m)) {
    # Calculate the fraction for the current set of values
    fm_fmcyp <- f_m[i] * f_m_CYP[i]
    denom <- 1 + ((k_inact[i] * I_u[i]) / (k_deg[i] *
                                             (K_iu[i] + I_u[i])
                                           )
                  )
    fraction <- fm_fmcyp / denom

    # Add the fraction to the total sum
    total_sum_firstpart <- total_sum_firstpart + fraction

    total_sum_secondpart <- total_sum_secondpart + fm_fmcyp
  }

  # Calculate the TDI equation
  result <- (f_guti / f_guti) * (1 / (total_sum_firstpart + (1- total_sum_secondpart)))

  # Return the result
  return(result)
}

induction <- function(f_guti,
                      f_gut,
                      f_m,
                      f_m_CYP,
                      I_u,
                      Emax_ind,
                      EC50_ind) {


  # Put all vectors in a list
  vectors <- list(f_guti, f_gut, f_m, f_m_CYP, I_u, Emax_ind, EC50_ind)

  # Check if all vectors have the same length
  if (!all(sapply(vectors, length) == length(vectors[[1]]))) {
    stop("All vectors must have the same length.")
  }

  #TODO(): The induction equation is confusing...
  # Which summation part is required?
  # 1. Summation of entire values before the + sign, remove all other summations
  # 2. Summation of entire values before the + sign, keep the sum of f_m
  # 3. Summation of the values after second * sign, keep the sum of f_m

  # Initialize the total sum for vector additions
  total_sum_firstpart <- 0
  total_sum_secondpart <- 0

  # Loop over the length of the vectors
  for (i in 1:length(f_m)) {
    # Calculate the fraction for the current set of values
    fm_fmcyp <- f_m[i] * f_m_CYP[i]
    denom <- 1 + ((k_inact[i] * I_u[i]) / (k_deg[i] *
                                             (K_iu[i] + I_u[i])
    )
    )
    fraction <- fm_fmcyp / denom

    # Add the fraction to the total sum
    total_sum_firstpart <- total_sum_firstpart + fraction

    total_sum_secondpart <- total_sum_secondpart + fm_fmcyp
  }



  # # Calculate the summation terms
  # sum_f_m_CYP <- sum(f_m_CYP)
  # sum_L_u_K_iu_Emax_ind_EC50_ind <- sum(L_u / (K_iu + Emax_ind / (EC50_ind + Emax_ind)))
  #
  # # Calculate the Induction equation
  # result <- (AUC_i / AUC) * (f_gi / f_gui) * (1 / (1 + sum_L_u_K_iu_Emax_ind_EC50_ind)) * (sum_f_m_CYP + (1 - sum_f_m_CYP))

  # Return the result
  return(result)
}


