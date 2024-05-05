# For the sum you would provide values in a vector

reversible_cyp_inhibition <- function(AUC_i, AUC, f_gi, f_gui, f_m_CYP, L_u, K_iu) {
  # Calculate the summation terms
  sum_f_m_CYP <- sum(f_m_CYP)
  sum_L_u_K_iu <- sum(L_u / K_iu)

  # Calculate the reversible CYP inhibition equation
  #(AUC_i / AUC)
  result <- (f_gi / f_gui) * (1 / (1 + sum_L_u_K_iu)) * (sum_f_m_CYP + (1 - sum_f_m_CYP))

  # Return the result
  return(result)
}

tdi <- function(AUC_i, AUC, f_gi, f_gui, f_m_CYP, L_u, K_iu, kinact, KI) {
  # Calculate the summation terms
  sum_f_m_CYP <- sum(f_m_CYP)
  sum_L_u_K_iu_kinact_KI <- sum(L_u / (K_iu + kinact / KI))

  # Calculate the TDI equation
  result <- (AUC_i / AUC) * (f_gi / f_gui) * (1 / (1 + sum_L_u_K_iu_kinact_KI)) * (sum_f_m_CYP + (1 - sum_f_m_CYP))

  # Return the result
  return(result)
}

induction <- function(AUC_i, AUC, f_gi, f_gui, f_m_CYP, L_u, K_iu, Emax_ind, EC50_ind) {
  # Calculate the summation terms
  sum_f_m_CYP <- sum(f_m_CYP)
  sum_L_u_K_iu_Emax_ind_EC50_ind <- sum(L_u / (K_iu + Emax_ind / (EC50_ind + Emax_ind)))

  # Calculate the Induction equation
  result <- (AUC_i / AUC) * (f_gi / f_gui) * (1 / (1 + sum_L_u_K_iu_Emax_ind_EC50_ind)) * (sum_f_m_CYP + (1 - sum_f_m_CYP))

  # Return the result
  return(result)
}


