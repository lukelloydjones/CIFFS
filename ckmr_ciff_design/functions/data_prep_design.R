# # Start set for testing
# samp.no.py  = rep(c(50, 90), 5) 
# years.samp = 10 
# prop_ju_py = 0.3  
# prop_sa_py = 0.2
# prop_ad_py = 0.5
# lethal = 0
# A - 84,  J - 28,  SA - 18 
data_prep_design <- function(samp.no.py, years.samp, prop_ju_py, prop_sa_py, prop_ad_py, lethal)
{
  # PARS
  #  samp.no.py  - number of individuals to sample from each year
  #  years.samp - number of years to sample
  #  prop_ju_py - proportion juveniles per year
  #  prop_sa_py - "" sub-adults ""
  #  prop_ad_py - "" adults ""
  #  lethal - 1 => Lethal sampling, 0 => mark-recapture
  # RETURNS
  #  saves the environment with the combinations in each Bju,Yad for CKMR design hessian and grad computations
  
  # ---------------------------------------
  # 
  # ---------------------------------------
  source("functions/prep_comb_kin.R")
  prob_vec <- c(prop_ju_py, prop_sa_py, prop_ad_py)
  #print(paste0("Preping kin and combinations for population sample ", repi))
  
  # ----------------------
  # Making of combinations
  # ----------------------
  # Read keystone sample
  pop_samp  <- read.table(paste0("keystone_sample/pop_samp_", 10, "_yrs_samp_", 
                            10, "_samp_no_py_", 150, ".txt"), 
                            header = T, stringsAsFactors = F) #### Careful here! 10 and 150 are hard coded. You're welcome future Luke.
  
  pop_samp <- pop_samp[pop_samp$year_samp > 2018-years.samp , ]
  
  unique(pop_samp$year_samp)
  
  # Partition by pars an array that stores data by year - within year partition by ju, sa, ad
  data_parts     <- ls()
  data_parts     <-  pop_samp %>% group_by(year_samp, life.sg.status)  
  year_acs_props <- data_parts %>% summarise(size = n())
  year_acs_props_ls <- year_acs_props %>% group_by(life.sg.status)
  min_year_acs_props_ls <- year_acs_props_ls %>% summarise(min(size))
  samp_size_upper_bound <- sum(min_year_acs_props_ls[,2])
  
  names(prob_vec) <- c("J", "SA", "A")
  names(samp.no.py) <- unique(pop_samp$year_samp)
  
  for (i in seq(1, dim(min_year_acs_props_ls)[1]))
  {
    ac <- as.character(min_year_acs_props_ls[i, 1])
    ac_num  <- prob_vec[ac] * max(samp.no.py)
    if (ac_num > min_year_acs_props_ls[i, 2]) 
    {
      stop(paste0("Unable to sample these proportions (violated by ",
                  ac, " - ", ac_num,
                  "). Rearrange to obtain samp.no.py less than or equal to ",
                  min_year_acs_props_ls[1, 1], " - ", min_year_acs_props_ls[1, 2], ", ",
                  min_year_acs_props_ls[2, 1], " - ", min_year_acs_props_ls[2, 2], ", ",
                  min_year_acs_props_ls[3, 1], " - ", min_year_acs_props_ls[3, 2], "."))
    }
  }
  #if (!all(prob_vec * max(samp.no.py) <= c(28, 18, 84))) {stop("Unable to sample these proportions. Rearrange to obtain samp.no.py less than or equal to J - 28,  SA - 18, A - 84 ")}
  
  parts_attr_grps <- attributes(data_parts)$groups
  parts_attr_inds <- attributes(data_parts)$groups$.rows
  # Bind up the data frame
  pop_samp_design <- data.frame(matrix(0, ncol = dim(pop_samp)[2]))
  colnames(pop_samp_design) <- colnames(pop_samp)
  for (i in seq(1, dim(parts_attr_grps)[1]))
  {
    #print(parts_attr_grps[i, c(1, 2)])
    attr_grps_i     <- as.character(parts_attr_grps[i, c(1, 2)])
    nos_this_year   <- floor(samp.no.py[attr_grps_i[1]] * prob_vec[attr_grps_i[2]]) # This could be problematic but ignore for now 11/12/2019
    data_parts_df_i <- pop_samp[parts_attr_inds[[i]], ] %>% sample_n(size = nos_this_year)
    pop_samp_design <- rbind(pop_samp_design, data_parts_df_i)
  }
  pop_samp_design <- pop_samp_design[-1, ]
  
  # Return comparison sets
  kin     <- new.env()
  com.kin <- prep_comb_kin(kin, pop_samp_design, lethal)
  return(com.kin)
}
