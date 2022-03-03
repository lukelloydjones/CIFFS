# Reading optimiser data produced by the HPC

setwd("/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_ciff_design/ckmr_design_hpc_out_var_props")

#optim_arr <- array(dim = c(5, 3, 90), dimnames = list(c(2014:2018), c("juv", "subad", "ad"), NULL))
optim_list <- list()
cv_arr <- array(dim = c(90, 3), dimnames = list(NULL, c("max_samps_tot", "rep", "min_cv")))

par <- read.table(paste0("design_300_1_3_par.txt"))
val <- read.table(paste0("design_300_1_3_val.txt"))

row_no = 1
for (max_samp_no in c(250, (3:10 * 100))) {
  for (rep in 1:10) {
    par <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_par.txt"))
    dimnames(par) <- list(c(2014:2018), c("juv", "subad", "ad"))
    val <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_val.txt"))
    
    optim_list[[row_no]] <- as.matrix(par)
    cv_arr[row_no, ] <- c(max_samp_no, rep, val[ , ])
    row_no = row_no + 1
  }
}

optim_df <- data.frame(cv_arr, check.names = F)
optim_df$pars <- optim_list
plot(optim_df$max_samps_tot, optim_df$min_cv)
plot(optim_df$year | optim_df$max_samps_tot, optim_df$p)


plot(optim_df$year | optim_df$max_samps_tot, optim_df$p)

xyplot(juv[, '2014'] ~ max_samps_tot,	optim_df)
xyplot(juv[, '2015'] ~ max_samps_tot,	optim_df)
xyplot(juv[, '2016'] ~ max_samps_tot,	optim_df)
xyplot(juv[, '2017'] ~ max_samps_tot,	optim_df)
xyplot(juv[, '2018'] ~ max_samps_tot,	optim_df)


xyplot(optim_arr[,'juv', ] ~ optim_df$max_samps_tot)
grid()

mean(optim_df$min_cv)
min(optim_df[which(optim_df$max_samps_tot==500), 'min_cv'])

write.csv(optim_df[, 1:3], "hpc_out_var_props_all", row.names = FALSE)
write.csv(aperm(optim_arr, c(3, 1, 2)), "hpc_out_var_props_all_mats")
optim_df$juv <- aperm(optim_arr, c(3, 1, 2))[, , 'juv']
optim_df$subad <- aperm(optim_arr, c(3, 1, 2))[, , 'subad']
optim_df$ad <- aperm(optim_arr, c(3, 1, 2))[, , 'ad']

write.csv(optim_df$juv, "hpc_out_var_props_juv", row.names = FALSE)
write.csv(optim_df$subad, "hpc_out_var_props_subad", row.names = FALSE)
write.csv(optim_df$ad, "hpc_out_var_props_ad", row.names = FALSE)
