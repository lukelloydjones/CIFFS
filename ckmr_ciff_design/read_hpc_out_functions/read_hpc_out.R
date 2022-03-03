# Reading optimiser data produced by the HPC

setwd("/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_ciff_design/ckmr_design_hpc_out")

optim_mat <- matrix(nrow=95, ncol=12)
dimnames(optim_mat) <- list(NULL, 
                             c("max_samp_no", "rep", 
                               2014:2018, "avg_samp_no",
                               "prop_juv", "prop_subad", "prop_ad", "min_cv"))

row_no = 1
for (max_samp_no in (2:20 * 10)) {
  for (rep in 1:5) {
    par <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_par.txt"))
    val <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_val.txt"))
    
    optim_mat[row_no, ] <- c(max_samp_no, rep, 
                              par[1:5, ], mean(par[1:5, ]),
                              (1000 - sum(par[6:7, ])), par[6:7, ], val[1, ])
    row_no = row_no + 1
  }
}

optim_df <- data.frame(optim_mat, check.names = F)
optim_df[ , 9:11] <- optim_df[ , 9:11]/1000
plot(optim_df$avg_samp_no, optim_df$min_cv)

write.csv(optim_df, "hpc_out_all", row.names = FALSE)