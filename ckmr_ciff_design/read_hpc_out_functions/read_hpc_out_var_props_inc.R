# Reading optimiser data produced by the HPC

setwd("/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_ciff_design/ckmr_design_hpc_out_var_props/design_increase_stable")
setwd("/Users/pet247/Dropbox/smp_summer_isabel_peters/ckmr_ciff_design/ckmr_design_hpc_out_var_props")

optim_arr <- array(dim = c(5, 3, 90), dimnames = list(c(2014:2018), c("juv", "subad", "ad"), NULL))
optim_list <- list()
cv_arr <- array(dim = c(90, 3), dimnames = list(NULL, c("max_samps_tot", "rep", "min_cv")))

row_no = 1
for (max_samp_no in c(250, (3:10 * 100))) {
  for (rep in 1:10) {
    # par <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_par.txt"))
    par <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_increase_par.txt"))
    # par <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_stab_par.txt"))
    dimnames(par) <- list(c(2014:2018), c("juv", "subad", "ad"))
    # val <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_val.txt"))
    val <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_increase_val.txt"))
    # val <- read.table(paste0("design_", max_samp_no, "_1_", rep, "_stab_val.txt"))
    
    optim_list[[row_no]] <- as.matrix(par)
    optim_arr[ , , row_no] <- as.matrix(par)
    cv_arr[row_no, ] <- c(max_samp_no, rep, val[ , ])
    row_no = row_no + 1
  }
}

optim_df <- data.frame(cv_arr, check.names = F)
optim_df$pars <- optim_list
plot(optim_df$max_samps_tot, optim_df$min_cv)
grid()

optim_df$juv <- aperm(optim_arr, c(3, 1, 2))[, , 'juv']
optim_df$subad <- aperm(optim_arr, c(3, 1, 2))[, , 'subad']
optim_df$ad <- aperm(optim_arr, c(3, 1, 2))[, , 'ad']

# xyplot(juv[, '2014'] ~ max_samps_tot,	optim_df)
# xyplot(juv[, '2015'] ~ max_samps_tot,	optim_df)
# xyplot(juv[, '2016'] ~ max_samps_tot,	optim_df)
# xyplot(juv[, '2017'] ~ max_samps_tot,	optim_df)
# xyplot(juv[, '2018'] ~ max_samps_tot,	optim_df)


pc <- c("mediumseagreen", "deepskyblue", "royalblue", "darkorchid")
#c("green3", "deepskyblue", "deeppink", "darkorchid")

par(mar = par()$mar + c(0,-0.5,0,-7))
avg_array = array(dim = c(5, 3, 9), dimnames = list(c(2014:2018), c("juv", "subad", "ad"), NULL))
par(mfrow=c(3, 3))
samp_ind = 1
for (samps in c(250, (3:10 * 100))) {
  d <- apply(simplify2array(optim_df[which(optim_df$max_samps_tot==samps), 'pars']), c(1, 2), mean)
  d <- d/sum(d)
  
  ### I don't know which is the more proper way to do the avg prop calc
  # d <- simplify2array(optim_df[which(optim_df$max_samps_tot==samps), 'pars'])
  # x <- apply(d, 3, sum)
  # for (i in 1:10) {
  #   d[ , , i] <- d[ , , i]/x[i]
  # }
  # d <- apply(d, c(1, 2), mean)
  
  plot(2014:2018, d[, 'juv'], col=pc[1], pch=19, ylim=c(0, 0.35),  bty='L',
       main=paste0("Maximum ", samps, " total samples"), xlab="Year", ylab="Proportion of total samples")
  grid(nx=NA, ny=NULL)
  lines(2014:2018, d[, 'juv'], col=pc[1], pch=19)
  points(2014:2018, d[, 'subad'], col=pc[2], pch=19)
  lines(2014:2018, d[, 'subad'], col=pc[2], pch=19)
  points(2014:2018, d[, 'ad'], col=pc[3], pch=19)
  lines(2014:2018, d[, 'ad'], col=pc[3], pch=19)
  points(2014:2018, apply(d, 1, sum), col=pc[4], pch=19)
  lines(2014:2018, apply(d, 1, sum), col=pc[4], pch=19)
  
  avg_array[ , , samp_ind] <- d
  samp_ind = samp_ind + 1
}

par(mfrow=c(1,1), mar = par()$mar + c(0,0.5,0,7))
overall_avg <- apply(avg_array, c(1, 2), mean)
plot(2014:2018, overall_avg[, 'juv'], col=pc[1], pch=19, ylim=c(0, 0.3), bty='L',
     main="Overall Average", xlab="Year", ylab="Proportion of total samples")
grid(nx=NA, ny=NULL)
lines(2014:2018, overall_avg[, 'juv'], col=pc[1])
points(2014:2018, overall_avg[, 'subad'], col=pc[2], pch=19)
lines(2014:2018, overall_avg[, 'subad'], col=pc[2])
points(2014:2018, overall_avg[, 'ad'], col=pc[3], pch=19)
lines(2014:2018, overall_avg[, 'ad'], col=pc[3])
points(2014:2018, apply(overall_avg, 1, sum), col=pc[4], pch=19)
lines(2014:2018, apply(overall_avg, 1, sum), col=pc[4])
legend(x = "right", legend = c("Juvenile", "Sub-adult", "Adult", "Total"), 
       col=pc, pch=19, lty=1, xpd = T, inset=c(-0.31,0), bty = "n")
par(mar = par()$mar + c(0,-0.5,0,-7))


overall_avg3 <- overall_avg

mean(optim_df$min_cv)
min(optim_df[which(optim_df$max_samps_tot==500), 'min_cv'])

write.csv(optim_df[, 1:3], "hpc_out_var_props_inc_all", row.names = FALSE)
write.csv(aperm(optim_arr, c(3, 1, 2)), "hpc_out_var_props_inc_all_mats")


write.csv(optim_df$juv, "hpc_out_var_props_juv", row.names = FALSE)
write.csv(optim_df$subad, "hpc_out_var_props_subad", row.names = FALSE)
write.csv(optim_df$ad, "hpc_out_var_props_ad", row.names = FALSE)
