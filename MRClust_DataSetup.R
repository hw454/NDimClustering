
dataname=out_dat

names(dataname)[names(dataname) == 'beta'] <- 'bx'
names(dataname)[names(dataname) == 'se'] <- 'bxse'
dataname$by = out_dat2$beta.outcome
dataname$byse = out_dat2$se.outcome
##############################################################
# estimates of the G-X associations and their standard errors
##############################################################
bx = dataname$bx
bxse = dataname$bxse
##############################################################
# estimates of the G-Y associations and their standard errors
##############################################################
by = dataname$by
byse = dataname$byse
###########################################
# ratio-estimates and their standard errors
###########################################
ratio_est = by/bx
ratio_est_se = byse/abs(bx)
#####################
# Names for the snps
#####################
snp_names = dataname$SNP
