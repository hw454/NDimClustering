remotes::install_github("mrcieu/ieugwasr")
library("ieugwasr")

# Set the input location and datafile.
input_dir <- "../NDimClustInputs/BMI_CAD/"
trait_input_csv <-  "TraitsData.csv"

get_study <- function(snp_list, id_list, trait_studies) {
  max_nsnp <- 0
  for (study_id in id_list){
    print(paste("associations with", study_id))
    assocs <- ieugwasr::associations(snp_list,
                                    id = study_id,
                                    proxies = 1,
                                    r2 = 0.1,
                                    align_alleles = 1,
                                    palindromes = 1,
                                    maf_threshold = 0.3)
    assoc_snps <- length(assocs)
    if (assoc_snps > max_nsnp) {
      max_nsnp <- assoc_snps
      study_max_entry <- trait_studies[trait_studies$id == study_id]
    }
  }
  out_list <- list("study" = study_max_entry, "assocs" = assocs)
  return(out_list)
}
set_data_row <- function(trait, study_id, snp_list,
                         df_dict, trait_info_df, assocs) {
  beta_df <- df_list["beta"]
  pval_df <- df_list["pval"]
  se_df <- df_list["se"]
  tstat_df <- df_list["tstat"]
  print(paste("Get snp matrix data col for ", trait))
  trait_snp_list <- assocs[assocs$id == study_id, "rsid"]
  # For each snp find the association score for the current trait.
  # Trait information received from assocs using gwas_id but stored in matrix
  # by trait label
  for (snp in trait_snp_list) {
    data_row <- row.names(assocs[(assocs$rsid == snp) &
                                 (assocs$id == study_id)])[1]
    beta_df[snp, trait]  <- assocs[data_row, "beta"]
    se_df[snp, trait]    <- assocs[data_row, "se"]
    pval_df[snp, trait]  <- assocs[data_row, "p"]
    tstat_df[snp, trait] <- assocs[data_row, "tstat"]
  }
  df_list <- list("beta" = beta_df,
             "se" = se_df,
             "pval" = pval_df,
             "tstat" = tstat_df)
  return(df_list)
}

# MAIN SCRIPT
#logging.basicConfig(filename = "OpenGwas_DataLoad.log",
#                    encoding = "utf-8",
#                    level = logging.DEBUG,
 #                   format = "%(asctime)s %(levelname)-8s %(message)s \n",
  #                  datefmt= "%Y-%m-%d %H:%M:%S",
   #                 filemode="w")
# Read in the desired traits.
trait_df <- read.csv(paste0(input_dir, trait_input_csv),
                        sep = ",",
                        header = TRUE)
print(trait_df)

# Set the column names for the trait dataframe
cols <- list("pheno_category", "phenotype", "n_complete_samples", "description",
      "variable_type", "source", "n_non_missing", "n_missing",
      "n_controls", "n_cases", "n_eff", "open_gwas_id",
      "keywords", "nsnp")
trait_info_df <- data.frame()

# Load GWAS info
print("GWAS info")
gwas_data <- ieugwasr::gwasinfo()
print(gwas_data)
# Find the outcome study
outcome_trait_row <- row.names(trait_df[trait_df$Type == "Outcome"])[1]
outcome_description <- replace(trait_df[outcome_trait_row,
                                  "Description"], ";", " ")
outcome_studies <- gwas_df[grepl(gwas_df["trait"],
                                 outcome_description,
                                 ignore.case = TRUE)]
# Choose the outcome study with the most SNPs
max_nsnp_o <- max(outcome_studies["nsnp"])
outcome_entry_row <- row.names(outcome_studies[
                                  outcome_studies$nsnp == max_nsnp_o])[1]
outcome_id <- outcome_studies[outcome_entry_row, "id"]
# Find the exposure snps
exposure_trait_row <- row.names(trait_df[trait_df$Type == "Exposure"])[1]
exposure_description <- replace(trait_df.at[exposure_trait_row,
                                  "Description"], ";", " ")
exposure_studies <- gwas_df[grepl(gwas_df["trait"],
                                  exposure_description, ignore.case = TRUE)]
# Choose the outcome study with the most SNPs
max_nsnp_e <- max(exposure_studies["nsnp"])
exposure_id <- row.names(exposure_studies[
                                exposure_studies$nsnp == max_nsnp_e])[1]
# Define the initial study by the top hits for the exposure
ic("tophits ", exposure_id)
init_study <- ieugwasr::tophits(id = exposure_id,
                          pval = 5e-8,
                          clump = 1,
                          r2 = 0.01,
                          force_server = FALSE,
                          access_token = "NULL")
# Initial SNP list is the SNPS in exposure study.
init_snp_list <- init_study$rsid
print("Number of snps in initial study", length(init_snp_list))
# Loop through each trait and locate the associations between the trait,
# outcome and the snp_list.
# Initialise the Dataframes with the column names as the trait labels
beta_df <- data.frame(row.names = init_snp_list)
se_df   <- data.frame(row.names = init_snp_list)
pval_df <- data.frame(row.names = init_snp_list)
tstat_df <- data.frame(row.names = init_snp_list)
data_list <- list("beta" = beta_df,
           "se" = se_df,
           "pval" = pval_df,
           "tstat" = tstat_df)
for (i in row.names(trait_df)) {
    if (trait_df[i, "Type"] == "Outcome") {
        pheno_cat <- "Outcome"
    } else if (trait_df[i, "Type"] == "Exposure") {
        pheno_cat <- "Exposure"
    } else {
        pheno_cat <- "Trait"
    }
    tl <- str(trait_df[i, "Phenotype"])
    t_entry <- trait_df[trait_df$Phenotype == tl]
    trait_description <- replace(t_entry[i, "Description"], ";", " ")
    trait_studies <- gwas_df[grepl(gwas_df["trait"],
                                    trait_description, ignore.case = False)]
    if (length(trait_studies) == 0) {
        print(paste("No matching studies found for trait", tl))
        next
    }
    # Find the overlap in SNPS between the outcome and the trait.
    # Use the GWAS with the largest overlap.
    out <- get_study(init_snp_list,
                                id_list = trait_studies$id,
                                trait_studies = trait_studies)
    study_max_entry <- out$study
    assocs <- out$assocs
    study_row <- row.names(study_max_entry)[1]
    max_nsnp  <- trait_studies[study_row, "nsnp"]
    study_description <- trait_studies[study_row, "trait"]
    study_id <- trait_studies[study_row, "id"]
    study_sample_size <- trait_studies[study_row, "sample_size"]
    study_controls <- trait_studies[study_row, "ncontrol"]
    study_description <- study_description.replace('"', "")
    study_description <- study_description.replace(";", "-")
    study_description <- study_description.replace(",", "-")
    new_trait <- data.frame("phenotype" = tl,
               "description" = study_description,
               "n_complete_samples" = study_sample_size,
               "variable_type" = "",
               "source" = "",
               "n_non_missing" = "",
               "n_missing" = "",
               "n_controls" = study_controls,
               "n_cases" = "",
               "n_eff" = "",
               "nsnp" = max_nsnp,
               "keywords" = replace(str(t_entry[i,"Keywords"]),";"," "),
               "pheno_category" = pheno_cat,
               "open_gwas_id" = study_id)
    trait_info_df <- rbind(trait_info_df, new_trait)
    # Get the trait id
    i <- i + 1
    # ADDFEATURE - With trait that's assigned add data to matrices
    study_row_ind <- row.names(trait_info_df)[length(trait_info_df)]
    data_list <- set_data_row(trait = tl,
                           study_id = study_id,
                           snp_list = init_snp_list,
                           df_list = data_list,
                           trait_info_df = trait_info_df,
                           assocs = assocs)
}
print(paste("trait df", trait_info_df))
# Save data to csv files
#------------------------------------
write.csv(trait_info_df, paste0(input_dir, "/trait_info_nfil.csv"), na = "NAN")
write.csv(beta_df, paste0(input_dir, "/unstdBeta_df.csv"), na = "NAN")
write.csv(se_df, paste0(input_dir, "/unstdSE_df.csv"), na = "NAN")
write.csv(pval_df, paste0(input_dir + "/pval_df.csv"), na = "NAN")
write.csv(tstat_df, paste0(input_dir + "/tstat_df.csv"), na = "NAN")
print(paste("number of points", nrow(beta_df)))