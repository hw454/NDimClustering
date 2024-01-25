import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging
import requests
from scipy.stats import t

# Set the input location and datafile.
exposure = "BMI"

# MAIN SCRIPT
logging.basicConfig(filename = "OpenGwas_DataLoad.log", 
                    encoding = "utf-8", 
                    level = logging.DEBUG,
                    format = "%(asctime)s %(levelname)-8s %(message)s \n",
                    datefmt= "%Y-%m-%d %H:%M:%S",
                    filemode="w")

# Set the column names for the trait dataframe
cols = ["pheno_category","phenotype","n_complete_samples","description",
      "variable_type","source","n_non_missing","n_missing",
      "n_controls","n_cases","n_eff","open_gwas_id",
      "keywords","nsnp"]
trait_info_df=pd.DataFrame(columns=cols)

# FIND THE STUDY BEST MATCHING THE EXPOSURE LABEL
# Load GWAS info
ic("GWAS info")
gwas_data = igd.gwasinfo()
gwas_df = pd.DataFrame.from_dict(gwas_data, orient="index")
# Find the exposure study
exposure_studies = gwas_df[gwas_df["trait"].str.contains(exposure, case = False)]
# Choose the exposure study with the most SNPs
max_nsnp_e = max(exposure_studies["nsnp"])
exposure_entry_row = exposure_studies[exposure_studies.nsnp == max_nsnp_e].index[0]
exposure_id = exposure_studies.at[exposure_entry_row,"id"]
# Define the initial study by the top hits for the exposure
ic("tophits ", exposure_id)
init_study = pd.DataFrame.from_dict(igd.tophits([exposure_id], 
                                              pval = 1e-18,
                                              clump = 1,
                                              r2 = 1e-3, 
                                              pop = "EUR",
                                              force_server = False, 
                                              access_token = "NULL"))
# Initial SNP list is the SNPS in exposure study.
init_snp_list = init_study.rsid
print("Number of snps in initial study", len(init_snp_list))
# Loop through each trait and locate the associations between the trait, outcome and the snp_list. 
# Initialise the Dataframes with the column names as the trait labels
phe_df = pd.DataFrame.from_dict(igd.phewas(init_snp_list)).sort_values(by='p', ascending=True)
study_ids = phe_df.id.unique()
trait_list = phe_df.trait.unique()

# Remove duplicates
id_trait_pairs = phe_df.loc[ : ,["id", "trait"]].unique()
duplicates = phe_df.loc[id_trait_pairs.trait.duplicated(keep = False)]
trait_rem_ids = []
for row in duplicates:
    trait = row.trait
    trait_dup = phe_df.loc[phe_df.trait == trait]
    trait_max_n = trait_dup.n.idemax()
    trait_keep_id = trait_dup.id[trait_max_n]
    rem_traits = trait_dup.drop(index = trait_keep_id)
    trait_rem_ids = trait_rem_ids + rem_traits.id
phe_df = phe_df.drop("id" in trait_rem_ids)

print(phe_df)

# Filter the traits by sample size
# Filter by high correlation to BMI r>0.75
for trait in trait_list:
    if trait.str.contains(exposure, case= False):
        pass
    #else:
        #rg = 
# Steiger-filtering. p-value>0.05/n_traits