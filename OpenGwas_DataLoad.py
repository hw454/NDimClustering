import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging
import requests
from scipy.stats import t

MAX_TRIES = 100

# Set the input location and datafile.
input_dir="../NDimClustInputs/BMI_CAD/"
trait_input_csv = "TraitsData.csv"

def get_study(snp_list, id_list):
  max_nsnp = 0
  for study_id in id_list:
    ic("associations with ", study_id)
    assocs = pd.DataFrame.from_dict(igd.associations(snp_list,
                                                 [study_id], 
                                                 proxies = 1, 
                                                 r2 = 0.5, 
                                                 align_alleles = 1, 
                                                 palindromes = 1, 
                                                 maf_threshold = 0.3)
                                                )
    assoc_snps=len(assocs)
    if assoc_snps>max_nsnp:
      max_nsnp=assoc_snps
      study_max_entry=trait_studies[trait_studies.id==study_id]
  return(study_max_entry, assocs)
def set_data_row(trait, study_id, snp_list, df_dict, trait_info_df, assocs):
  beta_df = df_dict["beta"]
  pval_df = df_dict["pval"]
  se_df = df_dict["se"]
  tstat_df = df_dict["tstat"]
  ic("Get snp matrix data col for ", trait)
  trait_snp_list = assocs[assocs.id == study_id].rsid
  # For each snp find the association score for the current trait.
  # Trait information received from assocs using gwas_id but stored in matrix
  # by trait label
  for snp in trait_snp_list:
    data_row = assocs[(assocs.rsid == snp) & (assocs.id == gwas_id)].index[0]
    ic(assocs.columns)
    beta_df.at[snp, trait]  = assocs.at[data_row,"beta"]
    se_df.at[snp, trait]    = assocs.at[data_row,"se"]
    pval_df.at[snp, trait]  = assocs.at[data_row,"p"]
    tstat_df.at[snp, trait] = assocs.at[data_row,"tstat"]
  df_dict = {"beta": beta_df,
             "se": se_df,
             "pval": pval_df,
             "tstat": tstat_df}
  return df_dict

# MAIN SCRIPT
logging.basicConfig(filename = "OpenGwas_DataLoad.log", 
                    encoding = "utf-8", 
                    level = logging.DEBUG,
                    format = "%(asctime)s %(levelname)-8s %(message)s \n",
                    datefmt= "%Y-%m-%d %H:%M:%S",
                    filemode="w")
# Read in the desired traits.
trait_df = pd.read_csv(input_dir + trait_input_csv,
                        sep = ",",
                        header = 0)

# Set the column names for the trait dataframe
cols = ["pheno_category","phenotype","n_complete_samples","description",
      "variable_type","source","n_non_missing","n_missing",
      "n_controls","n_cases","n_eff","open_gwas_id",
      "keywords","nsnp"]
trait_info_df=pd.DataFrame(columns=cols)

# Load GWAS info
ic("GWAS info")
gwas_data = igd.gwasinfo()
gwas_df = pd.DataFrame.from_dict(gwas_data, orient="index")
# Find the outcome study
outcome_trait_row = trait_df[trait_df.Type == "Outcome"].index[0]
outcome_description = trait_df.at[outcome_trait_row,
                                  "Description"].replace(';',' ')
outcome_studies = gwas_df[gwas_df["trait"].str.contains(outcome_description, case = False)]
# Choose the outcome study with the most SNPs
max_nsnp_o = max(outcome_studies["nsnp"])
outcome_entry_row = outcome_studies[outcome_studies.nsnp == max_nsnp_o].index[0]
outcome_id = outcome_studies.at[outcome_entry_row,"id"]
# Find the exposure snps
exposure_trait_row = trait_df[trait_df.Type == "Exposure"].index[0]
exposure_description = trait_df.at[exposure_trait_row,
                                  "Description"].replace(';',' ')
exposure_studies = gwas_df[gwas_df["trait"].str.contains(exposure_description, case = False)]
# Choose the outcome study with the most SNPs
max_nsnp_e = max(exposure_studies["nsnp"])
exposure_id = exposure_studies[exposure_studies.nsnp == max_nsnp_e].index[0]
# Define the initial study by the top hits for the exposure
ic("tophits ", exposure_id)
init_study = pd.DataFrame.from_dict(igd.tophits([exposure_id], 
                                              pval = 1e-7,
                                              clump = 0,
                                              r2 = 0.3, 
                                              force_server = False, 
                                              access_token = "NULL"))
# Initial SNP list is the SNPS in exposure study.
init_snp_list = init_study.rsid
print("Number of snps in initial study", len(init_snp_list))
# Loop through each trait and locate the associations between the trait, outcome and the snp_list. 
# Initialise the Dataframes with the column names as the trait labels
beta_df =pd.DataFrame(index=init_snp_list)
se_df   =pd.DataFrame(index=init_snp_list)
pval_df =pd.DataFrame(index=init_snp_list)
tstat_df=pd.DataFrame(index=init_snp_list)
data_dict = {"beta": beta_df,
           "se": se_df,
           "pval": pval_df,
           "tstat": tstat_df}
for i in trait_df.index:
    if (trait_df.at[i,"Type"] == "Outcome"): 
        pheno_cat = "Outcome"
    elif(trait_df.at[i,"Type"] == "Exposure"): 
        pheno_cat = "Exposure"
    else:
        pheno_cat = "Trait"
    tl = str(trait_df.at[i,"Phenotype"])
    t_entry = trait_df[trait_df.Phenotype==tl]
    trait_description = t_entry.at[i,"Description"].replace(';',' ')
    trait_studies = gwas_df[gwas_df["trait"].str.contains(trait_description, case = False)]
    if len(trait_studies) == 0: 
        print("No matching studies found for trait ",tl)
        continue
    # Find the overlap in SNPS between the outcome and the trait.
    # Use the GWAS with the largest overlap.
    study_max_entry, assocs = get_study(init_snp_list, 
                                id_list = trait_studies.id)
    study_row           =study_max_entry.index[0]
    max_nsnp            = trait_studies.at[study_row,"nsnp"]
    study_description   = trait_studies.at[study_row,"trait"]
    study_id            = trait_studies.at[study_row,"id"]
    study_sample_size   = trait_studies.at[study_row,"sample_size"]
    study_controls      = trait_studies.at[study_row,"ncontrol"]
    study_description = study_description.replace('"',"")
    study_description = study_description.replace(";","-")
    study_description = study_description.replace(",","-")
    new_trait=pd.DataFrame({"phenotype":tl,
               "description": study_description,
               "n_complete_samples": study_sample_size,
               "variable_type": "",
               "source": "",
               "n_non_missing": "",
               "n_missing": "",
               "n_controls": study_controls,
               "n_cases": "",
               "n_eff": "",
               "nsnp": max_nsnp,
               "keywords": [str(t_entry.at[i,"Keywords"]).replace(';',' ')],
               "pheno_category": pheno_cat,
               "open_gwas_id": study_id})
    trait_info_df=pd.concat([trait_info_df,new_trait],ignore_index=True)
    # Get the trait id
    #logging.info(new_trait)
    i+=1
    # ADDFEATURE - With trait that's assigned add data to matrices
    study_row_ind = trait_info_df.index[-1]
    data_dict = set_data_row(trait = tl,
                           study_id = study_id,
                           snp_list = init_snp_list,
                           df_dict = data_dict,
                           trait_info_df = trait_info_df,
                           assocs = assocs)
    
ic("trait df", trait_info_df)
# Save data to csv files
#------------------------------------
trait_info_df.to_csv(input_dir+"/trait_info_nfil.csv", na_rep="NAN")
data_dict.beta.to_csv(input_dir + "/unstdBeta_df.csv", na_rep = "NAN")
data_dict.se.to_csv(input_dir + "/unstdSE_df.csv", na_rep = "NAN")
data_dict.pval.to_csv(input_dir + "/pval_df.csv", na_rep = "NAN")
data_dict.tstat.to_csv(input_dir + "/tstat_df.csv", na_rep = "NAN")
ic("number of points", beta_df.shape[0])


