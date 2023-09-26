import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging
import requests
from scipy.stats import t

logging.basicConfig(filename='OpenGwas_DataLoad.log', 
                    encoding='utf-8', 
                    level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-8s %(message)s \n',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

input_dir='../NDimClustInputs/BP/'
trait_input_csv = "TraitsData.csv"
trait_df=pd.read_csv(input_dir+trait_input_csv,sep=',',header=0)
i=0
trait_list=[]
cols=["pheno_category","phenotype","n_complete_samples","description",
      "variable_type","source","n_non_missing","n_missing",
      "n_controls","n_cases","n_eff","open_gwas_id",
      "keywords","nsnp"]
trait_info_df=pd.DataFrame(columns=cols)
gwas_data = igd.gwasinfo()
gwas_df = pd.DataFrame.from_dict(gwas_data, orient='index')
# Find the outcome snps
outcome_trait_row=trait_df[trait_df.Type=='Outcome'].index
outcome_description=trait_df.at[i,'Description'].replace(';',' ')
outcome_studies = gwas_df[gwas_df['trait'].str.contains(outcome_description, case=False)]
# Choose the outcome study with the most SNPs
max_nsnp=max(outcome_studies['nsnp'])
outcome_entry_row=outcome_studies[outcome_studies.nsnp==max_nsnp].index[0]
outcome_id=outcome_studies.at[outcome_entry_row,'id']
init_study=pd.DataFrame.from_dict(igd.tophits([outcome_id], 
                                              pval=5e-15, clump=1, r2=0.1, 
                                              kb=10000, force_server=False, 
                                              access_token='NULL')).sort_values(by='p',
                                                                                 ascending=True)
init_snp_list=init_study.rsid
for i in trait_df.index:
    if (trait_df.at[i,'Type']=='Outcome'): 
        pheno_cat='Outcome'
    else:
        pheno_cat='Exposure'
    tl=str(trait_df.at[i,'Phenotype'])
    t_entry=trait_df[trait_df.Phenotype==tl]
    trait_description=t_entry.at[i,'Description'].replace(';',' ')
    trait_studies = gwas_df[gwas_df['trait'].str.contains(trait_description, case=False)]
    if len(trait_studies)==0: 
        print('No matching studies found for trait ',tl)
        continue
    max_nsnp=0
    # Find the overlap in SNPS between the outcome and the trait.
    # Use the GWAS with the largest overlap.
    for study_id in trait_studies.id:
        assocs = pd.DataFrame.from_dict(igd.associations(init_snp_list,
                                                 [study_id], 
                                                 proxies=1, 
                                                 r2=0.8, 
                                                 align_alleles=1, 
                                                 palindromes=1, 
                                                 maf_threshold=0.3)
                                                )
        assoc_snps=len(assocs)
        if assoc_snps>max_nsnp:
            max_nsnp=assoc_snps
            study_max_entry=trait_studies[trait_studies.id==study_id]
    study_row           =study_max_entry.index[0]
    study_description   =trait_studies.at[study_row,'trait']
    study_id            =trait_studies.at[study_row,'id']
    study_sample_size   =trait_studies.at[study_row,'sample_size']
    study_controls      =trait_studies.at[study_row,'ncontrol']
    new_trait=pd.DataFrame({'phenotype':tl,'description':study_description,
               "n_complete_samples":study_sample_size,
               "variable_type":"",
               "source":"","n_non_missing":"","n_missing":"",
               "n_controls":study_controls,
               "n_cases":"","n_eff":"",
               "nsnp":max_nsnp,
               "keywords":[str(t_entry.at[i,'Keywords']).replace(';',' ')],
               "pheno_category":pheno_cat,
               "open_gwas_id":trait_studies.at[study_row,'id']})
    trait_info_df=pd.concat([trait_info_df,new_trait],ignore_index=True)
    # Get the trait id
    #logging.info(new_trait)
    i+=1
    #init_study=igd.gwasinfo([t_entry.OpenGwas_ID[i]])
trait_info_df.to_csv(input_dir+'/trait_info_nfil.csv',na_rep='NAN')
ic(trait_info_df)
#ic(igd.tophits([outcome_id]))
init_study=pd.DataFrame.from_dict(igd.tophits([outcome_id], 
                                              pval=5e-15, clump=1, r2=0.1, 
                                              kb=10000, force_server=False, 
                                              access_token='NULL')).sort_values(by='p',
                                                                                 ascending=True)
#logging.info(init_study.head())
#ic(init_study.head())
init_snp_list=init_study.rsid
# ACCESS DATA VIA API
# ieugwaspy package doesn't appear to be working. Trially GWAS catalog API
#apiurl="https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"+trait_df.Trait_ID[0]
#response = requests.get(apiurl)
#print(apiurl)
#print(response)
#print(response.json())
#apiurl="https://www.ebi.ac.uk/gwas/rest/api/studies/"+trait_df.All_Study_IDs[1]
#response = requests.get(apiurl)
#print(apiurl)
#print(response)
#print(response.json())
# Create Dataframes for beta, se, pval and tstat
# Fill with incrementing axis.
#------------------------------------
# Initialise the Dataframes with the column names as the trait labels
beta_df =pd.DataFrame(index=init_snp_list)
se_df   =pd.DataFrame(index=init_snp_list)
pval_df =pd.DataFrame(index=init_snp_list)
tstat_df=pd.DataFrame(index=init_snp_list)
for i_tr in trait_info_df.index:
    gwas_id=trait_info_df.at[i_tr,'open_gwas_id']
    trait  =trait_info_df.at[i_tr,'phenotype']
    assocs = pd.DataFrame.from_dict(igd.associations(init_snp_list,
                                                 [gwas_id], 
                                                 proxies=1, 
                                                 r2=0.8, 
                                                 align_alleles=1, 
                                                 palindromes=1, 
                                                 maf_threshold=0.3)
                                                ).sort_values(by='p', ascending=True)
    snp_list = assocs[assocs.id==gwas_id].rsid
    beta_df[trait]   =pd.Series(dtype=float)
    se_df[trait]     =pd.Series(dtype=float)
    pval_df[trait]   =pd.Series(dtype=float)
    tstat_df[trait]  =pd.Series(dtype=float)
    # For each snp find the association score for the current trait.
    # Trait information received from assocs using gwas_id but stored in matrix
    # by trait label
    for snp in snp_list:
        data_row=assocs[(assocs.rsid==snp) & (assocs.id==gwas_id)].index[0]
        beta_df.at[snp,trait]    =assocs.at[data_row,'beta']
        se_df.at[snp,trait]      =assocs.at[data_row,'se']
        pval_df.at[snp,trait]    =assocs.at[data_row,'p']
beta_df.to_csv(input_dir+'/unstdBeta_df.csv',na_rep='NAN')
se_df.to_csv(input_dir+'/unstdSE_df.csv',na_rep='NAN')
pval_df.to_csv(input_dir+'/pval_df.csv',na_rep='NAN')
tstat_df.to_csv(input_dir+'/tstat_df.csv',na_rep='NAN')
print(beta_df)


