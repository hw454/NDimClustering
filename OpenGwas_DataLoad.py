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

class trait_info:
    def __init__(s,label):
        s.label=label
        s.description=""
        s.keywords=[label,]
        s.study_ids=[]
        s.best_study_id=""
        s.nsnp=0
    def __str__(s):
        label_str="Trait object "+s.label
        description_str="\nDescription: "+s.description
        id_str="\nStudy with most SNPs: "+s.best_study_id
        return label_str+description_str+id_str
    def add_description(s,text):
        "Add the text `text` to the current trait_info description"
        s.description=s.description+text
        return
    def clear_description(s):
        s.description=""
        return
    def add_keyword(s,k):
        s.keywords=s.keywords+[k]
        return
    def add_keywords_list(s,klist):
        s.keywords=s.keywords+klist
        return
    def add_study_ids(s,idlist):
         s.study_ids=s.study_ids+idlist
         return
    #def get_study_id(s,gwas_df):
    #     trait_studies = gwas_df[all(gwas_df['trait'].str.contains(w, case=False) for w in s.description)]
    #     s.add_study_ids(trait_studies.id)
    #     for id in trait_studies.id:
    #         N=trait_studies[trait_studies.id==id,'nsnp']
    #         if N>s.nsnp:
    #             s.best_study_id=id
    #     return s.study_ids
input_dir='../NDimClustInputs/'
trait_input_csv = "TraitsData.csv"
trait_df=pd.read_csv(input_dir+trait_input_csv,sep=';')
i=0
trait_list=[]
#data = igd.gwasinfo()
#print('Open gwas loaded')
#df = pd.DataFrame.from_dict(data, orient='index')
for tl in trait_df.Label:
    trait=trait_info(tl)
    t_entry=trait_df[trait_df.Label==tl]
    trait.add_description(t_entry.Description[i])
    # Get the trait id
    # trait.get_study_id()
    trait.add_study_ids(t_entry.All_Study_IDs[i].split(','))
    trait.add_keywords_list(t_entry.Keywords[i].split(','))
    trait.best_study_id=str(t_entry.OpenGwas_ID[i])
    trait_list=trait_list+[trait]
    logging.info(str(trait))
    ic(str(trait))
    i+=1
    #init_study=igd.gwasinfo([t_entry.OpenGwas_ID[i]])
init_study=pd.DataFrame.from_dict(igd.tophits([trait_df.OpenGwas_ID[0]], 
                                              pval=5e-15, clump=1, r2=0.1, 
                                              kb=10000, force_server=False, 
                                              access_token='NULL')).sort_values(by='p',
                                                                                 ascending=True)
logging.info(init_study.head())
ic(init_study.head())
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
gwas_ids=trait_df[trait_df.OpenGwas_ID.notnull()].OpenGwas_ID
ic(init_snp_list)
ic(gwas_ids)
assocs = pd.DataFrame.from_dict(igd.associations(init_snp_list,
                                                 gwas_ids, 
                                                 proxies=1, 
                                                 r2=0.8, 
                                                 align_alleles=1, 
                                                 palindromes=1, 
                                                 maf_threshold=0.3)
                                                ).sort_values(by='p', ascending=True)
ic(assocs.columns)
# Create Dataframes for beta, se, pval and tstat
# Fill with incrementing axis.
#------------------------------------
# Initialise the Dataframes with the column names as the trait labels
beta_df =pd.DataFrame(index=init_snp_list)
se_df   =pd.DataFrame(index=init_snp_list)
pval_df =pd.DataFrame(index=init_snp_list)
tstat_df=pd.DataFrame(index=init_snp_list)

for trait_id in gwas_ids:
    assocs = pd.DataFrame.from_dict(igd.associations(init_snp_list,
                                                 [trait_id], 
                                                 proxies=1, 
                                                 r2=0.8, 
                                                 align_alleles=1, 
                                                 palindromes=1, 
                                                 maf_threshold=0.3)
                                                ).sort_values(by='p', ascending=True)
    beta_df[trait_id]   =pd.Series(dtype=float)
    se_df[trait_id]     =pd.Series(dtype=float)
    pval_df[trait_id]   =pd.Series(dtype=float)
    tstat_df[trait_id]  =pd.Series(dtype=float)
    beta_df[trait_id]   =assocs[assocs.id==trait_id].beta
    se_df[trait_id]     =assocs[assocs.id==trait_id].se
    pval_df[trait_id]   =assocs[assocs.id==trait_id].p
    #degf=igd.gwasinfo([trait_id])['sample_size']-1
    #tstat_df[trait_id]  =[t.ppf(p,degf) for p in pval_df[trait_id]]
beta_df.to_csv(input_dir+'/unstdBeta_df.csv')
se_df.to_csv(input_dir+'/se_df.csv')
pval_df.to_csv(input_dir+'/pval_df.csv')
tstat_df.to_csv(input_dir+'/tstat_df.csv')


