import ieugwaspy as igd
import pandas as pd
from icecream import ic
import logging


class trait_info:
    def __init__(s,label):
        s.label=label
        s.description=""
        s.keywords=[label,]
        s.study_id=[]
        s.best_study_id=""
        s.nsnp=0
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
    def get_study_id(s,gwas_df):
        trait_studies = gwas_df[all(gwas_df['trait'].str.contains(w, case=False) for w in s.description)]
        s.add_study_ids(trait_studies.id)
        for id in trait_studies.id:
            N=trait_studies[trait_studies.id==id,'nsnp']
            if N>s.nsnp:
                s.best_study_id=id
        return s.study_ids
    def __str__(s):
        label_str="Trait object "+s.label+"\n"
        description_str="Description: "+s.description+"\n"
        id_str="Study with most SNPs: "+s.best_study_id+"\n"
        return label_str+description_str+id_str



trait_label_list=["LDL","SBP","WHR","EI","WFP"]#,"SBD","QOL"]
trait_description_list=["Low Density Lipoprotein",
                    "Systolic  Blood Pressure",
                    "Waist to hip ratio",
                    "Energy intake",
                    "Walking for pleasure"]#,
                    #"Sedentary behaviour duration",
                    #"Quality of life"]
trait_keywords_list =[["Lipoprotein","Cholesterol"],
                      ["Blood_Pressure","Cardio"],
                      ["Waist_Hip_Ratio","Body_Composition","Adiposity"],
                      ["Energy_Intake"],
                      ["Walking_pleasure"]]
#                      ,["Healthy_Movement_Behaviour","Sedentary"],
#                      ["Quality_of_life"]]
trait_ids = ["ieu-b-110","ieu-b-5075","ebi-a-GCST008052","ukb-b-7323","ukb-b-7337"]
#,"GCST006913","EFO_0011014"]
i=0
trait_list=[]
#data = igd.gwasinfo()
#print('Open gwas loaded')
#df = pd.DataFrame.from_dict(data, orient='index')
for tl in trait_label_list:
    trait=trait_info(tl)
    trait.add_description(trait_description_list[i])
    # Get the trait id
    # trait.get_study_id()
    trait.best_study_id=trait_ids[i]
    trait_list=trait_list+[trait]
    print(trait)
    i+=1
init_study=igd.gwasinfo(trait_ids[0])
print(init_study)
""" assocs = pd.DataFrame.from_dict(igd.associations(init_snp_list,
                                                 trait_ids[1:], 
                                                 proxies=1, 
                                                 r2=0.8, 
                                                 align_alleles=1, 
                                                 palindromes=1, 
                                                 maf_threshold=0.3)
                                                ).sort_values(by='p', ascending=True) """