#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Import packages
import numpy as np
import re
import pandas as pd
from scipy import stats
import sys

#Define a function to parse all binding scores from the file of random peptides
def parse_random_peptide_scores(filename) :
    random_peptide_scores = []
    with open(filename) as file :
        for line in file :
            line = line.rstrip().split()
            if "PEPLIST" in line and "Number" not in line:
                random_peptide_scores.append(float(line[11])) #Binding score
                
    return random_peptide_scores

#Define a function to create a dataframe containing each peptide with its binding score and rank
def parse_binding_scores(filename) :
    peptides = []
    binding_scores = []
    with open(filename) as file :
        for line in file :
            if line.isspace() == True :
                continue
            line = line.rstrip().split()
            if (bool(re.match("\d+", line[0]))) == True:
                peptides.append(line[2]) #Peptide
                binding_scores.append(float(line[11])) #Binding score
    
    data = {
        "Peptide": peptides,
        "Binding Score": binding_scores,
    }
    
    return pd.DataFrame(data)

#Define a function that transforms a list of binding scores into a list of ranks using random peptide scores
def transform_scores(binding_scores, random_peptide_scores) :
    ranks = []
    for score in binding_scores :
        ranks.append(100 - stats.percentileofscore(random_peptide_scores, score))
        
    return ranks

#Load list of allele names
allelefile = sys.argv[1]
alleles = pd.read_csv(allelefile)

#Load list of virus names
viruses = pd.read_csv("viruses.csv")


#Iterate through alleles
for allele in alleles :
    random_peptide_scores = parse_random_peptide_scores("Random Peptides/Random_peptides_"+allele+".out")
    
    #Iterate through viruses
    for virus in viruses :
        df_parsing = parse_binding_scores("NetMHCpan Output Files/"+virus+".fasta_"+allele+".output")
        df_parsing.insert(0, "Virus", virus)
        df_parsing['Computed Rank'] = transform_scores(list(df_parsing.iloc[:,2]), random_peptide_scores)
        df_parsing.to_csv("Dataframes/"+allele+".csv", mode='a', index=False)

