#!/usr/bin/env python

# Description: Please create a short python script that queries cbioportal: http://www.cbioportal.org/

# Input:
# A gene name

# Output:
# A table listing all distinct genomic variants (chr, start_position, end_position, reference_allele, variant_allele)
# found in the given gene, along with the count of unique samples in which each variant was found.  
# Please query all samples in cbioportal, and please exclude variants with undefined coordinates.

from bravado.client import SwaggerClient
import pandas as pd
import numpy as np

cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False})
# put description or remove
print(cbioportal)

### Input gene symbol/name
symbol = 'TNFRSF10B'  # TNFRSF10B
genedic = cbioportal.Genes.getGeneUsingGET(geneId = symbol ).result()
geneid = genedic['entrezGeneId']


### Define function: retrieve mutations in one sampleListId
colName = ['Gene','chr','startPosition','endPosition','referenceAllele','variantAllele','sampleId']    
def mutations_in_one_sampleList(mp, sl, id):
    # Retrieve mutation dictionary from one SampleListId
    mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        molecularProfileId = mp,
        sampleListId = sl,
        entrezGeneId = id,
        projection ='DETAILED'
    ).result()
    
    # Pull out specific values we need from mutations dictionary
    mydf = pd.DataFrame()
    for mut in mutations:
        mutVal = [mut.gene['hugoGeneSymbol'], mut.chr, mut.startPosition,mut.endPosition
                , mut.referenceAllele, mut.variantAllele, mut.sampleId]
        new = pd.DataFrame(mutVal).T
        mydf = mydf.append(new, ignore_index=True)
    return(mydf)

### Iterate through all sampleListId in all studies. Output mutations as a dataframe
studies = cbioportal.Studies.getAllStudiesUsingGET().result()
study_mol_profile = dict()
study_sample_list = dict()
mutTable = pd.DataFrame()

for study in studies: 
    ID = study.studyId
    # description something when it's slow, print dot or study id
    print('...')
    mol_prolile_list = cbioportal.Molecular_Profiles.getAllMolecularProfilesInStudyUsingGET(studyId=ID).result()
    for profile in mol_prolile_list: 
        if profile.datatype == "MAF": 
            mol_profile = profile.molecularProfileId
    sample_lists = cbioportal.Sample_Lists.getAllSampleListsInStudyUsingGET(studyId=ID, sortBy= 'category').result()
    study_mol_profile[ID] = mol_profile
    study_sample_list[ID] = sample_lists[0].sampleListId
    new = mutations_in_one_sampleList(mp = mol_profile, sl = sample_lists[0].sampleListId, id = geneid)
    mutTable = mutTable.append(new, ignore_index = True)

### Check how many mutations we obtained
# add description if you want to output this
print(mutTable.shape)
print(len(study_mol_profile))
print(len(study_sample_list))

### Remove mutations with undefined coordinates
mutTable.columns= colName
mutTable = mutTable.dropna(axis=0, subset= ['chr','startPosition','endPosition'] )

print(mutTable.shape)

### Count unique samples where each mutation was found. And output as CSV file.
mutCount = mutTable.groupby(['Gene','chr', 'startPosition','endPosition','referenceAllele','variantAllele']).nunique()
mutCount = mutCount['sampleId']
mutCount.to_csv('output_variants_unique_samples_count.csv', sep = ',')
# print some description!
mutCount.to_csv('output_variants.csv', sep = ',')

