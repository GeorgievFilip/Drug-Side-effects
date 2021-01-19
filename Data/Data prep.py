import urllib.request, urllib.error, urllib.parse
import json
import pandas as pd
from tqdm import tqdm
import pubchempy as pcp
import numpy as np
from sklearn.preprocessing import MultiLabelBinarizer


sider = pd.read_csv('Sider(v4).csv', dtype={'drug_rxnorn_id': object})

offsides = (
pd.read_csv('OFFSIDES.csv.gz', compression='gzip', error_bad_lines=False, low_memory=False, usecols = [
    'drug_rxnorn_id',
    'condition_concept_name',
    'condition_meddra_id'  
])
)

API_KEY = ""

def get_json(url):
    opener = urllib.request.build_opener()
    opener.addheaders = [('Authorization', 'apikey token=' + API_KEY)]
    return json.loads(opener.open(url).read())

def get_class(drug_id):
    search = None
    try:
        tree = get_json('http://data.bioontology.org/ontologies/MEDDRA/classes/http%3A%2F%2Fpurl.bioontology.org%2Fontology%2FMEDDRA%2F{}/tree'.format(drug_id))
        for i in range(len(tree)):
            if str(drug_id) in str(tree[i]):
                search = tree[i]['prefLabel']
    except:
        pass
    return search


offsides_condition = offsides.drop_duplicates(subset=['condition_meddra_id'], keep='first')
offsides_condition['condition_concept_class'] = offsides_condition['condition_meddra_id'].apply(get_class)


offsides_condition = offsides_condition[offsides_condition['condition_concept_class'].map(offsides_condition['condition_concept_class'].value_counts()) > 1]
offsides = offsides.merge(offsides_condition[['condition_meddra_id', 'condition_concept_class']])
offsides = offsides.dropna()


unique_drugs = offsides[['drug_rxnorn_id']].drop_duplicates(keep='first')
unique_drugs['CanonicalSMILES'] = unique_drugs['drug_rxnorn_id'].apply(lambda x: pcp.Compound.from_cid(x).canonical_smiles)
unique_drugs.replace(to_replace = 'CC(C)(C(=O)O)OC1=CC=C(C=C1)CCNC(=O)C2C=CC(=[ClH])C=C2', value = 'CN1C(S(=O)(=O)CCC1=O)C2=CC=C(C=C2)Cl', inplace = True)
offsides = offsides.merge(unique_drugs)


off = offsides[['drug_rxnorn_id', 'condition_concept_class', 'CanonicalSMILES']]
off.reset_index(drop=True, inplace=True)


labels = pd.DataFrame(off.groupby('drug_rxnorn_id')['condition_concept_class'].apply(lambda x: tuple(set(x.values))))
labels.reset_index(level=0, inplace=True)
mlb = MultiLabelBinarizer()
mlb.fit(labels['condition_concept_class'])


off = pd.concat([labels, pd.DataFrame(data=mlb.fit_transform(labels['condition_concept_class']), columns=mlb.classes_)], axis= 1).merge(off.drop_duplicates(subset=['drug_rxnorn_id'], keep='first')[['drug_rxnorn_id','CanonicalSMILES']])
off.drop('condition_concept_class', axis=1, inplace=True)

df = pd.concat([off,sider])
df = df.drop_duplicates(subset=['drug_rxnorn_id'], keep='first', ignore_index=True)
col = df.pop("CanonicalSMILES")
df.insert(1, col.name, col)

# Final dataset
df.to_csv('offsides&sider_labels.csv',index=False)

splits = np.array_split(df, 5)
for d in range(len(splits)):
    with open('smiles{}.smi'.format(d+1), 'w') as f:
        for i in splits[d]['CanonicalSMILES']:
            f.write("{}\n".format(i))