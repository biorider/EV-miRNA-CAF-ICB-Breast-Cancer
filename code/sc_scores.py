import pandas as pd, numpy as np
from scipy import stats
d=pd.read_csv('/home/claude/data/02_scrnaseq_caf/table5_caf_metabolic_scores.csv')
rows=[]
for c in ['Glycolysis','OXPHOS','FAO','CAF_Activation','miRNA_Target_Suppression','Warburg_Index']:
    a=d[d.group=='Cancer_CAF'][c]; b=d[d.group=='Normal_Fibroblast'][c]
    u,p=stats.mannwhitneyu(a,b,alternative='two-sided')
    rbc=2*u/(len(a)*len(b))-1
    rows.append([c,a.mean(),a.std(),b.mean(),b.std(),a.mean()-b.mean(),rbc,p])
o=pd.DataFrame(rows,columns=['score','mean_CAF','sd_CAF','mean_NF','sd_NF','diff','rank_biserial','p']); print(o.to_string()); o.to_csv('out/sc_pathway_scores.csv',index=False)
print(d.group.value_counts())
