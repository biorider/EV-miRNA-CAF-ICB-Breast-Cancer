import pandas as pd, numpy as np, sys
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from scipy import stats
d=pd.read_csv(sys.argv[1])
print(d.shape, d.vital_status.value_counts().to_dict())
d=d.dropna(subset=['os_months','os_event'])
print('n with OS',len(d),'events',int(d.os_event.sum()))
print('groups',d.Signature_Group.value_counts().to_dict())
print(d.groupby('Signature_Group').os_event.sum().to_dict())
h=d[d.Signature_Group=='High']; l=d[d.Signature_Group=='Low']
lr=logrank_test(h.os_months,l.os_months,h.os_event,l.os_event); print('logrank median split p=%.3f'%lr.p_value)
for g,s in [('High',h),('Low',l)]:
    k=KaplanMeierFitter().fit(s.os_months,s.os_event); print(g,'median OS',k.median_survival_time_)
mt=multivariate_logrank_test(d.os_months,d.Signature_Tertile,d.os_event); print('tertile logrank p=%.3f'%mt.p_value, d.Signature_Tertile.value_counts().to_dict())
# stage
def st(s):
    if not isinstance(s,str): return np.nan
    s=s.replace('Stage ','')
    for k,v in [('IV',4),('III',3),('II',2),('I',1)]:
        if s.startswith(k): return v
    return np.nan
d['stage_num']=d.ajcc_stage.map(st)
dd=d.dropna(subset=['stage_num'])
r,p=stats.spearmanr(dd.CAF_Signature_Score,dd.stage_num); print('stage spearman rho=%.3f p=%.4f n=%d'%(r,p,len(dd)))
print(dd.groupby('stage_num').CAF_Signature_Score.agg(['count','median','mean']).round(3))
kw=stats.kruskal(*[g.CAF_Signature_Score for _,g in dd.groupby('stage_num')]); print('KW p=%.4f'%kw.pvalue)
r,p=stats.spearmanr(d.CAF_Signature_Score,d.age_at_diagnosis,nan_policy='omit'); print('age rho=%.3f p=%.3g'%(r,p))
# Cox
c=d.dropna(subset=['stage_num','age_at_diagnosis']).copy()
c['score_per_SD']=(c.CAF_Signature_Score-c.CAF_Signature_Score.mean())/c.CAF_Signature_Score.std()
c['stage_III_IV']=(c.stage_num>=3).astype(int); c['age_per10']=c.age_at_diagnosis/10
out=[]
for name,cols in [('univariable',['score_per_SD']),('adjusted',['score_per_SD','age_per10','stage_III_IV'])]:
    cph=CoxPHFitter().fit(c[['os_months','os_event']+cols],'os_months','os_event')
    s=cph.summary
    for v in cols:
        out.append([name,v,s.loc[v,'exp(coef)'],s.loc[v,'exp(coef) lower 95%'],s.loc[v,'exp(coef) upper 95%'],s.loc[v,'p']])
    print(name,'n=',len(c),'events',int(c.os_event.sum()),'C-index %.3f'%cph.concordance_index_)
    # PH test
    from lifelines.statistics import proportional_hazard_test
    print(proportional_hazard_test(cph,c[['os_months','os_event']+cols],time_transform='rank').summary[['p']].round(3).to_dict())
o=pd.DataFrame(out,columns=['model','variable','HR','lo95','hi95','p']); print(o.round(4)); o.to_csv('/home/claude/re/out/tcga_cox.csv',index=False)
