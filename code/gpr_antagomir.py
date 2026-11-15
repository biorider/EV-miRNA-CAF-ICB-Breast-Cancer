import cobra, pandas as pd, numpy as np, re, ast, itertools, json, sys, warnings
warnings.filterwarnings('ignore')
D='/home/claude/data/'; O='/home/claude/re/out/'
de=pd.read_csv(D+'02_scrnaseq_caf/table4_caf_de_results.csv').set_index('Gene')
lfc=de.log2FC.to_dict()
net=pd.read_csv(D+'01_ev_mirna_targets/table2_mirna_target_network.csv')
gm=pd.read_csv('/home/claude/ext/HumanGEM_genes.tsv',sep='\t',dtype=str)
ent2sym={}
for _,r in gm.iterrows():
    if isinstance(r.geneEntrezID,str):
        for e in r.geneEntrezID.split(';'): ent2sym[e.strip()]=r.geneSymbols
m=cobra.io.load_matlab_model('/home/claude/ext/Recon3D_301.mat')
def sym(gid): return ent2sym.get(gid.split('.')[0])
# --- GPR evaluation on measured genes only
def tok(rule):
    return re.findall(r'\(|\)|and|or|[^\s()]+',rule)
def parse(tokens):
    def expr(i):
        node,i=term(i); items=[node]; op=None
        while i<len(tokens) and tokens[i] in('and','or'):
            o=tokens[i]; n,i=term(i+1)
            if op is None: op=o
            if o!=op:  # mixed precedence: and binds tighter -> rebuild
                op='mixed'
            items.append((o,n))
        if len(items)==1: return node,i
        # build with and-precedence
        groups=[[items[0]]]
        for o,n in items[1:]:
            if o=='and': groups[-1].append(n)
            else: groups.append([n])
        ors=[('and',g) if len(g)>1 else g[0] for g in groups]
        return (('or',ors) if len(ors)>1 else ors[0]),i
    def term(i):
        t=tokens[i]
        if t=='(':
            n,i=expr(i+1); return n,i+1
        return ('g',t),i+1
    n,_=expr(0); return n
def ev(node,score):
    if node[0]=='g':
        s=sym(node[1]); return score.get(s) if s else None
    vals=[ev(c,score) for c in node[1]]; vals=[v for v in vals if v is not None]
    if not vals: return None
    return min(vals) if node[0]=='and' else max(vals)
def flat(node):
    if node[0]=='g': return [node[1]]
    return list(itertools.chain(*[flat(c) for c in node[1]]))
def ev_mean(node,score):  # original 'mixed=mean' heuristic
    gs=[sym(g) for g in flat(node)]; v=[score[g] for g in gs if g in score]
    return np.mean(v) if v else None
PW={'Glycolysis':['Glycolysis/gluconeogenesis'],'OXPHOS':['Oxidative phosphorylation'],'TCA cycle':['Citric acid cycle'],
    'Fatty acid oxidation':['Fatty acid oxidation'],'Fatty acid synthesis':['Fatty acid synthesis'],
    'Pentose phosphate':['Pentose phosphate pathway'],'Glutamate/glutaminolysis':['Glutamate metabolism']}
rx={}
for r in m.reactions:
    if not r.gene_reaction_rule: continue
    rx[r.id]=(r.subsystem,parse(tok(r.gene_reaction_rule)),r)
lac=[rid for rid,(s,n,r) in rx.items() if any(mm.id=='lac_L[e]' for mm in r.metabolites)]
def pathway_scores(score,mode='tree'):
    res={}; cov={}
    groups={k:[rid for rid,(s,n,r) in rx.items() if s in v] for k,v in PW.items()}; groups['Lactate transport (MCT)']=lac
    for k,rids in groups.items():
        vals=[];genes=set()
        for rid in rids:
            n=rx[rid][1]; v=ev(n,score) if mode=='tree' else ev_mean(n,score)
            if v is not None:
                vals.append(v); genes|={sym(g) for g in flat(n) if sym(g) in score}
        res[k]=np.mean(vals) if vals else np.nan; cov[k]=(len(rids),len(vals),','.join(sorted(genes)))
    return res,cov
base,cov=pathway_scores(lfc)
cv=pd.DataFrame([(k,a,b,c) for k,(a,b,c) in cov.items()],columns=['Pathway','Reactions_with_GPR','Reactions_scored','Measured_genes_contributing'])
cv['CAF_score']=cv.Pathway.map(base)
basem,_=pathway_scores(lfc,'mean'); cv['CAF_score_meanrule']=cv.Pathway.map(basem)
print(cv.to_string()); cv.to_csv(O+'recon3d_pathway_coverage.csv',index=False)
WI=lambda s: s['Glycolysis']-s['OXPHOS']
print('Warburg index (tree GPR) %.3f ; mean-rule %.3f'%(WI(base),WI(basem)))
# leave-one-gene-out
lo=[]
for g in lfc:
    sc={k:v for k,v in lfc.items() if k!=g}; s,_=pathway_scores(sc); lo.append((g,WI(s)))
lo=pd.DataFrame(lo,columns=['left_out','WI']).sort_values('WI'); lo.to_csv(O+'warburg_leave_one_out.csv',index=False)
print('LOO range',lo.WI.min().round(3),lo.WI.max().round(3)); print(lo.head(4).to_string(),lo.tail(4).to_string())
# --- miRNA-informed antagomiR model
w={'miRTarBase+miRDB':1.0,'miRTarBase':1.0,'miRDB':0.5}
net['w']=net.Evidence.map(w)
S=net.groupby('Target_Gene').w.sum().to_dict(); K=net.groupby('Target_Gene').miRNA.nunique().to_dict()
targets_in_panel=sorted(g for g in S if g in lfc)
print('targets measured',len(targets_in_panel),targets_in_panel)
rows=[]
for beta in [0.1,0.2,0.3,0.5,1.0]:
    for eps in [0.5,0.7,0.9,1.0]:
        sc=dict(lfc)
        for g in targets_in_panel:
            if lfc[g]<0:
                mcomp=min(-lfc[g],beta*S[g])  # miRNA-attributable repression (log2)
                sc[g]=lfc[g]+eps*mcomp
        s,_=pathway_scores(sc)
        rows.append(dict(beta=beta,eps=eps,Glycolysis=s['Glycolysis'],OXPHOS=s['OXPHOS'],WI=WI(s),
             rescue_pct=100*(WI(base)-WI(s))/WI(base) if WI(base)!=0 else np.nan))
am=pd.DataFrame(rows); am.to_csv(O+'antagomir_sensitivity.csv',index=False); print(am.round(3).to_string())
# per-gene rescued values at reference scenario
ref=[]
for g in targets_in_panel:
    mc=min(-lfc[g],0.3*S[g]) if lfc[g]<0 else 0
    ref.append(dict(gene=g,n_miRNAs=K[g],evidence_weight=S[g],log2FC_CAF=lfc[g],miRNA_attributable=-mc,rescued_log2FC_eps07=lfc[g]+0.7*mc))
pd.DataFrame(ref).to_csv(O+'antagomir_pergene_reference.csv',index=False); print(pd.DataFrame(ref).round(3).to_string())
# legacy uniform scaling for comparison
s,_=pathway_scores({k:0.4*v for k,v in lfc.items()}); print('legacy uniform 0.4 scaling WI %.3f'%WI(s))
json.dump(dict(WI_tree=WI(base),WI_mean=WI(basem),scores=base,scores_mean=basem),open(O+'gpr_summary.json','w'),indent=1,default=float)
