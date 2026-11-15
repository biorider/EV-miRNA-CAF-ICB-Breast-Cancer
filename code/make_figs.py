import pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, json
from matplotlib.lines import Line2D
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
  'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':600,'pdf.fonttype':42,'svg.fonttype':'none'})
D='/home/claude/data/'; O='/home/claude/re/out/'; F='/home/claude/re/figs/'
OK=dict(blue='#0072B2',orange='#E69F00',green='#009E73',verm='#D55E00',purple='#CC79A7',sky='#56B4E9',yellow='#F0E442',grey='#7f7f7f',black='#000000')
CAT={'CAF Activation':OK['verm'],'Glycolysis':OK['orange'],'miRNA Target':OK['green'],'Immune Ligand':OK['purple'],'NF-κB / Signaling':OK['blue'],
     'OXPHOS':OK['sky'],'EMT':OK['grey'],'Fatty Acid Oxidation':OK['yellow'],'Lipid Synthesis':'#999999','Pentose Phosphate':'#bbbbbb'}
def lab(ax,s,x=-0.12,y=1.03): ax.text(x,y,s,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')
def save(fig,name):
    for ext in ['png','tiff','pdf']:
        kw=dict(bbox_inches='tight'); 
        if ext=='tiff': kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(F+f'{name}.{ext}',**kw)
    plt.close(fig)
# ---------------- Figure 1
net=pd.read_csv(D+'01_ev_mirna_targets/table2_mirna_target_network.csv')
de=pd.read_csv(D+'02_scrnaseq_caf/table4_caf_de_results.csv')
cnt=net.groupby('Target_Gene').miRNA.nunique(); hubs=cnt[cnt>=2].sort_values(ascending=True).index.tolist()
mirs=sorted(net.miRNA.unique(),key=lambda m:-net[net.miRNA==m].shape[0])
fig=plt.figure(figsize=(7.2,4.6)); gs=fig.add_gridspec(1,2,width_ratios=[1.25,1],wspace=0.6)
ax=fig.add_subplot(gs[0]); ecol={'miRTarBase':OK['blue'],'miRDB':OK['orange'],'miRTarBase+miRDB':OK['verm']}
for _,r in net[net.Target_Gene.isin(hubs)].iterrows():
    ax.scatter(mirs.index(r.miRNA),hubs.index(r.Target_Gene),s=38,c=ecol[r.Evidence],edgecolors='k',linewidths=0.3,zorder=3)
ax.set_xticks(range(len(mirs))); ax.set_xticklabels(mirs,rotation=60,ha='right')
ax.set_yticks(range(len(hubs))); ax.set_yticklabels([f'{g} ({cnt[g]})' for g in hubs])
ax.grid(color='#e5e5e5',lw=0.5); ax.set_axisbelow(True); ax.set_xlim(-0.6,len(mirs)-0.4)
ax.set_title('Hub targets (≥2 EV miRNAs)',loc='left')
ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc=c,mec='k',mew=0.3,label=k.replace('+',' + ')) for k,c in ecol.items()],
          loc='upper center',bbox_to_anchor=(0.5,-0.32),ncol=3,frameon=False,title='Evidence source')
lab(ax,'A',-0.28)
ax=fig.add_subplot(gs[1]); s=de[de.Significant].sort_values('log2FC')
y=np.arange(len(s)); cols=[CAT.get(c,OK['grey']) for c in s.Category]
ax.hlines(y,0,s.log2FC,color=cols,lw=1.6); ax.scatter(s.log2FC,y,c=cols,s=28,edgecolors='k',linewidths=0.3,zorder=3)
ax.axvline(0,color='k',lw=0.6); [ax.axvline(v,color='#bbbbbb',ls='--',lw=0.5) for v in (-0.25,0.25)]
ax.set_yticks(y); ax.set_yticklabels(s.Gene); ax.set_xlabel('log$_2$FC (CAF vs normal fibroblast)')
ax.set_title('Significant DE genes (18 of 59)',loc='left')
cats=list(dict.fromkeys(s.Category))
ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc=CAT[c],mec='k',mew=0.3,label=c) for c in cats],loc='lower right',frameon=False)
lab(ax,'B',-0.25); save(fig,'Figure1')
# ---------------- Figure 2
sc=pd.read_csv(D+'02_scrnaseq_caf/table5_caf_metabolic_scores.csv'); es=pd.read_csv(O+'sc_pathway_scores.csv').set_index('score')
fig=plt.figure(figsize=(7.4,3.6)); gs=fig.add_gridspec(1,2,width_ratios=[1.6,1],wspace=0.32)
ax=fig.add_subplot(gs[0]); sc_cols=['Glycolysis','OXPHOS','FAO','CAF_Activation','Warburg_Index']
pos=[];data=[];colors=[]
for i,c in enumerate(sc_cols):
    for j,(g,col) in enumerate([('Normal_Fibroblast',OK['blue']),('Cancer_CAF',OK['verm'])]):
        data.append(sc[sc.group==g][c].values); pos.append(i*2.6+j*0.95); colors.append(col)
vp=ax.violinplot(data,positions=pos,widths=0.85,showextrema=False)
for b,c in zip(vp['bodies'],colors): b.set_facecolor(c); b.set_alpha(0.45); b.set_edgecolor('none')
bp=ax.boxplot(data,positions=pos,widths=0.22,showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
for b in bp['boxes']: b.set_facecolor('white'); b.set_linewidth(0.6)
for i,c in enumerate(sc_cols):
    m=sc.groupby('group')[c].mean(); ax.scatter([i*2.6,i*2.6+0.95],[m['Normal_Fibroblast'],m['Cancer_CAF']],marker='D',s=10,c='k',zorder=4)
    ax.text(i*2.6+0.47,3.55,f"r$_{{rb}}$={es.loc[c,'rank_biserial']:+.2f}",ha='center',fontsize=6.5)
ax.set_xticks([i*2.6+0.47 for i in range(len(sc_cols))]); ax.set_xticklabels(['Glycolysis','OXPHOS','FAO','CAF\nactivation','Warburg\nindex'])
ax.set_ylabel('Per-cell score (mean log1p expression)'); ax.set_ylim(-1.2,3.9)
ax.legend(handles=[Line2D([],[],marker='s',ls='',mfc=OK['blue'],alpha=.6,label='Normal fibroblast'),Line2D([],[],marker='s',ls='',mfc=OK['verm'],alpha=.6,label='CAF'),Line2D([],[],marker='D',ls='',mfc='k',ms=3,label='mean')],loc='upper center',bbox_to_anchor=(0.5,-0.14),frameon=False,ncol=3)
ax.set_title('Per-cell pathway scores',loc='left'); lab(ax,'A',-0.1)
ax=fig.add_subplot(gs[1])
for _,r in de.iterrows():
    c=CAT.get(r.Category,OK['grey']); ax.scatter(r.pct_normal,r.pct_cancer,s=22 if r.Significant else 10,c=c,edgecolors='k' if r.Significant else 'none',linewidths=0.3,alpha=0.9)
for g in ['CXCL1','IL6','LDHA','MYC','VEGFA','NFKB1','PTEN','COL1A1','ACTA2','POSTN','FOXO3','SIRT1']:
    r=de[de.Gene==g].iloc[0]; off={'PTEN':(-22,-8),'NFKB1':(2,4),'VEGFA':(3,-8),'FOXO3':(-26,4)}.get(g,(3,-2)); ax.annotate(g,(r.pct_normal,r.pct_cancer),xytext=off,textcoords='offset points',fontsize=6)
ax.plot([0,100],[0,100],'k--',lw=0.5); ax.set_xlim(0,100); ax.set_ylim(0,100)
ax.set_xlabel('% normal fibroblasts detecting gene'); ax.set_ylabel('% CAFs detecting gene'); ax.set_title('Detection rates (59-gene panel)',loc='left'); lab(ax,'B',-0.25)
save(fig,'Figure2')
# ---------------- Figure 3
cov=pd.read_csv(O+'recon3d_pathway_coverage.csv'); lo=pd.read_csv(O+'warburg_leave_one_out.csv'); am=pd.read_csv(O+'antagomir_sensitivity.csv')
g=json.load(open(O+'gpr_summary.json'))
fig=plt.figure(figsize=(7.2,3.3)); gs=fig.add_gridspec(1,3,width_ratios=[1.25,1,1.05],wspace=0.55)
ax=fig.add_subplot(gs[0]); c2=cov.dropna(subset=['CAF_score']).sort_values('CAF_score')
ax.barh(range(len(c2)),c2.CAF_score,color=[OK['verm'] if v>0 else OK['blue'] for v in c2.CAF_score],height=0.6)
ax.set_yticks(range(len(c2))); ax.set_yticklabels([f"{p}\n(n={int(n)} rxn; {len(str(gg).split(','))} genes)" for p,n,gg in zip(c2.Pathway,c2.Reactions_scored,c2.Measured_genes_contributing)],fontsize=6)
ax.axvline(0,color='k',lw=0.6); ax.set_xlabel('GPR-mapped activity score\n(mean log$_2$FC, CAF vs NF)'); ax.set_title('Recon3D GPR mapping',loc='left'); lab(ax,'A',-0.75)
ax=fig.add_subplot(gs[1]); l2=lo.sort_values('WI'); show=pd.concat([l2.head(4),l2.tail(4)])
ax.barh(range(len(show)),show.WI,color=OK['grey'],height=0.6); ax.axvline(g['WI_tree'],color=OK['verm'],lw=1,ls='--',label=f"all genes ({g['WI_tree']:.3f})")
ax.set_yticks(range(len(show))); ax.set_yticklabels(['− '+x for x in show.left_out]); ax.set_xlabel('Warburg index\n(Glycolysis − OXPHOS)')
ax.set_title('Leave-one-gene-out',loc='left'); ax.legend(frameon=False,loc='upper left',fontsize=6); lab(ax,'B',-0.45)
ax=fig.add_subplot(gs[2]); pv=am.pivot(index='beta',columns='eps',values='rescue_pct').sort_index(ascending=False)
im=ax.imshow(pv.values,cmap='viridis',vmin=0,vmax=100,aspect='auto')
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]): ax.text(j,i,f"{pv.values[i,j]:.0f}",ha='center',va='center',fontsize=6.5,color='w' if pv.values[i,j]<60 else 'k')
ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels([f"{x:.0%}" for x in pv.columns]); ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels([f"{x:g}" for x in pv.index])
ax.set_xlabel('AntagomiR efficiency (ε)'); ax.set_ylabel('Max repression per\nevidence unit, β (log$_2$)'); ax.set_title('Predicted rescue of\nWarburg index (%)',loc='left')
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.ax.tick_params(labelsize=6); lab(ax,'C',-0.5); save(fig,'Figure3')
# ---------------- Figure 4
t=pd.read_csv(D+'04_module5_signature/table9_tcga_brca_signature_scores.csv'); cox=pd.read_csv(O+'tcga_cox.csv')
def st(s):
    if not isinstance(s,str): return np.nan
    s=s.replace('Stage ','')
    for k,v in [('IV',4),('III',3),('II',2),('I',1)]:
        if s.startswith(k): return v
t['stage']=t.ajcc_stage.map(st)
fig,axs=plt.subplots(2,2,figsize=(7.2,5.6)); plt.subplots_adjust(wspace=0.35,hspace=0.5)
ax=axs[0,0]; ax.hist(t.CAF_Signature_Score,bins=40,color=OK['sky'],edgecolor='white',lw=0.3); ax.axvline(t.CAF_Signature_Score.median(),color='k',ls='--',lw=0.8)
ax.set_xlabel('34-gene composite score (mean z)'); ax.set_ylabel('Patients'); ax.set_title(f'TCGA-BRCA (n = {len(t)} with OS data)',loc='left'); lab(ax,'A')
ax=axs[0,1]; d=[t[t.stage==k].CAF_Signature_Score for k in [1,2,3,4]]
bp=ax.boxplot(d,showfliers=False,patch_artist=True,widths=0.55,medianprops=dict(color='k'))
for b in bp['boxes']: b.set_facecolor(OK['sky']); b.set_alpha(0.7)
for i,x in enumerate(d): ax.scatter(np.random.default_rng(i).normal(i+1,0.07,len(x)),x,s=2,c='k',alpha=0.25)
ax.set_xticklabels([f'{s}\n(n={len(x)})' for s,x in zip(['I','II','III','IV'],d)]); ax.set_xlabel('AJCC stage'); ax.set_ylabel('Score')
ax.set_title('Spearman ρ = 0.073, p = 0.025; Kruskal–Wallis p < 0.001',loc='left',fontsize=7.5); lab(ax,'B')
ax=axs[1,0]
for gname,col in [('High',OK['verm']),('Low',OK['blue'])]:
    s=t[t.Signature_Group==gname]; KaplanMeierFitter(label=f'{gname} (n={len(s)}, events={int(s.os_event.sum())})').fit(s.os_months,s.os_event).plot_survival_function(ax=ax,ci_show=True,color=col,lw=1)
ax.set_xlabel('Months'); ax.set_ylabel('Overall survival probability'); ax.set_ylim(0,1.02); ax.legend(frameon=False,loc='lower left'); ax.set_title('Median split; log-rank p = 0.60',loc='left'); lab(ax,'C')
ax=axs[1,1]; lbl={'score_per_SD':'Score/SD','age_per10':'Age/10 y','stage_III_IV':'Stage III–IV'}
rows=[(f"{lbl[r.variable]} ({'uni' if r.model=='univariable' else 'adj'})",r.HR,r.lo95,r.hi95,r.p) for r in cox.itertuples()]
for i,(n,h,l,u,p) in enumerate(rows[::-1]):
    ax.errorbar(h,i,xerr=[[h-l],[u-h]],fmt='s',color=OK['verm'] if 'Signature' in n else 'k',ms=4,capsize=2,lw=0.8)
    ax.text(5.1,i,f'{h:.2f} ({l:.2f}–{u:.2f})\np={p:.2g}',va='center',fontsize=6)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]],fontsize=6.5); ax.axvline(1,color='k',lw=0.6,ls='--')
ax.set_xlim(0.5,5); ax.set_xticks([1,2,3,4]); ax.set_xlabel('Hazard ratio (95% CI)')
ax.set_title('Cox models (n = 955, 102 events)',loc='left'); lab(ax,'D'); save(fig,'Figure4')
# ---------------- Figure 5
m=pd.read_csv(O+'icb_pergene_merged.csv')
fig=plt.figure(figsize=(7.4,4.2)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.05],wspace=0.42)
ax=fig.add_subplot(gs[0]); texts=[]
for r in m.itertuples():
    col=OK['verm'] if r.replicated_FDR05 else (OK['orange'] if r.replicated_FDR10 else '#bbbbbb')
    ax.scatter(r.rho_pembro,r.rho_durva,s=26,c=col,marker='o' if r.in_34gene_signature else '^',edgecolors='k',linewidths=0.3,zorder=3)
    if r.replicated_FDR10 or r.gene in ('PTEN','ACTA2','COL1A1','SIRT1','TRAF6','FOXO3'):
        texts.append(ax.text(r.rho_pembro,r.rho_durva,r.gene,fontsize=6))
ax.axhline(0,color='k',lw=0.5); ax.axvline(0,color='k',lw=0.5); ax.plot([-0.6,0.6],[-0.6,0.6],ls=':',c='grey',lw=0.6)
ax.set_xlim(-0.6,0.6); ax.set_ylim(-0.6,0.65); from adjustText import adjust_text; adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',lw=0.3,color='grey')); ax.set_xlabel('Spearman ρ with pCR — pembrolizumab arm\n(GSE194040, n = 69)'); ax.set_ylabel('Spearman ρ with pCR — durvalumab/olaparib arm\n(GSE173839, n = 71)')
ax.set_title('Cross-cohort concordance (ρ$_s$ = 0.83)',loc='left')
ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc=OK['verm'],mec='k',mew=.3,label='FDR<5% both'),Line2D([],[],marker='o',ls='',mfc=OK['orange'],mec='k',mew=.3,label='FDR<10% both'),
   Line2D([],[],marker='o',ls='',mfc='#bbbbbb',mec='k',mew=.3,label='not replicated'),Line2D([],[],marker='o',ls='',mfc='w',mec='k',label='in 34-gene signature'),Line2D([],[],marker='^',ls='',mfc='w',mec='k',label='panel only')],
   loc='upper left',frameon=False,fontsize=6); lab(ax,'A',-0.2)
ax=fig.add_subplot(gs[1]); mm=m[(m.padj_pembro<0.25)|(m.padj_durva<0.25)].copy(); mm['avg']=(mm.rho_pembro+mm.rho_durva)/2; mm=mm.sort_values('avg')
for j,(col,pc,name) in enumerate([('rho_pembro','padj_pembro','Pembrolizumab'),('rho_durva','padj_durva','Durvalumab/olaparib')]):
    for i,r in enumerate(mm.itertuples()):
        rho=getattr(r,col); q=getattr(r,pc)
        ax.scatter(j,i,s=max(8,-np.log10(q)*45),c=[plt.cm.RdBu_r((rho+0.6)/1.2)],edgecolors='k' if q<0.05 else 'none',linewidths=0.6)
ax.set_xticks([0,1]); ax.set_xticklabels(['Pembrolizumab\n(GSE194040)','Durvalumab/olaparib\n(GSE173839)']); ax.set_xlim(-0.6,1.6)
ax.set_yticks(range(len(mm))); ax.set_yticklabels([g+('' if s else ' †') for g,s in zip(mm.gene,mm.in_34gene_signature)])
sm=plt.cm.ScalarMappable(cmap='RdBu_r',norm=plt.Normalize(-0.6,0.6)); cb=fig.colorbar(sm,ax=ax,fraction=0.05,pad=0.02); cb.set_label('Spearman ρ with pCR'); cb.ax.tick_params(labelsize=6)
ax.legend(handles=[plt.scatter([],[],s=-np.log10(q)*45,c='grey',label=f'FDR {q}') for q in (0.2,0.05,0.001)],loc='upper center',bbox_to_anchor=(0.45,-0.16),ncol=3,frameon=False,fontsize=6)
ax.set_title('Genes with FDR<25% in ≥1 cohort',loc='left'); lab(ax,'B',-0.25); save(fig,'Figure5')
print('done')
