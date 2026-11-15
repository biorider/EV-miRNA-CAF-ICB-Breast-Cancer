exec(open('/home/claude/re/scripts/common.py').read())
de=pd.read_csv(D+'02_scrnaseq_caf/table4_caf_de_results.csv')
# ---------------- Figure 2 (R1b: CAF activation on its own axis; effect sizes moved into tick labels)
sc=pd.read_csv(D+'02_scrnaseq_caf/table5_caf_metabolic_scores.csv'); es=pd.read_csv(O+'sc_pathway_scores.csv').set_index('score')
GR=[('Normal_Fibroblast',OK['blue']),('Cancer_CAF',OK['verm'])]
NAMES={'Glycolysis':'Glycolysis','OXPHOS':'OXPHOS','FAO':'FAO','Warburg_Index':'Warburg\nindex','CAF_Activation':'CAF\nactivation'}
def violins(ax,cols,step=2.9):
    pos=[];data=[];colors=[]
    for i,c in enumerate(cols):
        for j,(g,col) in enumerate(GR):
            data.append(sc[sc.group==g][c].values); pos.append(i*step+j*0.95); colors.append(col)
    vp=ax.violinplot(data,positions=pos,widths=0.85,showextrema=False)
    for b,c in zip(vp['bodies'],colors): b.set_facecolor(c); b.set_alpha(0.45); b.set_edgecolor('none')
    bp=ax.boxplot(data,positions=pos,widths=0.22,showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.8),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
    for b in bp['boxes']: b.set_facecolor('white'); b.set_linewidth(0.6)
    for i,c in enumerate(cols):
        m=sc.groupby('group')[c].mean(); ax.scatter([i*step,i*step+0.95],[m['Normal_Fibroblast'],m['Cancer_CAF']],marker='D',s=10,c='k',zorder=4)
    ax.set_xticks([i*step+0.47 for i in range(len(cols))])
    ax.set_xticklabels([NAMES[c].replace('\n',' ')+'\n' + f"r$_{{rb}}$ = {es.loc[c,'rank_biserial']:+.2f}".replace('-','−') for c in cols],fontsize=6.5)
    ax.set_xlim(-0.7,(len(cols)-1)*step+1.65)
fig=plt.figure(figsize=(7.8,3.7)); gs=fig.add_gridspec(1,3,width_ratios=[1.5,0.42,1.05],wspace=0.5)
ax=fig.add_subplot(gs[0]); violins(ax,['Glycolysis','OXPHOS','FAO','Warburg_Index'])
ax.axhline(0,color='#cccccc',lw=0.5,zorder=0)
ax.set_ylabel('Per-cell score (mean log1p expression)'); ax.set_ylim(-1.15,2.8)
ax.set_title('Per-cell metabolic pathway scores',loc='left'); lab(ax,'A',-0.13)
ax2=fig.add_subplot(gs[1]); violins(ax2,['CAF_Activation']); ax2.set_ylim(-0.1,5.5); ax2.set_ylabel('Per-cell score',labelpad=1)
ax2.set_title('CAF activation',loc='center')
ax.legend(handles=[Line2D([],[],marker='s',ls='',mfc=OK['blue'],alpha=.6,label='Normal fibroblast (n = 8,000)'),Line2D([],[],marker='s',ls='',mfc=OK['verm'],alpha=.6,label='CAF (n = 8,000)'),Line2D([],[],marker='D',ls='',mfc='k',ms=3,label='mean')],loc='upper center',bbox_to_anchor=(0.75,-0.2),frameon=False,ncol=3)
ax=fig.add_subplot(gs[2])
for _,r in de.iterrows():
    c=CAT.get(r.Category,OK['grey']); ax.scatter(r.pct_normal,r.pct_cancer,s=22 if r.Significant else 10,c=c,edgecolors='k' if r.Significant else 'none',linewidths=0.3,alpha=0.9)
for g in ['CXCL1','IL6','LDHA','MYC','VEGFA','NFKB1','PTEN','COL1A1','ACTA2','POSTN','FOXO3','SIRT1']:
    r=de[de.Gene==g].iloc[0]; off={'PTEN':(-22,-8),'NFKB1':(2,4),'VEGFA':(3,-8),'FOXO3':(-26,4)}.get(g,(3,-2)); ax.annotate(g,(r.pct_normal,r.pct_cancer),xytext=off,textcoords='offset points',fontsize=6)
ax.plot([0,100],[0,100],'k--',lw=0.5); ax.set_xlim(0,100); ax.set_ylim(0,100)
ax.set_xlabel('% normal fibroblasts detecting gene'); ax.set_ylabel('% CAFs detecting gene'); ax.set_title('Detection rates (59-gene panel)',loc='left'); lab(ax,'B',-0.34)
save(fig,'Figure2')
# ---------------- Figure 3 (R1b: reference line labelled below the axis, not over bars)
cov=pd.read_csv(O+'recon3d_pathway_coverage.csv'); lo=pd.read_csv(O+'warburg_leave_one_out.csv'); am=pd.read_csv(O+'antagomir_sensitivity.csv')
g=json.load(open(O+'gpr_summary.json'))
fig=plt.figure(figsize=(7.2,3.5)); gs=fig.add_gridspec(1,3,width_ratios=[1.25,1,1.05],wspace=0.55)
ax=fig.add_subplot(gs[0]); c2=cov.dropna(subset=['CAF_score']).sort_values('CAF_score')
ax.barh(range(len(c2)),c2.CAF_score,color=[OK['verm'] if v>0 else OK['blue'] for v in c2.CAF_score],height=0.6)
ax.set_yticks(range(len(c2))); ax.set_yticklabels([f"{p}\n(n={int(n)} rxn; {len(str(gg).split(','))} gene{'s' if len(str(gg).split(','))>1 else ''})" for p,n,gg in zip(c2.Pathway,c2.Reactions_scored,c2.Measured_genes_contributing)],fontsize=6)
ax.axvline(0,color='k',lw=0.6); ax.set_xlabel('GPR-mapped activity score\n(mean log$_2$FC, CAF vs NF)'); ax.set_title('Recon3D GPR mapping',loc='left'); lab(ax,'A',-0.75)
ax=fig.add_subplot(gs[1]); l2=lo.sort_values('WI'); show=pd.concat([l2.head(4),l2.tail(4)])
ax.barh(range(len(show)),show.WI,color=OK['grey'],height=0.6); ax.axvline(g['WI_tree'],color=OK['verm'],lw=1,ls='--')
ax.set_yticks(range(len(show))); ax.set_yticklabels(['− '+x for x in show.left_out]); ax.set_xlabel('Warburg index\n(Glycolysis − OXPHOS)')
ax.set_title('Leave-one-gene-out',loc='left'); lab(ax,'B',-0.45)
ax.legend(handles=[Line2D([],[],color=OK['verm'],ls='--',lw=1,label=f"All genes retained (WI = {g['WI_tree']:.3f})".replace('-','−'))],
          loc='upper center',bbox_to_anchor=(0.5,-0.2),frameon=False,fontsize=6.5)
ax=fig.add_subplot(gs[2]); pv=am.pivot(index='beta',columns='eps',values='rescue_pct').sort_index(ascending=False)
im=ax.imshow(pv.values,cmap='viridis',vmin=0,vmax=100,aspect='auto')
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]): ax.text(j,i,f"{pv.values[i,j]:.0f}",ha='center',va='center',fontsize=6.5,color='w' if pv.values[i,j]<60 else 'k')
ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels([f"{x:.0%}" for x in pv.columns]); ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels([f"{x:g}" for x in pv.index])
ax.set_xlabel('AntagomiR efficiency (ε)'); ax.set_ylabel('Max repression per\nevidence unit, β (log$_2$)'); ax.set_title('Predicted rescue of\nWarburg index (%)',loc='left')
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.03); cb.ax.tick_params(labelsize=6); lab(ax,'C',-0.5); save(fig,'Figure3')
print('ok')
