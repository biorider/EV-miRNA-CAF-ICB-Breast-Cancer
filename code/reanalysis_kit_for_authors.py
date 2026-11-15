"""
Reanalysis kit for the R1 revision (EV-miRNA / CAF / ICB manuscript)
=====================================================================
Purpose: run the four analyses that could NOT be re-executed in the revision
sandbox because raw GEO / GDC / CellxGene data were not reachable there.
Run this in the Biomni environment (or any Python 3.10+ env with internet access).

    pip install pandas numpy scipy scikit-learn lifelines statsmodels GEOparse cellxgene-census scanpy

Outputs are written to ./R1_reanalysis_outputs/. Each block is independent; set the
RUN_* flags below. Every block prints the exact sentences/values needed to update the
manuscript (search the revised manuscript for "[R1-KIT]" notes in the response letter).

A. I-SPY2 re-scoring with the EXACT 34-gene signature (Table 5) instead of the legacy
   32-gene panel, plus (i) a CAF-restricted sub-score, (ii) adjustment for a stromal-
   abundance proxy, (iii) HR/HER2 stratification, (iv) logistic models with arm x score
   interaction (formal test of ICB specificity), and (v) bootstrap 95% CIs for AUC.
B. TCGA-BRCA: per-gene Cox models (all patients, censoring respected) replacing the
   deceased-only Spearman analysis; stromal/purity-adjusted Cox; PAM50 stratification.
C. scRNA-seq: donor-level pseudobulk differential expression (edgeR-like / Wilcoxon on
   donor means) to remove cell-level pseudoreplication; dataset/assay-stratified check;
   CAF subtype (myCAF/iCAF) composition of the cancer vs normal fibroblast pools.
"""
import os, numpy as np, pandas as pd
from scipy import stats
OUT = "R1_reanalysis_outputs"; os.makedirs(OUT, exist_ok=True)
RUN_A, RUN_B, RUN_C = True, True, True

SIG34 = ["ACTA2","ALDOA","BAK1","CDKN1B","COL1A1","COL1A2","CXCL1","CXCL12","ENO1","FAP","FN1",
         "FOXO3","GLS","GPD1L","HOXD10","IL6","IRAK1","ISCU","LATS2","LDHA","MYC","NFKB1","PDCD4",
         "PKM","POSTN","PTEN","S100A4","SIRT1","STAT3","TGFB2","TRAF6","VEGFA","ZEB1","ZEB2"]
CAF_ONLY = ["ACTA2","COL1A1","COL1A2","FAP","POSTN","S100A4","FN1"]
# Stromal-abundance proxy: canonical pan-fibroblast genes NOT in the signature
STROMA_PROXY = ["DCN","LUM","PDGFRA","PDGFRB","COL3A1","COL6A3","SFRP2","MMP2"]

def zmean(expr, genes):
    g = [x for x in genes if x in expr.index]
    z = expr.loc[g].sub(expr.loc[g].mean(1), axis=0).div(expr.loc[g].std(1) + 1e-9, axis=0)
    return z.mean(0), g

def boot_auc(y, s, n=2000, seed=1):
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed); a = []
    y = np.asarray(y); s = np.asarray(s)
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) == 2: a.append(roc_auc_score(y[i], s[i]))
    return roc_auc_score(y, s), np.percentile(a, [2.5, 97.5])

# ------------------------------------------------------------------------------- A
if RUN_A:
    import GEOparse, statsmodels.formula.api as smf
    # NOTE: GSE194040 expression is a supplementary matrix (ComBat-adjusted, gene symbols).
    # Download it manually from GEO if GEOparse does not expose it, then point EXPR_FILE to it.
    COHORTS = {
        "GSE194040": dict(EXPR_FILE=None, icb_arm="Pembrolizumab", ctrl_arm="Paclitaxel"),
        "GSE173839": dict(EXPR_FILE=None, icb_arm="Durvalumab/Olaparib", ctrl_arm="Control"),
    }
    rows = []
    for gse, cfg in COHORTS.items():
        g = GEOparse.get_GEO(gse, destdir=OUT, silent=True)
        meta = g.phenotype_data
        # --- harmonise clinical columns (inspect meta.columns once and adapt these 4 lines) ---
        arm  = meta.filter(regex="arm|treatment", axis=1).iloc[:, 0].astype(str)
        pcr  = meta.filter(regex="pcr", axis=1).iloc[:, 0].astype(str).str.extract(r"(\d)")[0].astype(float)
        hr   = meta.filter(regex="hr", axis=1).iloc[:, 0].astype(str).str.extract(r"(\d)")[0].astype(float)
        her2 = meta.filter(regex="her2", axis=1).iloc[:, 0].astype(str).str.extract(r"(\d)")[0].astype(float)
        expr = pd.read_csv(cfg["EXPR_FILE"], sep="\t", index_col=0) if cfg["EXPR_FILE"] else g.pivot_samples("VALUE")
        expr = expr.groupby(level=0).mean()
        df = pd.DataFrame({"arm": arm, "pCR": pcr, "HR": hr, "HER2": her2})
        df["sig34"], used = zmean(expr, SIG34)
        df["caf_only"], _ = zmean(expr, CAF_ONLY)
        df["stroma"], _ = zmean(expr, STROMA_PROXY)
        df["sig34_resid"] = smf.ols("sig34 ~ stroma", df.dropna()).fit().resid
        df["ICB"] = df.arm.str.contains(cfg["icb_arm"], case=False).astype(int)
        df = df.dropna(subset=["pCR"])
        df.to_csv(f"{OUT}/{gse}_scores.csv")
        print(gse, "signature genes found:", len(used), "/34")
        for armname, sub in [("ICB", df[df.ICB == 1]), ("Control", df[df.ICB == 0])]:
            for score in ["sig34", "caf_only", "sig34_resid"]:
                a, b = sub[sub.pCR == 1][score], sub[sub.pCR == 0][score]
                p = stats.mannwhitneyu(a, b).pvalue
                auc, ci = boot_auc(sub.pCR, -sub[score])
                rows.append([gse, armname, score, len(a), len(b), p, auc, *ci])
            for lab, s2 in [("TNBC", sub[(sub.HR == 0) & (sub.HER2 == 0)]), ("HR+HER2-", sub[(sub.HR == 1) & (sub.HER2 == 0)])]:
                if s2.pCR.nunique() == 2 and len(s2) > 10:
                    auc, ci = boot_auc(s2.pCR, -s2.sig34)
                    rows.append([gse, f"{armname}:{lab}", "sig34", int(s2.pCR.sum()), int((1 - s2.pCR).sum()),
                                 stats.mannwhitneyu(s2[s2.pCR == 1].sig34, s2[s2.pCR == 0].sig34).pvalue, auc, *ci])
        # Formal ICB-specificity test (interaction), adjusted for HR status
        m = smf.logit("pCR ~ sig34 * ICB + HR", df.dropna(subset=["HR"])).fit(disp=0)
        print(gse, "interaction sig34:ICB  OR=%.2f  p=%.3f" % (np.exp(m.params["sig34:ICB"]), m.pvalues["sig34:ICB"]))
        # Per-gene Spearman (ICB arm) for the 34 genes, BH-FDR
        from statsmodels.stats.multitest import multipletests
        icb = df[df.ICB == 1]; pg = []
        for gname in used:
            r, p = stats.spearmanr(expr.loc[gname, icb.index], icb.pCR); pg.append([gname, r, p])
        pg = pd.DataFrame(pg, columns=["gene", "rho", "p"]); pg["padj"] = multipletests(pg.p, method="fdr_bh")[1]
        pg.sort_values("p").to_csv(f"{OUT}/{gse}_pergene_sig34.csv", index=False)
    pd.DataFrame(rows, columns=["cohort", "arm", "score", "n_pCR", "n_nonpCR", "wilcoxon_p", "AUC", "AUC_lo", "AUC_hi"]).to_csv(f"{OUT}/ISPY2_summary.csv", index=False)

# ------------------------------------------------------------------------------- B
if RUN_B:
    from lifelines import CoxPHFitter
    from statsmodels.stats.multitest import multipletests
    # Expect: TCGA-BRCA log2(TPM+1) gene x sample matrix + clinical table with os_months, os_event,
    # age, stage, PAM50, and an ESTIMATE stromal score / tumour purity (e.g. from TCGA PanCanAtlas).
    expr = pd.read_csv("TCGA_BRCA_log2TPM.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv("TCGA_BRCA_clinical_with_purity_PAM50.tsv", sep="\t", index_col=0)
    common = clin.index.intersection(expr.columns); expr, clin = expr[common], clin.loc[common]
    clin["sig34"], _ = zmean(expr, SIG34)
    res = []
    for gname in [x for x in SIG34 if x in expr.index]:
        d = clin[["os_months", "os_event", "age", "stage34", "purity"]].copy()
        d["gene"] = (expr.loc[gname] - expr.loc[gname].mean()) / expr.loc[gname].std()
        d = d.dropna()
        cu = CoxPHFitter().fit(d[["os_months", "os_event", "gene"]], "os_months", "os_event")
        ca = CoxPHFitter().fit(d, "os_months", "os_event")
        res.append([gname, cu.hazard_ratios_["gene"], cu.summary.loc["gene", "p"], ca.hazard_ratios_["gene"], ca.summary.loc["gene", "p"]])
    res = pd.DataFrame(res, columns=["gene", "HR_uni", "p_uni", "HR_adj_age_stage_purity", "p_adj"])
    res["padj_uni"] = multipletests(res.p_uni, method="fdr_bh")[1]; res["padj_adj"] = multipletests(res.p_adj, method="fdr_bh")[1]
    res.sort_values("p_uni").to_csv(f"{OUT}/TCGA_pergene_cox.csv", index=False)
    r, p = stats.spearmanr(clin.sig34, clin.purity, nan_policy="omit")
    print("TCGA signature vs tumour purity: rho=%.2f p=%.2g" % (r, p))
    for sub, d in clin.groupby("PAM50"):
        d = d[["os_months", "os_event", "sig34"]].dropna()
        if d.os_event.sum() >= 10:
            c = CoxPHFitter().fit(d, "os_months", "os_event")
            print("PAM50", sub, "HR/SD=%.2f p=%.3f" % (c.hazard_ratios_["sig34"], c.summary.loc["sig34", "p"]))

# ------------------------------------------------------------------------------- C
if RUN_C:
    import cellxgene_census, scanpy as sc
    GENES = SIG34 + ["LDHA","HK2","PFKFB3","SLC2A1","SLC16A3","NDUFA1","NDUFB5","NDUFS1","COX5A","UQCRC1","ATP5F1A","SDHA",
                     "PDGFRA","DCN","LUM","CFD","C3","CXCL12","IL6","MYH11","RGS5","FOS","JUN","EGR1","ATF3","HSPA1A"]
    with cellxgene_census.open_soma(census_version="2025-11-08") as cen:
        ad = cellxgene_census.get_anndata(cen, "Homo sapiens",
            obs_value_filter="tissue_general == 'breast' and cell_type in ['fibroblast','fibroblast of mammary gland'] and is_primary_data == True",
            var_value_filter=f"feature_name in {sorted(set(GENES))}",
            obs_column_names=["donor_id", "dataset_id", "disease", "assay", "suspension_type"])
    ad.var_names = ad.var.feature_name
    ad.obs["group"] = np.where(ad.obs.disease.str.contains("breast", case=False) & ~ad.obs.disease.eq("normal"), "CAF",
                               np.where(ad.obs.disease.eq("normal"), "NF", "other"))
    ad = ad[ad.obs.group != "other"].copy()
    print(ad.obs.groupby(["group"]).agg(n_cells=("donor_id", "size"), n_donors=("donor_id", "nunique"), n_datasets=("dataset_id", "nunique")))
    # Pseudobulk: sum counts per donor, CPM-normalise, compare donors (Wilcoxon on donors) -> Squair et al. 2021
    X = pd.DataFrame(ad.X.toarray(), index=ad.obs_names, columns=ad.var_names)
    pb = X.groupby(ad.obs.donor_id).sum(); lib = X.sum(1).groupby(ad.obs.donor_id).sum()
    cpm = np.log2(pb.div(lib, axis=0) * 1e6 + 1)
    grp = ad.obs.groupby("donor_id").group.first().loc[cpm.index]
    rows = []
    for gname in cpm.columns:
        a, b = cpm.loc[grp == "CAF", gname], cpm.loc[grp == "NF", gname]
        rows.append([gname, a.mean() - b.mean(), stats.mannwhitneyu(a, b).pvalue, len(a), len(b)])
    from statsmodels.stats.multitest import multipletests
    pbde = pd.DataFrame(rows, columns=["gene", "log2FC_pseudobulk", "p", "n_CAF_donors", "n_NF_donors"])
    pbde["padj"] = multipletests(pbde.p, method="fdr_bh")[1]
    pbde.sort_values("p").to_csv(f"{OUT}/scRNA_pseudobulk_DE.csv", index=False)
    # Dissociation stress check: are CXCL1/IL6/MYC co-varying with immediate-early genes in NF?
    nf = X[ad.obs.group == "NF"]
    ieg = np.log1p(nf[["FOS", "JUN", "EGR1", "ATF3"]]).mean(1)
    for gname in ["CXCL1", "IL6", "MYC", "VEGFA", "NFKB1"]:
        print("NF: %s vs IEG score rho=%.2f" % (gname, stats.spearmanr(np.log1p(nf[gname]), ieg)[0]))
