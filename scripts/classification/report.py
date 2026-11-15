"""HTML report for a tumor vs. normal classification project:
input -> feature selection -> classification (one section per model) -> clustering -> enrichment -> copy number.

Run through Snakemake (`script:` directive, rule `report`). Models, clustering methods and cutoffs come
from the project's config.yaml; every number is computed from the tables in results/.
"""
import os
import sys

import pandas as pd

# Shared report module: <repository>/scripts/html_report.py (this script lives in scripts/classification/).
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, prm, cfg = snakemake.input, snakemake.output, snakemake.params, snakemake.config  # noqa: F821
NEG, POS = cfg["classes"]["negative"], cfg["classes"]["positive"]
P = "plots"


def read_tsv(path, **kw):
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def by_stem(paths, suffix):
    """{id: path} for files named <id><suffix>."""
    return {os.path.basename(p)[: -len(suffix)]: p for p in paths if p.endswith(suffix)}


def describe_model(mcfg):
    n = mcfg["n_features"]
    feats = "all genes" if n == "all" else f"top {n} training-set BSS/WSS genes"
    method = mcfg["method"]
    if method == "knn":
        return f"k-nearest neighbours, k = {mcfg['k']}, {feats} (z-scored with training statistics)"
    if method == "random_forest":
        bal = "class-balanced " if mcfg.get("balanced") else ""
        return f"{bal}random forest, {mcfg['ntree']} trees, {feats}"
    if method == "glmnet":
        a = mcfg["alpha"]
        kind = "LASSO" if a == 1 else ("ridge" if a == 0 else f"elastic net (alpha = {a})")
        bal = ", class-weighted" if mcfg.get("balanced") else ""
        return (f"{kind} logistic regression{bal}, penalty tuned by patient-grouped CV "
                f"({mcfg.get('s', 'lambda.1se')}), {feats}")
    if method == "svm":
        kernel = mcfg.get("kernel", "linear")
        bal = ", class-weighted" if mcfg.get("balanced", True) else ""
        return f"SVM, {kernel} kernel, cost = {mcfg.get('cost', 1)}{bal}, {feats}"
    if method == "pamr":
        return (f"nearest shrunken centroids, {mcfg.get('prior', 'balanced')} prior, threshold tuned by "
                f"patient-grouped CV, {feats}")
    if method == "xgboost":
        return (f"gradient-boosted trees, {mcfg.get('nrounds', 200)} rounds, depth {mcfg.get('max_depth', 3)}, "
                f"eta {mcfg.get('eta', 0.1)}, {feats}")
    if method == "dlda":
        return f"diagonal linear discriminant analysis, {mcfg.get('prior', 'balanced')} prior, {feats}"
    return f"{method}, {feats}"


samples = pd.read_csv(inp.samples)
split = pd.read_csv(inp.split)
top = read_tsv(inp.top_genes)
comparison = read_tsv(inp.comparison)
reactome = read_tsv(inp.reactome)
cn = read_tsv(inp.copy_number)
confusion = by_stem(inp.confusion, "_confusion_matrix.tsv")
gene_tables = by_stem(inp.genes, "_genes.tsv")
contingency = by_stem(inp.contingency, "_contingency.tsv")
models = cfg["models"]

n_class = samples["class"].value_counts()
best = comparison.sort_values(["balanced_accuracy", "auc"], ascending=False).iloc[0]
n_top = cfg["feature_selection"]["n_top"]

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
r.kpis([
    (f"{len(samples):,}", "RNA-seq samples"),
    (f"{n_class.get(POS, 0):,} / {n_class.get(NEG, 0):,}", f"{POS} / {NEG}"),
    (len(models), "models compared"),
    (f"{best['balanced_accuracy']:.3f}", f"best balanced accuracy ({best['model']})"),
    (f"{len(reactome):,}", f"Reactome pathways (padj < {cfg['enrichment']['padj_cutoff']})"),
])

# ---- Data -----------------------------------------------------------------------------
r.h2("Data and design")
types = samples.groupby(["type_code", "class"]).size().reset_index(name="samples")
r.md(f"""
TCGA {cfg['dataset']} level-3 data from the Broad GDAC Firehose: RSEM-normalized RNA-seq gene expression
(analysed as `log2(RSEM + 1)`, {len(samples):,} samples) and somatic copy-number segments from Affymetrix
SNP 6.0 arrays. The class comes from the sample-type code of the TCGA barcode: codes below 10 are
**{POS}** (01 primary, 06 metastatic), codes 10 and above are **{NEG}** (11 solid tissue normal).

The classes are imbalanced (about {n_class.get(POS, 0) / max(n_class.get(NEG, 1), 1):.0f} tumors per normal),
so models are compared on sensitivity, specificity, balanced accuracy and AUC, not accuracy alone.
""")
r.table(types, "RNA-seq samples by TCGA sample-type code.")
r.md("""
**Pipeline**

1. **prepare**: parses the expression table into a samples × genes matrix.
2. **split**: patient-level, stratified train/test split; genes ranked by BSS/WSS on the training set only.
3. **feature_selection**: descriptive BSS/WSS ranking on the full cohort and heatmap.
4. **classify / classify_with_genes**: one job per model in `config.yaml`, with the method's code in `models/<method>.R`; outputs are named after the model id. Methods with per-gene weights also write a gene table.
5. **compare_models**: metrics table, bar chart and ROC curves of all models.
6. **clustering**: one job per method (k-means, PAM) on the top full-cohort genes.
7. **enrichment**: Reactome over-representation of the top genes.
8. **copy_number**: per-sample copy-number alteration burden.
""")

# ---- Feature selection ------------------------------------------------------------------
r.h2("Feature selection")
up = int((top["log2_fold_change"] > 0).sum()) if len(top) else 0
r.md(f"""
Genes are ranked by the ratio of between-group to within-group sum of squares (BSS/WSS; Dudoit et al.,
2002). Of the top {n_top} genes on the full cohort, **{up}** are higher and **{n_top - up}** lower in
{POS} samples. This ranking is descriptive (heatmap, clustering, enrichment); the classifiers use a
ranking computed on the training set only.
""")
r.figure(f"{P}/feature_selection/heatmap_top_genes.png",
         f"Top {n_top} genes, all {NEG} samples plus an equal number of random {POS} samples "
         "(row z-scores clipped to ±3).")
cols = ["rank", "symbol", "bss_wss", "log2_fold_change"]
r.table(top[cols], f"Top {n_top} genes on the full cohort (positive log2FC = higher in {POS}).",
        max_rows=25, collapsed=True)

# ---- Classification -----------------------------------------------------------------------
r.h2("Classification")
sets = split.groupby(["set", "class"]).size().unstack(fill_value=0).reset_index()
r.md(f"""
Patients were split {1 - cfg['split']['test_fraction']:.0%}/{cfg['split']['test_fraction']:.0%} into
training and test sets, stratified by whether they contributed a normal sample, so a patient's tumor and
normal never fall on different sides. Gene selection and standardization were fit on the training set.
""")
r.table(sets, "Samples per set.")
r.table(pd.DataFrame([{"model": m, "description": describe_model(c)} for m, c in models.items()]),
        "Models (`models:` in config.yaml).")
show = comparison[["model", "n_features", "accuracy", "sensitivity", "specificity", "balanced_accuracy", "auc"]]
r.table(show, f"Test-set performance ({POS} is the positive class).", digits=4)
r.figure(f"{P}/classification/model_comparison.png", "Test-set metrics per model.")
r.figure(f"{P}/classification/roc_curves.png",
         "ROC curves on the test set (kNN scores are neighbour vote shares; k = 1 gives a single point).")

for m, mcfg in models.items():
    r.h3(m)
    row = comparison[comparison["model"] == m].iloc[0]
    cm = read_tsv(confusion[m])
    fp = int(cm.query("actual == @NEG and predicted == @POS")["N"].sum())
    fn = int(cm.query("actual == @POS and predicted == @NEG")["N"].sum())
    r.md(f"""
{describe_model(mcfg)}. Balanced accuracy **{row['balanced_accuracy']:.3f}**, AUC **{row['auc']:.3f}**;
{fn} {POS} sample(s) called {NEG} and {fp} {NEG} sample(s) called {POS}.
""")
    r.figure(f"{P}/classification/{m}_confusion_matrix.png", f"Confusion matrix of `{m}` on the test set.",
             width="420px", code=[f"method_{mcfg['method']}", "confusion_matrix"])
    if m in gene_tables:
        gt = read_tsv(gene_tables[m])
        r.md(f"**{len(gt)}** gene(s) with a non-zero weight.")
        r.figure(f"{P}/classification/{m}_genes.png", f"Top genes of `{m}` by absolute weight.",
                 width="520px", code=["genes"])
        cols = [c for c in ["symbol", "weight", "direction"] if c in gt.columns]
        r.table(gt[cols], f"Per-gene weights of `{m}`.", max_rows=20, collapsed=True, digits=4)

# ---- Clustering -----------------------------------------------------------------------------
r.h2("Clustering")
r.md(f"""
Unsupervised clustering with k = {cfg['clustering']['k']} on the z-scored top {n_top} full-cohort genes.
Agreement is the share of samples whose cluster matches their class under the best one-to-one mapping
of clusters to classes.
""")
for method in cfg["clustering"]["methods"]:
    ct = read_tsv(contingency[method])
    wide = ct.pivot(index="cluster", columns="class", values="N").fillna(0).astype(int)
    agree = None
    if wide.shape == (2, 2):
        d = wide.values
        agree = max(d[0, 0] + d[1, 1], d[0, 1] + d[1, 0]) / d.sum()
    r.h3(method)
    if agree is not None:
        r.md(f"Agreement with the {POS}/{NEG} labels: **{agree:.1%}**.")
    r.table(wide, f"`{method}` clusters vs. class.", index=True)
    r.figure(f"{P}/clustering/{method}_clusplot.png",
             f"`{method}` clusters on the first two principal components, points colored by true class.",
             width="620px", code=["clusplot"])

# ---- Enrichment --------------------------------------------------------------------------------
r.h2("Reactome enrichment")
r.md(f"""
Over-representation of the top {n_top} genes in Reactome pathways (ReactomePA), with every gene on the
platform as the universe and Benjamini–Hochberg adjustment. **{len(reactome)}** pathways have
adjusted p < {cfg['enrichment']['padj_cutoff']}.
""")
r.figure(f"{P}/enrichment/reactome_dotplot.png", "Top enriched Reactome pathways.")
r.figure(f"{P}/enrichment/reactome_emapplot.png", "Enrichment map: pathways linked by shared genes.")
if len(reactome):
    r.table(reactome[["Description", "GeneRatio", "BgRatio", "p.adjust", "geneID"]],
            "Enriched Reactome pathways.", max_rows=30, collapsed=True)

# ---- Copy number ---------------------------------------------------------------------------------
r.h2("Copy number")
thr = cfg["copy_number"]["gain_loss_threshold"]
fga = cn.groupby("class")["frac_genome_altered"].agg(["count", "median"]).reset_index()
fga.columns = ["class", "samples", "median fraction of genome altered"]
r.md(f"""
Somatic copy-number segments ({len(cn):,} samples, germline CNVs removed). The fraction of genome altered
is the share of the segmented genome with |log2 ratio| > {thr}. In the copy-number data, normals are mostly
matched blood samples (code 10).
""")
r.table(fga, "Copy-number alteration burden by class.")
r.figure(f"{P}/copy_number/fraction_genome_altered.png", "Fraction of genome altered by tissue type.")

# ---- Appendix --------------------------------------------------------------------------------------
outputs = [
    ("data/samples.csv", "Sample table (barcode, patient, sample-type code, class)"),
    ("data/train_test_split.csv", "Train/test assignment of every sample"),
    ("results/feature_selection/top_genes_full_cohort.tsv", "Top genes on the full cohort"),
    ("results/classification/bss_wss_ranking_train.tsv", "Gene ranking on the training set"),
    ("results/classification/<model>_{predictions,metrics,confusion_matrix}.tsv", "Per-model test-set outputs"),
    ("results/classification/<model>_genes.tsv", "Per-gene weights (methods that have them)"),
    ("results/classification/model_comparison.tsv", "Metrics of every model"),
    ("results/clustering/<method>_{assignments,contingency}.tsv", "Per-method clustering outputs"),
    ("results/enrichment/reactome.tsv", "Reactome enrichment"),
    ("results/copy_number/per_sample.tsv", "Per-sample copy-number summary"),
    ("data/models/<model>.rds", "Fitted model objects (not versioned)"),
]
r.appendix(outputs, env_files=list(inp.envs), config_files=list(inp.config), session_file=inp.session)
r.write(out.html)
