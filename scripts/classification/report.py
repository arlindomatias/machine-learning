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
GROUP = cfg.get("group_name", "group")
FEAT = cfg.get("feature_name", "feature")
MODULES = cfg.get("modules", {})


def enabled(module):
    return MODULES.get(module, {}).get("enabled", False)


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
    feats = f"all {FEAT}s" if n == "all" else f"top {n} training-set BSS/WSS {FEAT}s"
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
        return (f"{kind} logistic regression{bal}, penalty tuned by {GROUP}-grouped CV "
                f"({mcfg.get('s', 'lambda.1se')}), {feats}")
    if method == "svm":
        kernel = mcfg.get("kernel", "linear")
        bal = ", class-weighted" if mcfg.get("balanced", True) else ""
        return f"SVM, {kernel} kernel, cost = {mcfg.get('cost', 1)}{bal}, {feats}"
    if method == "pamr":
        bal = "class-balanced " if mcfg.get("balanced", True) else ""
        return f"{bal}nearest shrunken centroids, threshold tuned by {GROUP}-grouped CV, {feats}"
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
boot = read_tsv(inp.bootstrap)
pairwise = read_tsv(inp.pairwise)
errors = read_tsv(inp.errors)
features = pd.read_csv(inp.features)
reactome = read_tsv(inp.reactome) if enabled("enrichment") else pd.DataFrame()
cn = read_tsv(inp.copy_number) if enabled("copy_number") else pd.DataFrame()
confusion = by_stem(inp.confusion, "_confusion_matrix.tsv")
gene_tables = by_stem(inp.genes, "_genes.tsv")
contingency = by_stem(inp.contingency, "_contingency.tsv")
models = cfg["models"]

n_class = samples["class"].value_counts()
best = comparison.sort_values(["balanced_accuracy", "auc"], ascending=False).iloc[0]
n_top = cfg["feature_selection"]["n_top"]

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
kpis = [
    (f"{len(samples):,}", "samples"),
    (f"{n_class.get(POS, 0):,} / {n_class.get(NEG, 0):,}", f"{POS} / {NEG}"),
    (len(models), "models compared"),
    (f"{best['balanced_accuracy']:.3f}", f"best balanced accuracy ({best['model']})"),
]
if enabled("enrichment"):
    kpis.append((f"{len(reactome):,}", f"Reactome pathways (padj < {MODULES['enrichment']['padj_cutoff']})"))
r.kpis(kpis)

# ---- Data -----------------------------------------------------------------------------
r.h2("Data and design")
icfg = cfg["input"]
big, small = (POS, NEG) if n_class.get(POS, 0) >= n_class.get(NEG, 0) else (NEG, POS)
ratio = n_class.get(big, 0) / max(n_class.get(small, 1), 1)
r.md(cfg.get("description", ""))
r.md(f"""
**Input**: format `{icfg['format']}`, transform `{icfg.get('transform', 'none')}`; {len(samples):,} samples ×
{len(features):,} {FEAT}s after the prepare step, in {samples['group'].nunique():,} {GROUP}s.
The classes are {"imbalanced" if ratio >= 1.5 else "roughly balanced"} (about {ratio:.1f} {big} per {small}),
so models are compared on sensitivity, specificity, balanced accuracy, MCC and AUC, not accuracy alone.
""")
by = [c for c in cfg.get("report", {}).get("sample_columns", []) if c in samples.columns] + ["class"]
r.table(samples.groupby(by).size().reset_index(name="samples"), f"Samples by {', '.join(by)}.")
steps = [
    f"**prepare**: reads the raw input with `inputs/{icfg['format']}.R` and writes the expression matrix, sample table and {FEAT} table.",
    f"**split**: {GROUP}-level, stratified train/test split; {FEAT}s ranked by BSS/WSS on the training set only.",
    "**feature_selection**: descriptive BSS/WSS ranking on the full cohort and heatmap.",
    f"**classify / classify_with_genes**: one job per model in `config.yaml`, with the method's code in `models/<method>.R`; outputs are named after the model id. Methods with per-{FEAT} weights also write a {FEAT} table.",
    f"**compare_models**: metrics with {GROUP} bootstrap intervals, pairwise McNemar and DeLong tests, ROC and calibration curves, and misclassified samples.",
]
if enabled("clustering"):
    steps.append(f"**clustering**: one job per method ({', '.join(MODULES['clustering']['methods'])}) on the top full-cohort {FEAT}s.")
if enabled("enrichment"):
    steps.append(f"**enrichment**: Reactome over-representation of the top {FEAT}s.")
if enabled("copy_number"):
    steps.append("**copy_number**: per-sample copy-number alteration burden.")
r.md("**Pipeline**\n\n" + "\n".join(f"{i}. {t}" for i, t in enumerate(steps, 1)))

# ---- Feature selection ------------------------------------------------------------------
r.h2("Feature selection")
up = int((top["log2_fold_change"] > 0).sum()) if len(top) else 0
r.md(f"""
{FEAT.capitalize()}s are ranked by the ratio of between-group to within-group sum of squares (BSS/WSS; Dudoit et al.,
2002). Of the top {n_top} {FEAT}s on the full cohort, **{up}** are higher and **{n_top - up}** lower in
{POS} samples. This ranking is descriptive (heatmap, clustering, enrichment); the classifiers use a
ranking computed on the training set only.
""")
r.figure(f"{P}/feature_selection/heatmap_top_genes.png",
         f"Top {n_top} {FEAT}s, every {small} sample plus an equal number of random {big} samples "
         "(row z-scores clipped to ±3).")
cols = ["rank", "symbol", "bss_wss", "log2_fold_change"]
r.table(top[cols], f"Top {n_top} {FEAT}s on the full cohort (positive log2FC = higher in {POS}).",
        max_rows=25, collapsed=True)

# ---- Classification -----------------------------------------------------------------------
r.h2("Classification")
sets = split.groupby(["set", "class"]).size().unstack(fill_value=0).reset_index()
if "set" in samples.columns:
    split_text = "The train/test split was predefined in the sample table."
else:
    split_text = (f"{GROUP.capitalize()}s were split {1 - cfg['split']['test_fraction']:.0%}/"
                  f"{cfg['split']['test_fraction']:.0%} into training and test sets, stratified by whether they "
                  f"contain a {NEG} sample.")
r.md(f"""
{split_text} All samples of a {GROUP} fall on the same side. {FEAT.capitalize()} selection and standardization
were fit on the training set.
""")
r.table(sets, "Samples per set.")
r.table(pd.DataFrame([{"model": m, "description": describe_model(c)} for m, c in models.items()]),
        "Models (`models:` in config.yaml).")
ev = cfg["evaluation"]
level = f"{ev['ci_level']:.0%}"


def with_ci(model, metric, digits=3):
    row = boot[(boot["model"] == model) & (boot["metric"] == metric)].iloc[0]
    return f"{row['estimate']:.{digits}f} [{row['lower']:.{digits}f}, {row['upper']:.{digits}f}]"


r.h3("Performance")
r.md(f"""
All metrics are on the test set, with **{POS}** as the positive class.

- **Sensitivity / specificity**: share of {POS} / {NEG} samples called correctly. **Balanced accuracy** is their mean.
- **PPV / NPV**: share of {POS} / {NEG} calls that are correct. With few {NEG} samples, NPV is the more fragile of the two.
- **MCC** (Matthews correlation): uses all four cells of the confusion matrix. It is 1 for a perfect classifier and 0 at chance, and it is not inflated by the class imbalance.
- **ROC AUC**: ranking quality of the scores. **PR AUC ({NEG})** is the area under the precision-recall curve with {NEG}, the rare class, as the class of interest (chance level = its share of the test set, {(split.query("set == 'test'")["class"] == NEG).mean():.3f}).
- **Brier score / log loss**: quality of the probabilities themselves (lower is better). Log loss clips probabilities to [{ev['log_loss_eps']}, {1 - ev['log_loss_eps']}].

Intervals are {level} percentile intervals from {ev['bootstrap']:,} bootstrap resamples of test-set **{GROUP}s** (all samples of a
{GROUP} stay together). Every model is evaluated on the same resamples. Sensitivity and specificity also have
Wilson intervals in `model_comparison.tsv`.
""")
cols = ["model", "n_features", "sensitivity", "specificity", "balanced_accuracy", "ppv", "npv", "mcc",
        "auc", "pr_auc_negative", "brier", "log_loss"]
r.table(comparison[cols], f"Test-set performance ({POS} is the positive class).", digits=4)
ci_metrics = {"balanced_accuracy": "balanced accuracy", "mcc": "MCC", "auc": "ROC AUC",
              "npv": "NPV", "brier": "Brier score"}
r.table(pd.DataFrame([{"model": m, **{label: with_ci(m, k) for k, label in ci_metrics.items()}} for m in models]),
        f"Point estimate [{level} {GROUP} bootstrap interval].")
r.figure(f"{P}/classification/model_comparison.png",
         f"Test-set metrics per model with {level} {GROUP} bootstrap intervals.")
r.figure(f"{P}/classification/roc_curves.png",
         "ROC curves on the test set (kNN scores are neighbour vote shares; k = 1 gives a single point).")

r.h3("Pairwise comparisons")
n_pairs = len(pairwise)
sig_m = int((pairwise["mcnemar_p_holm"] < 0.05).sum())
sig_d = int((pairwise["delong_p_holm"] < 0.05).sum())
r.md(f"""
Every pair of models is compared on the same test samples. The **exact McNemar test** asks whether one model's errors are
systematically different from the other's, using only the samples where exactly one of the two is right. The
**DeLong test** compares the two AUCs, accounting for their correlation. P-values are Holm-adjusted across the
{n_pairs} pairs. **{sig_m}** pair(s) differ by McNemar and **{sig_d}** by DeLong after adjustment.
""")
top = pairwise.sort_values("mcnemar_p").head(15)[
    ["model_a", "model_b", "errors_a", "errors_b", "only_a_correct", "only_b_correct", "mcnemar_p",
     "mcnemar_p_holm", "auc_a", "auc_b", "delong_p", "delong_p_holm"]]
r.table(top, "The 15 pairs with the smallest McNemar p-values (all pairs in `pairwise_tests.tsv`).",
        digits=4, collapsed=True)

r.h3("Calibration")
r.md(f"""
For each bin of predicted P({POS}) ({ev['calibration_bins']} equal-width bins), the observed share of {POS} samples.
A well-calibrated model lies on the diagonal. Most predictions fall in the two extreme bins, so the middle bins
hold only a few samples (point size). kNN with k = 1 gives only 0 or 1.
""")
r.figure(f"{P}/classification/calibration.png", "Calibration of the predicted probabilities on the test set.")

r.h3("Misclassified samples")
if len(errors):
    always = errors[errors["n_models_wrong"] == len(models)]
    r.md(f"""
**{len(errors)}** test sample(s) are misclassified by at least one model and **{len(always)}** by every model.
Samples that every model gets wrong may have a questionable label or composition (for example a tumor sample
with a high share of normal tissue) rather than reveal a weakness of any one model.
""")
    r.table(errors[["sample_id", "actual", "n_models_wrong", "models_wrong", "mean_prob_positive"]],
            "Test samples misclassified by at least one model.", digits=3)
    r.figure(f"{P}/classification/misclassified_samples.png",
             "Predictions of every model for the misclassified samples (red = wrong).")
else:
    r.md("No test sample is misclassified by any model.")

for m, mcfg in models.items():
    r.h3(m)
    row = comparison[comparison["model"] == m].iloc[0]
    cm = read_tsv(confusion[m])
    fp = int(cm.query("actual == @NEG and predicted == @POS")["N"].sum())
    fn = int(cm.query("actual == @POS and predicted == @NEG")["N"].sum())
    r.md(f"""
{describe_model(mcfg)}. Balanced accuracy **{with_ci(m, 'balanced_accuracy')}**, MCC **{with_ci(m, 'mcc')}**,
AUC **{row['auc']:.3f}**, Brier score **{row['brier']:.4f}**;
{fn} {POS} sample(s) called {NEG} and {fp} {NEG} sample(s) called {POS}.
""")
    r.figure(f"{P}/classification/{m}_confusion_matrix.png", f"Confusion matrix of `{m}` on the test set.",
             width="420px", code=[f"method_{mcfg['method']}", "confusion_matrix"])
    if m in gene_tables:
        gt = read_tsv(gene_tables[m])
        r.md(f"**{len(gt)}** {FEAT}(s) with a non-zero weight.")
        r.figure(f"{P}/classification/{m}_genes.png", f"Top {FEAT}s of `{m}` by absolute weight.",
                 width="520px", code=["genes"])
        cols = [c for c in ["symbol", "weight", "direction"] if c in gt.columns]
        r.table(gt[cols], f"Per-{FEAT} weights of `{m}`.", max_rows=20, collapsed=True, digits=4)

# ---- Clustering -----------------------------------------------------------------------------
if enabled("clustering"):
    r.h2("Clustering")
    r.md(f"""
    Unsupervised clustering with k = {MODULES['clustering']['k']} on the z-scored top {n_top} full-cohort {FEAT}s.
    Agreement is the share of samples whose cluster matches their class under the best one-to-one mapping
    of clusters to classes.
    """)
    for method in MODULES["clustering"]["methods"]:
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
if enabled("enrichment"):
    r.h2("Reactome enrichment")
    r.md(f"""
    Over-representation of the top {n_top} {FEAT}s in Reactome pathways (ReactomePA), with every measured {FEAT} with an Entrez id as the
    universe and Benjamini–Hochberg adjustment. **{len(reactome)}** pathways have
    adjusted p < {MODULES['enrichment']['padj_cutoff']}.
    """)
    r.figure(f"{P}/enrichment/reactome_dotplot.png", "Top enriched Reactome pathways.")
    r.figure(f"{P}/enrichment/reactome_emapplot.png", "Enrichment map: pathways linked by shared genes.")
    if len(reactome):
        r.table(reactome[["Description", "GeneRatio", "BgRatio", "p.adjust", "geneID"]],
                "Enriched Reactome pathways.", max_rows=30, collapsed=True)

# ---- Copy number ---------------------------------------------------------------------------------
if enabled("copy_number"):
    r.h2("Copy number")
    thr = MODULES["copy_number"]["gain_loss_threshold"]
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
    ("data/samples.csv", "Sample table (sample_id, group, class and extra columns)"),
    ("data/features.csv", f"{FEAT.capitalize()} table (feature_id, symbol, entrez)"),
    ("data/train_test_split.csv", "Train/test assignment of every sample"),
    ("results/feature_selection/top_genes_full_cohort.tsv", f"Top {FEAT}s on the full cohort"),
    ("results/classification/bss_wss_ranking_train.tsv", f"{FEAT.capitalize()} ranking on the training set"),
    ("results/classification/<model>_{predictions,metrics,confusion_matrix}.tsv", "Per-model test-set outputs"),
    ("results/classification/<model>_genes.tsv", f"Per-{FEAT} weights (methods that have them)"),
    ("results/classification/model_comparison.tsv", "Metrics of every model (with Wilson intervals)"),
    ("results/classification/bootstrap_ci.tsv", f"{GROUP.capitalize()} bootstrap intervals of every metric"),
    ("results/classification/pairwise_tests.tsv", "McNemar and DeLong tests for every pair of models"),
    ("results/classification/misclassified_samples.tsv", "Test samples misclassified by at least one model"),
    *([("results/clustering/<method>_{assignments,contingency}.tsv", "Per-method clustering outputs")]
      if enabled("clustering") else []),
    *([("results/enrichment/reactome.tsv", "Reactome enrichment")] if enabled("enrichment") else []),
    *([("results/copy_number/per_sample.tsv", "Per-sample copy-number summary")] if enabled("copy_number") else []),
    ("data/models/<model>.rds", "Fitted model objects (not versioned)"),
]
r.appendix(outputs, env_files=list(inp.envs), config_files=list(inp.config), session_file=inp.session)
r.write(out.html)
