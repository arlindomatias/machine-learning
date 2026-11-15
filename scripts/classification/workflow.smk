# ============================================================
# Shared Snakemake workflow: binary classification from an expression matrix, plus optional
# clustering, enrichment and copy-number modules. A project folder holds only config.yaml and a
# Snakefile:
#
#   configfile: "config.yaml"
#   SCRIPTS = "../scripts/classification"      # this folder, relative to the project
#   include: f"{SCRIPTS}/workflow.smk"
#   rule all:
#       input: rules.all_outputs.input
#
# Every model in `models:` and every enabled module gets its own outputs, named after its id.
# Raw inputs are read from $DATA/<input.dir>/ (never stored in the repository).
# `script:` and `conda:` paths below are relative to this file; input/output paths are relative
# to the project folder.
# ============================================================
import os

ENV = "envs/ml.yaml"
ENV_REPORT = "envs/report.yaml"
RAW = os.path.join(os.environ["DATA"], config["input"]["dir"])

MODELS = list(config["models"])
MODULES = config.get("modules", {})


def enabled(module):
    return MODULES.get(module, {}).get("enabled", False)


def has_gene_weights(mcfg):
    """Methods whose fitted model gives one weight per feature (written to <model>_genes.tsv)."""
    if mcfg["method"] == "svm":
        return mcfg.get("kernel", "linear") == "linear"
    return mcfg["method"] in {"random_forest", "glmnet", "pamr", "xgboost", "dlda"}


GENE_MODELS = [m for m in MODELS if has_gene_weights(config["models"][m])]
PLAIN_MODELS = [m for m in MODELS if m not in GENE_MODELS]
CLUSTERING = MODULES["clustering"]["methods"] if enabled("clustering") else []


def one_of(ids):
    return "|".join(ids) if ids else "^$"


def raw_input_files():
    """Raw files of the input adapter, by role (expression, samples, features)."""
    icfg = config["input"]
    files = {"expression": os.path.join(RAW, icfg["expression"])}
    for role in ("samples", "features"):
        if isinstance(icfg.get(role), dict) and "file" in icfg[role]:
            files[role] = os.path.join(RAW, icfg[role]["file"])
    return files


MODULE_OUTPUTS = {
    "clustering": expand("results/clustering/{method}_assignments.tsv", method=CLUSTERING),
    "enrichment": ["results/enrichment/reactome.tsv"],
    "copy_number": ["results/copy_number/per_sample.tsv"],
}


# Every final output. Snakemake takes the default target from the project's Snakefile, which
# therefore defines `rule all` with `input: rules.all_outputs.input`.
rule all_outputs:
    input:
        "results/classification/model_comparison.tsv",
        *[f for m, files in MODULE_OUTPUTS.items() if enabled(m) for f in files],
        "reports/sessionInfo.txt",
        "reports/report.html",


# ------------------------------------------------------------
# 01: raw input -> the three files of the input contract (adapter inputs/<format>.R)
# ------------------------------------------------------------
rule prepare:
    input:
        **raw_input_files(),
        adapter=f"{SCRIPTS}/inputs/{config['input']['format']}.R",
    output:
        expression="data/expression.rds",
        samples="data/samples.csv",
        features="data/features.csv",
    conda:
        ENV
    log:
        "logs/prepare.log",
    script:
        "01_prepare.R"


# ------------------------------------------------------------
# 02: group-level train/test split + feature ranking on the training set
# ------------------------------------------------------------
rule split:
    input:
        expression=rules.prepare.output.expression,
        samples=rules.prepare.output.samples,
        features=rules.prepare.output.features,
    output:
        split="data/train_test_split.csv",
        ranking="results/classification/bss_wss_ranking_train.tsv",
    conda:
        ENV
    log:
        "logs/split.log",
    script:
        "02_split.R"


# ------------------------------------------------------------
# 03: descriptive feature ranking on the full cohort + heatmap
# ------------------------------------------------------------
rule feature_selection:
    input:
        expression=rules.prepare.output.expression,
        samples=rules.prepare.output.samples,
        features=rules.prepare.output.features,
    output:
        top_genes="results/feature_selection/top_genes_full_cohort.tsv",
        heatmap="plots/feature_selection/heatmap_top_genes.png",
    conda:
        ENV
    log:
        "logs/feature_selection.log",
    script:
        "03_feature_selection.R"


# ------------------------------------------------------------
# 04: one job per model. Every method runs the same script, which sources the method's file in
# models/; methods with per-feature weights also write <model>_genes.{tsv,png}.
# ------------------------------------------------------------
def method_file(wildcards):
    return f"{SCRIPTS}/models/{config['models'][wildcards.model]['method']}.R"


CLASSIFY_INPUT = dict(
    expression=rules.prepare.output.expression,
    features=rules.prepare.output.features,
    split=rules.split.output.split,
    ranking=rules.split.output.ranking,
    method=method_file,
)
CLASSIFY_OUTPUT = dict(
    model="data/models/{model}.rds",
    predictions="results/classification/{model}_predictions.tsv",
    metrics="results/classification/{model}_metrics.tsv",
    confusion="results/classification/{model}_confusion_matrix.tsv",
    confusion_plot="plots/classification/{model}_confusion_matrix.png",
)


rule classify:
    input:
        **CLASSIFY_INPUT,
    output:
        **CLASSIFY_OUTPUT,
    wildcard_constraints:
        model=one_of(PLAIN_MODELS),
    conda:
        ENV
    log:
        "logs/classification/{model}.log",
    script:
        "04_classification.R"


rule classify_with_genes:
    input:
        **CLASSIFY_INPUT,
    output:
        **CLASSIFY_OUTPUT,
        genes="results/classification/{model}_genes.tsv",
        genes_plot="plots/classification/{model}_genes.png",
    wildcard_constraints:
        model=one_of(GENE_MODELS),
    conda:
        ENV
    log:
        "logs/classification/{model}.log",
    script:
        "04_classification.R"


# ------------------------------------------------------------
# 05: comparison of every model: metrics, group bootstrap intervals, pairwise tests (McNemar,
# DeLong), calibration and misclassified samples
# ------------------------------------------------------------
rule compare_models:
    input:
        metrics=expand("results/classification/{model}_metrics.tsv", model=MODELS),
        predictions=expand("results/classification/{model}_predictions.tsv", model=MODELS),
    output:
        comparison="results/classification/model_comparison.tsv",
        bootstrap="results/classification/bootstrap_ci.tsv",
        pairwise="results/classification/pairwise_tests.tsv",
        errors="results/classification/misclassified_samples.tsv",
        comparison_plot="plots/classification/model_comparison.png",
        roc_plot="plots/classification/roc_curves.png",
        calibration_plot="plots/classification/calibration.png",
        errors_plot="plots/classification/misclassified_samples.png",
    conda:
        ENV
    log:
        "logs/compare_models.log",
    script:
        "05_compare_models.R"


# ------------------------------------------------------------
# 06 (module clustering): one job per clustering method
# ------------------------------------------------------------
rule clustering:
    input:
        expression=rules.prepare.output.expression,
        samples=rules.prepare.output.samples,
        top_genes=rules.feature_selection.output.top_genes,
    output:
        assignments="results/clustering/{method}_assignments.tsv",
        contingency="results/clustering/{method}_contingency.tsv",
        clusplot="plots/clustering/{method}_clusplot.png",
    wildcard_constraints:
        method=one_of(CLUSTERING),
    conda:
        ENV
    log:
        "logs/clustering/{method}.log",
    script:
        "06_clustering.R"


# ------------------------------------------------------------
# 07 (module enrichment): Reactome enrichment of the top features
# ------------------------------------------------------------
rule enrichment:
    input:
        features=rules.prepare.output.features,
        top_genes=rules.feature_selection.output.top_genes,
    output:
        table="results/enrichment/reactome.tsv",
        dotplot="plots/enrichment/reactome_dotplot.png",
        emapplot="plots/enrichment/reactome_emapplot.png",
    conda:
        ENV
    log:
        "logs/enrichment.log",
    script:
        "07_enrichment.R"


# ------------------------------------------------------------
# 08 (module copy_number, TCGA Firehose only): copy-number segments
# ------------------------------------------------------------
rule copy_number:
    input:
        copy_number=os.path.join(RAW, MODULES.get("copy_number", {}).get("file", "")),
    output:
        per_sample="results/copy_number/per_sample.tsv",
        fga_plot="plots/copy_number/fraction_genome_altered.png",
    conda:
        ENV
    log:
        "logs/copy_number.log",
    script:
        "08_copy_number.R"


rule session_info:
    output:
        session="reports/sessionInfo.txt",
    conda:
        ENV
    log:
        "logs/session_info.log",
    script:
        "session_info.R"


# ------------------------------------------------------------
# Report: HTML summary of the whole analysis (plots/ + results/). Sections of disabled modules
# are left out.
# ------------------------------------------------------------
def report_input(_):
    d = dict(
        samples=rules.prepare.output.samples,
        features=rules.prepare.output.features,
        split=rules.split.output.split,
        top_genes=rules.feature_selection.output.top_genes,
        comparison=rules.compare_models.output.comparison,
        bootstrap=rules.compare_models.output.bootstrap,
        pairwise=rules.compare_models.output.pairwise,
        errors=rules.compare_models.output.errors,
        confusion=expand("results/classification/{model}_confusion_matrix.tsv", model=MODELS),
        genes=expand("results/classification/{model}_genes.tsv", model=GENE_MODELS),
        contingency=expand("results/clustering/{method}_contingency.tsv", method=CLUSTERING),
        plots=[
            rules.feature_selection.output.heatmap,
            *expand("plots/classification/{model}_confusion_matrix.png", model=MODELS),
            *expand("plots/classification/{model}_genes.png", model=GENE_MODELS),
            rules.compare_models.output.comparison_plot,
            rules.compare_models.output.roc_plot,
            rules.compare_models.output.calibration_plot,
            rules.compare_models.output.errors_plot,
            *expand("plots/clustering/{method}_clusplot.png", method=CLUSTERING),
        ],
        session=rules.session_info.output.session,
        config=["config.yaml"],
        envs=[f"{SCRIPTS}/envs/ml.yaml", f"{SCRIPTS}/envs/report.yaml"],
        scripts=[f"{SCRIPTS}/report.py", f"{SCRIPTS}/../html_report.py"],
        code=[
            f"{SCRIPTS}/03_feature_selection.R",
            f"{SCRIPTS}/04_classification.R",
            *sorted({f"{SCRIPTS}/models/{m['method']}.R" for m in config["models"].values()}),
            f"{SCRIPTS}/05_compare_models.R",
            *([f"{SCRIPTS}/06_clustering.R"] if enabled("clustering") else []),
            *([f"{SCRIPTS}/07_enrichment.R"] if enabled("enrichment") else []),
            *([f"{SCRIPTS}/08_copy_number.R"] if enabled("copy_number") else []),
        ],
    )
    if enabled("enrichment"):
        d["reactome"] = rules.enrichment.output.table
        d["plots"] += [rules.enrichment.output.dotplot, rules.enrichment.output.emapplot]
    if enabled("copy_number"):
        d["copy_number"] = rules.copy_number.output.per_sample
        d["plots"] += [rules.copy_number.output.fga_plot]
    return d


rule report:
    input:
        unpack(report_input),
    output:
        html="reports/report.html",
    params:
        title=f"{config['dataset']}: {config['classes']['positive']} vs. {config['classes']['negative']} classification",
        subtitle=config["title"],
    conda:
        ENV_REPORT
    log:
        "logs/report.log",
    script:
        "report.py"
