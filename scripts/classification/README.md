# Binary classification workflow (shared)

Snakemake workflow and R scripts for binary classification from an expression-like matrix (genes, probes,
proteins): feature ranking, a configurable set of models, model comparison, and optional clustering,
enrichment and copy-number modules. A project folder holds only its `config.yaml` and a short `Snakefile`.
Everything specific to a dataset is in the config.

Used by: [`../../brca_tumor_normal_classification/`](../../brca_tumor_normal_classification).

## Starting a new project

1. Put the raw files under `$DATA/<dir>/` (never in the repository).
2. Create `<project>/Snakefile`:

   ```python
   configfile: "config.yaml"

   SCRIPTS = "../scripts/classification"


   include: f"{SCRIPTS}/workflow.smk"


   rule all:
       input:
           rules.all_outputs.input,
   ```

   `rule all` must live in the project's Snakefile: Snakemake takes the default target from the main Snakefile,
   not from included files.
3. Copy a `config.yaml` (e.g. the BRCA one) and edit `dataset`, `title`, `description`, `input`, `classes` and
   `modules`. Models, evaluation and the split usually stay as they are.
4. Run it:

   ```bash
   conda activate genomic-variants        # Snakemake 9
   cd <project>
   snakemake -n --use-conda               # dry-run
   snakemake --use-conda --cores 2
   ```

## Input contract

The `prepare` rule (`01_prepare.R` + an adapter in `inputs/`) is the only dataset-specific step. It writes three
files, and every later step reads only these three:

**`data/expression.rds`**: numeric matrix.
- Rows are samples, named by `sample_id`, in the order of `samples.csv`.
- Columns are features, named by `feature_id`.
- Values are on the analysis scale (e.g. log2). No `NA`, `NaN` or `Inf`.

**`data/samples.csv`**: one row per matrix row.

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique id, equal to the matrix row names. |
| `group` | yes | Unit that must not be split between train and test (patient, donor). Equal to `sample_id` when samples are independent. |
| `class` | yes | `classes.negative` or `classes.positive`. Samples of other classes are dropped by `prepare`, which logs how many. |
| `set` | no | `train` / `test`. If present, it fixes the split (predefined or external validation set). |
| others | no | Kept unchanged (e.g. `batch`, `tissue`). `report.sample_columns` in the config breaks the sample counts down by them. |

**`data/features.csv`**: one row per matrix column.

| Column | Required | Meaning |
|---|---|---|
| `feature_id` | yes | Equal to the matrix column names. |
| `symbol` | yes | Label used in tables and plots (`feature_id` when there is no annotation). |
| `entrez` | yes, may be empty | Entrez Gene id. Only the enrichment module uses it. |

`validate_inputs()` (in `00_common.R`) checks all of the above at the end of `prepare` and stops with a list of
every violation. It also requires exactly two classes. Multiclass problems are not supported: the metrics and
several methods assume a positive and a negative class.

### Sample order

The random train/test split depends on the order of the samples (the order of their groups, to be exact). The order is
the one the adapter returns: the file order for `firehose_rsem`, and the sample sheet order for `matrix`. The same input
files always give the same split. Reordering the sample sheet gives a different, equally valid split.

## Input formats (`input.format`)

Each format is one adapter file, `inputs/<format>.R`. It defines `load_input(icfg, files, classes)` and returns
`list(X, samples, features)` on the raw scale. `01_prepare.R` then maps the class labels (`classes.map`), drops
samples of other classes, applies `input.transform`, validates and writes.

| Format | Input | Samples and features |
|---|---|---|
| `firehose_rsem` | Broad GDAC Firehose level-3 RSEM gene expression file (`expression:`). | From the TCGA barcode: `group` = participant, `class` from the sample-type code (< 10 positive), `type_code` kept. Features from the `SYMBOL\|ENTREZ` ids. |
| `matrix` | Any delimited matrix (`expression:`; first column = row ids) with `orientation: features_x_samples` (default, GEO style) or `samples_x_features`. | A sample sheet (`samples:`) and optionally a feature annotation (`features:`), with a column mapping. |

`input.transform`: `log2p1` (log2(x + 1), for RSEM/TPM-like values) or `none` (data already on a log scale,
e.g. microarrays).

Example `matrix` input for a GEO microarray:

```yaml
input:
  dir: GSE12345                 # $DATA/GSE12345/
  format: matrix
  transform: none
  expression: GSE12345_matrix.tsv
  orientation: features_x_samples
  samples:
    file: samples.tsv
    columns: {sample_id: geo_accession, group: patient, class: tissue}   # group is optional
  features:                     # optional
    file: GPL570.tsv
    columns: {feature_id: ID, symbol: Gene Symbol, entrez: ENTREZ_GENE_ID}

classes:
  negative: Normal
  positive: Tumor
  map: {"normal breast tissue": Normal, "breast tumor": Tumor}
```

Matrix samples that are not in the sheet are dropped and logged. Sheet rows that are not in the matrix are an error.
Entrez ids that map to several genes (`1 /// 2`) are set to empty.

To support a new format, add `inputs/<format>.R` with `load_input()`. Nothing else changes.

## Config reference

| Key | Meaning |
|---|---|
| `dataset`, `title` | Report title and subtitle. |
| `description` | Markdown for the report's data section. |
| `group_name`, `feature_name` | Words used in the report for a group (e.g. `patient`) and a feature (e.g. `gene`). |
| `input` | See above. |
| `classes` | `negative`, `positive`, optional `map` (raw label -> class) and `colors`. |
| `feature_selection.n_top` | Size of the descriptive full-cohort ranking (heatmap, clustering, enrichment). |
| `split.test_fraction` | Share of groups in the test set. |
| `models` | One entry per model: `method`, `n_features` (number of top training-set BSS/WSS features, or `all`) and method options. The options are documented at the top of `models/<method>.R`. |
| `evaluation` | `ci_level`, `bootstrap`, `log_loss_eps`, `calibration_bins`. |
| `modules` | `clustering`, `enrichment` (`organism`: a ReactomePA organism), `copy_number` (TCGA Firehose segments only). Each has `enabled`. Disabled modules are left out of the workflow and the report. |
| `report.sample_columns` | Extra `samples.csv` columns for the sample table in the report. |
| `seed` | Random seed for every step. |

## Models

Each method is one file, `models/<method>.R`, that defines
`fit_model(X_train, y_train, X_test, mcfg, classes, groups)`. The contract is documented in `04_classification.R`,
which does everything shared by all models: data, feature selection, metrics and figures.

Available methods: `knn`, `random_forest`, `glmnet`, `svm`, `pamr`, `xgboost` and `dlda`.

To add a method:
1. Write `models/<method>.R`.
2. If it returns per-feature weights, add it to `has_gene_weights()` in `workflow.smk`.
3. Add its package to `envs/ml.yaml` and regenerate `envs/ml.linux-64.pin.txt`. Install the new package with
   `--freeze-installed` on top of the current pin, so existing versions do not change.

## Files

| File | Role |
|---|---|
| `workflow.smk` | Rules. Included by each project's Snakefile. |
| `00_common.R` | Sourced by every R step: logging, config, input contract (`validate_inputs()`, readers), BSS/WSS, metrics. |
| `01_prepare.R`, `inputs/` | Raw input to the contract files. `inputs/tcga.R` holds TCGA helpers shared with `08_copy_number.R`. |
| `02_split.R` | Group-level stratified split and training-set feature ranking. |
| `03_feature_selection.R` | Full-cohort ranking and heatmap. |
| `04_classification.R`, `models/` | One model per job. |
| `05_compare_models.R` | Metrics, bootstrap intervals, pairwise tests, ROC, calibration, misclassified samples. |
| `06_clustering.R`, `07_enrichment.R`, `08_copy_number.R` | Optional modules. |
| `report.py` | HTML report (uses `../html_report.py`). |
| `envs/` | Conda environments `ml` (R) and `report` (Python), with pin files. |
