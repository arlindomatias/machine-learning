# TCGA-BRCA: tumor vs. normal classification

Machine learning on the TCGA breast cancer cohort (TCGA-BRCA), configured in `config.yaml` and run with the
shared scripts in [`../scripts/classification/`](../scripts/classification):

- **Feature selection.** Genes are ranked by the between/within-group sum-of-squares ratio (BSS/WSS).
- **Classification.** Tumor vs. normal tissue from RNA-seq expression.
  - Every model listed under `models:` in `config.yaml` is trained and evaluated separately.
  - Each model's outputs are named after its id.
- **Clustering.** k-means and PAM on the top genes.
- **Enrichment.** Reactome over-representation of the top genes.
- **Copy number.** Somatic copy-number alteration burden.

**Report:** [arlindomatias.github.io/machine-learning/brca_tumor_normal_classification/reports/report.html](https://arlindomatias.github.io/machine-learning/brca_tumor_normal_classification/reports/report.html). It covers every step, has one section per model, and shows the plot code.

## Data

Broad GDAC Firehose level-3 files (stddata), stored outside the repository in `$DATA/TCGA-BRCA/`:

| File | Content |
|---|---|
| `BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt` | RSEM-normalized gene expression, 20,531 genes × 1,212 samples (1,100 tumor, 112 normal) |
| `BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg` | Somatic copy-number segments (log2 ratio), Affymetrix SNP 6.0, germline CNVs removed |

The class comes from the TCGA barcode: sample-type codes `01`–`09` are tumors and `10`–`19` are normals.

## Usage

```bash
conda activate genomic-variants        # Snakemake 9
snakemake -n --use-conda               # dry-run
snakemake --use-conda --cores 2
```

R runs in `envs/ml.yaml` (env `ml`), pinned by `envs/ml.linux-64.pin.txt`. The report uses `envs/report.yaml`.

### Adding or changing a model

Add an entry under `models:` in `config.yaml`:

```yaml
models:
  knn_k3:
    method: knn            # knn | random_forest | glmnet | svm | pamr | xgboost | dlda
    k: 3
    n_features: 100        # top training-set BSS/WSS genes, or "all"
```

Snakemake then runs only the new model and rebuilds the comparison and the report. Existing models are left untouched.
The options of each method are documented at the top of its file in
[`../scripts/classification/models/`](../scripts/classification/models).

### Adding a method

Each method is one file, `../scripts/classification/models/<method>.R`, that defines
`fit_model(X_train, y_train, X_test, mcfg, classes, groups)`. It returns the fitted object, the test-set predictions and
positive-class scores, and optionally per-gene weights. The contract is documented in `04_classification.R`, which does
everything shared by all models (data, gene selection, metrics, figures). If the method has per-gene weights, add it to
`has_gene_weights()` in the `Snakefile`. Add its package to `envs/ml.yaml` and regenerate the pin file.

## Pipeline

| Rule | Script | Main outputs |
|---|---|---|
| `prepare` | `01_prepare.R` | `data/expression_log2.rds`, `data/samples.csv` |
| `split` | `02_split.R` | `data/train_test_split.csv`, `results/classification/bss_wss_ranking_train.tsv` |
| `feature_selection` | `03_feature_selection.R` | `results/feature_selection/top_genes_full_cohort.tsv`, `plots/feature_selection/heatmap_top_genes.png` |
| `classify`, `classify_with_genes` (one job per model) | `04_classification.R` + `models/<method>.R` | `results/classification/<model>_{predictions,metrics,confusion_matrix}.tsv`, `plots/classification/<model>_confusion_matrix.png`, `data/models/<model>.rds`; methods with per-gene weights also write `<model>_genes.{tsv,png}` |
| `compare_models` | `05_compare_models.R` | `results/classification/{model_comparison,bootstrap_ci,pairwise_tests,misclassified_samples}.tsv`, `plots/classification/{model_comparison,roc_curves,calibration,misclassified_samples}.png` |
| `clustering` (one job per method) | `06_clustering.R` | `results/clustering/<method>_{assignments,contingency}.tsv`, `plots/clustering/<method>_clusplot.png` |
| `enrichment` | `07_enrichment.R` | `results/enrichment/reactome.tsv`, `plots/enrichment/reactome_{dotplot,emapplot}.png` |
| `copy_number` | `08_copy_number.R` | `results/copy_number/per_sample.tsv`, `plots/copy_number/fraction_genome_altered.png` |
| `session_info`, `report` | `session_info.R`, `report.py` | `reports/sessionInfo.txt`, `reports/report.html` |

`00_common.R` is sourced by every R step and handles logging, config, sample metadata and BSS/WSS.
Logs are written to `logs/`.

## Methods

- **Expression scale.** `log2(RSEM + 1)`.
- **Train/test split.** 70/30 by patient, stratified by whether the patient has a normal sample.
  A patient's tumor and normal never fall on different sides.
- **Training-set fitting only.** Gene selection (BSS/WSS) and kNN standardization are fit on the training set.
- **Evaluation.** On the test set (tumor = positive), settings under `evaluation:` in `config.yaml`:
  - Classification: sensitivity, specificity, balanced accuracy, PPV, NPV and MCC. These are used because the classes
    are imbalanced (~10 tumors per normal).
  - Ranking: ROC AUC, plus PR AUC with the normal (rare) class as the class of interest.
  - Probabilities: Brier score, log loss and calibration curves.
  - Uncertainty: patient-level bootstrap intervals (the same resamples for every model) and Wilson intervals for
    sensitivity and specificity.
  - Model vs. model: exact McNemar and DeLong tests for every pair, Holm-adjusted.
  - Error analysis: the test samples each model misclassifies.
- **Descriptive steps.** The full-cohort ranking is used only by the heatmap, clustering and enrichment.
