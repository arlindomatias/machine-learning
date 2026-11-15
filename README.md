# Machine learning

Snakemake + R workflows applying machine learning to public cancer genomics data.

| Project | Dataset | Role |
|---|---|---|
| [`brca_tumor_normal_classification/`](brca_tumor_normal_classification) | TCGA-BRCA | Tumor vs. normal classification from RNA-seq (kNN, random forest), clustering, Reactome enrichment, copy number ([report](https://arlindomatias.github.io/machine-learning/brca_tumor_normal_classification/reports/report.html)) |

## Layout

Projects follow the same layout as the `bulk-rna-seq` repository:

| Folder | Contents |
|---|---|
| `<project>/data/` | Inputs and intermediate objects (`*.rds` and fitted models are not versioned) |
| `<project>/results/` | Final results, **tables only** (`.tsv`/`.csv`) |
| `<project>/plots/` | All figures |
| `<project>/envs/` | Conda environments (`envs/<name>.yaml` with `name: <name>`, plus `<name>.linux-64.pin.txt`) |
| `<project>/reports/` | Self-contained `report.html` (built in Python) and `sessionInfo.txt` |
| `<project>/logs/`, `cache/` | Rule logs and downloaded annotation (not versioned) |
| `scripts/` | Shared scripts: `html_report.py` (copied from `bulk-rna-seq/scripts/`, keep identical) and one folder per pipeline (`scripts/classification/`) |

- **Heavy raw files never go into the repository:** they live under `$DATA/<dataset>/`, with
  `$DATA=/mnt/sda2/HD/Datasets`. Snakefiles read the path from the environment.
- Dataset-specific values (inputs, classes, models, cutoffs, seed) go in each project's `config.yaml`.
- Every rule has a `log:`.
- Reports are versioned and published with GitHub Pages (branch `main`, root, `.nojekyll`):
  `https://arlindomatias.github.io/machine-learning/<project>/reports/report.html`.
- Plot code is wrapped in `# <figure name="...">` / `# </figure>` markers, so the report can show it.

## Running

```bash
conda activate genomic-variants        # Snakemake 9
cd <project>
snakemake -n --use-conda               # dry-run
snakemake --use-conda --cores 2
```
