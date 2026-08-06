# Supervised STELLAR encoder label transfer for ExSeq cell-type annotation

This repository contains the code used to transfer broad cell-type labels from an annotated spatial reference dataset to an unlabeled ExSeq tissue. It uses the graph encoder published with **STELLAR** and trains that encoder with a standard supervised cross-entropy objective on the reference labels.

> **Important methodological distinction**  
> This is a **supervised adaptation of the STELLAR encoder**. It does **not** run the original open-world/semi-supervised STELLAR training objective, pseudo-labeling procedure, or novel-cell-type discovery step. Therefore, the output classes are restricted to the classes represented in the annotated reference.

The original STELLAR method and software are described in:

> Brbić M, Cao K, Hickey JW, et al. Annotation of spatially resolved single-cell data with STELLAR. *Nature Methods*. 2022;19:1411–1418. doi:10.1038/s41592-022-01651-8  
> Original software: https://github.com/snap-stanford/stellar

## Repository structure

```text
.
├── run_label_transfer.py                  # command-line entry point
├── stellar_transfer/                      # reusable pipeline modules
│   ├── config.py                          # YAML configuration
│   ├── data.py                            # input parsing and gene alignment
│   ├── graph.py                           # KNN/radius graph construction
│   ├── model.py                           # STELLAR Encoder import and training
│   ├── pipeline.py                        # end-to-end workflow and outputs
│   └── plotting.py                        # diagnostic spatial plots
├── configs/
│   └── manuscript_config_template.yaml    # analysis settings from the manuscript workflow
├── examples/
│   ├── toy_reference.csv                  # small synthetic annotated reference
│   ├── toy_target_counts.csv              # small synthetic target count table
│   ├── toy_target_points.csv              # small synthetic target point metadata
│   ├── toy_config.yaml                    # runnable toy configuration
│   └── generate_toy_data.py               # regenerates the synthetic data
├── tests/test_core.py
├── requirements.txt
├── environment.yml
└── CITATION.cff
```

## Reference dataset

### Biological role

The reference is an **annotated spatial cell-by-gene table**. Each reference cell supplies three types of information:

1. a gene-expression vector used as node features;
2. an XY position used to construct the spatial graph; and
3. an existing cell-type annotation used as the supervised training label.

The reference should be generated independently of the target tissue and should contain representative cells from every class that may be assigned to the target. Because this implementation is closed-set, a target cell cannot be assigned to a biologically new class that is absent from the reference.

### Required reference columns

The reference CSV contains one row per cell.

| Column | Required | Description |
|---|---:|---|
| `cell_id` | yes | Unique cell identifier. A non-numeric first column can also be interpreted as the ID column. |
| `global_x`, `global_y` | yes | Spatial cell coordinates used to build the graph. Units must be consistent within the reference. |
| label column | yes | Raw annotation column, configured as `reference.label_column`; in the manuscript workflow this is `auto_label_v5`. |
| confidence column | conditional | Used when `reference.keep_confidence` is not `null`; in the manuscript workflow this is `auto_confidence_v5`. |
| gene columns | yes | Numeric raw count columns, one column per gene. Gene names are matched to the target case-insensitively. |
| `total_counts` | optional | Used for reference filtering when present; otherwise it is calculated from the selected gene columns. |
| `region_name`, `Z`, `fov` | optional | Retained as metadata but not used as model features. |

For maximum reproducibility, list the gene columns explicitly under `reference.gene_columns`. When this field is `null`, the loader treats numeric columns as genes after excluding known metadata columns and any columns listed under `reference.non_gene_columns`.

### Reference preprocessing in the manuscript configuration

The supplied manuscript template performs the following operations:

- removes background IDs ending in `_0`;
- retains reference annotations with high or medium confidence;
- excludes unresolved/unknown labels;
- harmonizes detailed labels into the broad classes `Excitatory_Neuron`, `GABA_Neuron`, `Astro-Epen`, `Immune`, `Vascular`, and `OPC-Oligo`;
- retains cells with at least 40 total transcripts;
- removes classes represented by fewer than 10 retained cells;
- caps each reference class at 250 cells using a fixed random seed;
- normalizes every cell to a total of 10,000 and applies `log(1+x)` transformation.

The exact raw-to-broad label mapping is stored in the YAML configuration rather than hard-coded in Python.

### Included small reference

`examples/toy_reference.csv` is a deterministic **synthetic** reference provided only to demonstrate file formatting and execution. It contains artificial marker-enriched counts for six broad classes. It is not derived from the study, must not be used for biological interpretation, and does not replace the full study reference.

To regenerate all toy files:

```bash
python examples/generate_toy_data.py
```

## Target inputs

Two target files are required.

### 1. Target count table

`target_count_csv` is a genes-by-cells CSV:

```text
gene,C1_1,C1_2,C1_3,...
Slc17a7,12,3,0,...
Gad1,0,18,1,...
...
```

The first column must be named `gene`. Remaining columns are cell IDs. Duplicate gene names and duplicate normalized cell IDs are summed.

### 2. Target point metadata

`target_points_csv` contains one or more spatial points per cell and must include:

| Column | Required | Description |
|---|---:|---|
| `cell` | yes | Cell assignment for the point/transcript. |
| `global_x`, `global_y` | yes | Point coordinates. Cell centroids are calculated as the mean coordinates. |
| `region_name` | optional | Used only for regional summaries and the plotting-only CA/DG excitatory color split. |
| `Z` | optional | Used to select a central Z slab for RNA-position convex-hull visualization. |
| `fov` | optional | Retained in the prediction table. |

Cell IDs such as `C4_123` and `C04-123` are normalized to `C04_123` before counts and metadata are matched.

## Graph construction

A separate spatial graph is built for the reference and target. The default manuscript configuration uses a 3-nearest-neighbor graph (`graph.mode: knn`, `graph.k: 3`) and adds reciprocal edges (`make_undirected: true`). Radius graphs are also supported.

The upstream STELLAR `GraphDataset` receives expression matrices with shape `[cells, shared genes]`, integer reference labels, and edge lists for the reference and target. Only genes shared between both datasets are used.

## Model training

The code imports `GraphDataset` and `Encoder` directly from a local checkout of the original STELLAR repository. The encoder consists of the published STELLAR input projection, GraphSAGE layer, and normalized linear classifier. It is optimized on the annotated reference graph using `CrossEntropyLoss`.

Manuscript-template training settings:

```yaml
training:
  seed: 1
  target_sum: 10000
  epochs: 60
  learning_rate: 0.001
  weight_decay: 0.0
```

`graph.partition_parts: 4` reproduces the partitioned-reference training setting used in the analysis script. Set it to `1` for full-graph training when the graph is small or the installed PyTorch Geometric build does not provide a graph-partitioning backend.

The reported `reference_training_accuracy` is a **training-set sanity check**, not an independent validation score.

## Installation

### 1. Clone this repository and the original STELLAR code

```bash
git clone <THIS-REPOSITORY-URL>
cd supervised-stellar-exseq-label-transfer
mkdir -p external
git clone https://github.com/snap-stanford/stellar.git external/stellar
```

The original STELLAR repository is not copied into this archive; it remains a separately cited dependency under its own license.

### 2. Create an environment

The original STELLAR authors reported compatibility with Python 3.8, PyTorch 1.9.1, and PyTorch Geometric 2.0. A matching environment template is supplied:

```bash
conda env create -f environment.yml
conda activate supervised-stellar-transfer
```

Alternatively, install PyTorch and PyTorch Geometric using builds appropriate for the local CUDA/CPU setup, then install the remaining packages:

```bash
pip install -r requirements.txt
```

## Run the synthetic example

After cloning upstream STELLAR into `external/stellar`:

```bash
python run_label_transfer.py --config examples/toy_config.yaml
```

The example uses full-graph training and writes results to `outputs/toy_run/`.

## Run the manuscript workflow

Copy and edit the template:

```bash
cp configs/manuscript_config_template.yaml configs/my_run.yaml
```

Update the five paths at the top of the new YAML file, then run:

```bash
python run_label_transfer.py --config configs/my_run.yaml
```

No study-specific absolute paths are embedded in the Python source.

## Main outputs

Each run creates an output directory containing:

| Output | Description |
|---|---|
| `target_predictions_supervised_stellar_encoder.csv` | Cell IDs, coordinates, predicted labels, confidence, and per-class probabilities. |
| `target_prediction_counts.csv` | Predicted broad-class counts. |
| `prediction_counts_by_region_diagnostic_only.csv` | Region-by-predicted-class diagnostic table. |
| `reference_training_check.csv` | Reference true/predicted labels and confidence; training sanity check only. |
| `shared_genes.txt` | Ordered genes used as model features. |
| `label_mapping.json` | Broad label-to-integer mapping. |
| `training_history.csv` | Epoch loss and reference training accuracy. |
| `model_state_dict.pt` | Trained encoder weights. |
| `config_used.yaml` | Exact configuration copied into the result directory. |
| `run_summary.json` | Data dimensions, graph statistics, parameters, package versions, and prediction counts. |
| `spatial_predictions.*` | Spatial centroid prediction plot. |
| `spatial_prediction_confidence.png` | Spatial maximum-class-probability plot. |
| `spatial_predictions_rna_hulls.*` | Convex hulls of RNA positions; a visualization aid, not true cell morphology. |

The optional `CA_GLUT`/`DG_GLUT` split is stored in `plot_label` and is used only for plotting. The actual model prediction remains `Excitatory_Neuron` in `pred_label`.

## Suggested Methods wording

> Broad cell-type labels were transferred from an annotated spatial reference to target ExSeq tissues using a supervised adaptation of the STELLAR graph encoder. Reference and target cells were represented by log-normalized expression values of genes shared between the two datasets, and separate undirected spatial 3-nearest-neighbor graphs were constructed from cell coordinates. The STELLAR encoder was trained on the annotated reference graph for 60 epochs using cross-entropy loss, a learning rate of 0.001, and a fixed random seed of 1. The trained encoder was then applied to the target graph to obtain broad cell-type predictions and per-class probabilities. This supervised implementation used the published STELLAR encoder architecture but did not use the original semi-supervised/open-world objective or novel-class discovery procedure.

Adjust this paragraph when the final configuration differs from the supplied manuscript template.

## Reproducibility checks

Run the non-STELLAR unit tests with:

```bash
pytest
```

The tests verify ID normalization, normalization behavior, graph construction, toy reference parsing, target parsing, and shared-gene alignment. A complete model run additionally requires the upstream STELLAR checkout and a compatible PyTorch Geometric installation.

## Citation and licensing

Please cite both the study repository and the original STELLAR publication. `CITATION.cff` contains the STELLAR DOI and a placeholder software citation entry that should be updated with the final repository URL, complete author list, and release identifier.

Before public release, add a license selected by the study authors/institution. The original STELLAR dependency has its own license and should not be relicensed by this repository.
