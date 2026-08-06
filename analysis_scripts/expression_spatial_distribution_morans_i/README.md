# Expression, spatial distribution, and Moran’s I (MATLAB)

MATLAB analysis scripts for comparing **WT** vs **5xFAD** spatial transcriptomics data (ExSeq and Xenium): differential expression by region, Moran’s I spatial autocorrelation, and ExSeq spatial distribution / localization along a registered hippocampus reference map.

## Requirements

- **MATLAB** (R2020b or later recommended)
- **Bioinformatics Toolbox** — required for `rnaseqde` in the expression scripts (default statistical test)
- **Parquet support** — required only for `CompareExpressionBetweenGroups_Xenium_FullSections_AllGenes.m` (`parquetread` / related functions)

Transcript data files are **not** bundled. Place your files under `data/` as described below, then edit each script’s EXAMPLE sample list to match your filenames.

## Folder layout

```text
expression_spatial_distribution_morans_i/
  README.md
  data/
    exseq/                    # ExSeq transcript CSVs
    xenium/                   # Xenium transcript CSVs
    xenium_full_sections/     # Xenium full-section Parquet files
  results/                    # created/used by scripts (default outputs)
  CompareExpressionBetweenGroups_ExSeq.m
  CompareExpressionBetweenGroups_Xenium.m
  CompareExpressionBetweenGroups_Xenium_FullSections_AllGenes.m
  CompareMoranBetweenGroups_ExSeq.m
  CompareMoranBetweenGroups_Xenium.m
  ExSeqSpatialDistributionXY.m
  exseq_spatial_calls.csv     # gene–region calls for spatial screening
  hippocampus_reference_map.png
```

Paths are resolved relative to each script’s location via `fileparts(mfilename('fullpath'))`. Edit the **USER SETTINGS** block at the top of a script if your data live elsewhere.

## Quick start

1. Copy your transcript files into the matching `data/` subfolder.
2. Open the script you want to run.
3. In **USER SETTINGS**, confirm `inputFolder` / `outputFolder` (defaults usually need no change).
4. Replace the **EXAMPLE sample list** filenames with your file basenames (keep the `name` / `group` / `files` structure).
5. Run the script from MATLAB (current folder can be anywhere; paths are script-relative).

### Sample list rules (all scripts)

- Each sample needs `.name`, `.group`, and `.files` (cell array of basenames in `inputFolder`).
- `.group` must be exactly `"WT"` or `"FAD"` (5xFAD animals use `"FAD"`).
- Keep an **equal number** of WT and FAD samples (scripts error otherwise).
- At least **two** samples per group are required for statistical testing.
- Multiple files under one sample are merged (expression / Moran) or reassembled (spatial registration).

## Input schemas

### Expression scripts (ExSeq / Xenium CSV)

Required columns:

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `region_name` | Anatomical region label |
| `cell_type` | Optional; required only if `useCellTypes = true` |

One row per transcript (or molecule). Counts are aggregated inside the script.

### Expression script (Xenium full sections)

Same columns as above, stored in **Parquet** files (read by row groups).

### Moran scripts (ExSeq / Xenium)

Required columns (aliases accepted where noted in the script):

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `region_name` | Region label (`region` / `hippocampal_region` also accepted) |
| `global_x_um` | X coordinate in µm (`x_um` / `x` also accepted) |
| `global_y_um` | Y coordinate in µm (`y_um` / `y` also accepted) |
| `Z_um` | Z coordinate in µm (`global_z_um` / `z_um` / `z` also accepted) |

### ExSeq spatial distribution

- Transcript CSVs: same spatial columns as the Moran scripts.
- Calls table (`exseq_spatial_calls.csv` by default): columns `analysis`, `gene`, `region`, where `analysis` is `expression` or `moran`.
- Reference image: `hippocampus_reference_map.png` by default.

## Scripts

### `CompareExpressionBetweenGroups_ExSeq.m`

Differential expression by hippocampal region for ExSeq.

| Setting | Default |
|---------|---------|
| Input | `data/exseq/` |
| Output | `results/exseq_expression/` |
| Example design | 6 WT + 6 FAD animals |

Main outputs include raw and analysis count tables, filter diagnostics, and `all_WT_vs_5xFAD_region_results.xlsx` / significant-result tables.

### `CompareExpressionBetweenGroups_Xenium.m`

Same pipeline for Xenium hippocampus CSVs.

| Setting | Default |
|---------|---------|
| Input | `data/xenium/` |
| Output | `results/xenium_expression/` |
| Example design | 2 WT + 2 FAD animals |

### `CompareExpressionBetweenGroups_Xenium_FullSections_AllGenes.m`

Expression comparison across **all mapped brain regions** from full-section Parquet files (large inputs; gene-name diagnostics are disabled for speed).

| Setting | Default |
|---------|---------|
| Input | `data/xenium_full_sections/` |
| Output | `results/xenium_full_sections_expression/` |
| Example design | 2 WT + 2 FAD animals |

### `CompareMoranBetweenGroups_ExSeq.m`

Moran’s I per gene × region × animal (10 µm voxels), then WT vs FAD testing.

| Setting | Default |
|---------|---------|
| Input | `data/exseq/` |
| Output | `results/exseq_moran/` |
| Example design | 6 WT + 6 FAD animals |

### `CompareMoranBetweenGroups_Xenium.m`

Moran’s I per **section**, then sections are averaged to animals in `averageSectionsToAnimals` before group testing.

| Setting | Default |
|---------|---------|
| Input | `data/xenium/` |
| Output | `results/xenium_moran/` |
| Example design | 2 WT + 2 FAD animals after averaging |

If you change section/animal naming, update both the EXAMPLE sample list and `averageSectionsToAnimals` so section columns map to the animal-level names in `sampleNames`.

### `ExSeqSpatialDistributionXY.m`

Registers each ExSeq animal to `hippocampus_reference_map.png`, then tests genotype-associated shifts along reference X/Y for gene–region pairs listed in `exseq_spatial_calls.csv`.

| Setting | Default |
|---------|---------|
| Input | `data/exseq/` |
| Output | `results/exseq_spatial_distribution/` |
| Calls | `exseq_spatial_calls.csv` |
| Reference | `hippocampus_reference_map.png` |

Call from MATLAB, for example:

```matlab
ExSeqSpatialDistributionXY();   % uses defaults above
```

Optional arguments: `outputFolder`, `inputFolder`, `referenceImagePath`, `callsCsvPath`, bootstrap replicate count, and a flag to generate significant-pair plots.

**Multi-section orientation:** `sectionOrientation` and `sectionVerticalAdjustment` default to no correction. If your sections need mirroring, rotation, or vertical offsets before whole-animal reassembly, edit those helpers for your sample names and section indices.

## Notes

- Group labels in code are `"WT"` and `"FAD"` even when referring to 5xFAD biologically.
- Analysis options (normalization, filters, FDR cutoffs, statistical tests) remain editable in each script’s settings blocks; defaults match the published analysis choices where applicable.
- No proprietary absolute paths are required; everything defaults under this folder.
