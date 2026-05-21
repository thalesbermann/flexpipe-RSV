# flexpipe-RSV

Nextstrain pipeline for genomic epidemiology of Respiratory Syncytial Virus (RSV-A and RSV-B). This tool is derived from [flexpipe](https://github.com/InstitutoTodosPelaSaude/flexpipe) and adapted for RSV, providing a complete workflow for fetching, curating, subsampling, and phylogenetically analysing RSV genomes.

This repository contains all essential files to generate RSV Nextstrain builds for both subtypes A and B. Using this pipeline, users can perform genomic epidemiology analyses, visualize phylogeographic results, and track RSV spread based on genomic data and associated metadata.

## Getting Started

To run this pipeline, see the instructions available in the original [flexpipe repository](https://github.com/InstitutoTodosPelaSaude/flexpipe), which covers Unix CLI navigation, installation of a Nextstrain environment with conda/mamba, and a step-by-step tutorial on generating a Nextstrain build (preparing, aligning, and visualizing genomic data).

---

## Builds

| Build | Subtype | Data source | Nextclade / ViralQC dataset |
|-------|---------|-------------|----------------------------|
| RSV_A | RSV-A | Pathoplexus (primary) + NCBI | ViralQC (`rsv-a`) |
| RSV_B | RSV-B | Pathoplexus (primary) + NCBI | ViralQC (`rsv-b`) |

Each build lives in its own subdirectory (`RSV_A/` and `RSV_B/`) with independent `config/`, `data/`, `ingest/`, `phylogenetic/`, and `scripts/` folders.

---

## Pipeline Overview

```
fetch_pathoplexus  (or fetch_ncbi)
    └── merge_local_sequences  (ITpS sequences + Pathoplexus/NCBI)
            └── viralqc        (BLAST + Nextclade QC + clade assignment)
                    └── curate_qc  (normalisation, dedup, filters)
                            └── prepare  (subsampling)
                                    ├── coordinates  (geocoding → latlongs.tsv)
                                    ├── generate_name2hue  (colour palette)
                                    └── colours  (colour_scheme.tsv)
                                            └── [phylogenetic/Snakefile]
                                                    align → mask → tree → refine
                                                    → ancestral → translate → traits
                                                    → clades → export → auspice/results.json
```

---

## Stage 1 — Ingest

Sequences and metadata are fetched from **Pathoplexus** by default (with NCBI as fallback), using organism-specific endpoints. The `data_source` parameter in `config.yaml` selects the active source — only one is used per run.

| Parameter | RSV_A | RSV_B |
|-----------|-------|-------|
| `data_source` | `pathoplexus` | `pathoplexus` |
| `pathoplexus.organism` | `rsv-a` | `rsv-b` |
| `pathoplexus.min_completeness` | 0.70 | 0.70 |
| `ncbi.taxid` (fallback) | 12234 | 12233 |
| `ncbi.genome_size` (bp) | 15225 | 15222 |
| `ncbi.min_length` | 70% | 70% |
| `ncbi.max_length` | 110% | 110% |

Local ITpS sequences (in `data/new_sequences.fasta` + `data/metadata.xlsx`) are merged with Pathoplexus sequences at this stage.

---

## Stage 2 — QC and Curation

**ViralQC** (BLAST + Nextclade) assigns genome quality grades (A–D) and clade labels to every sequence. The `curate.py` script then:

- Renames and standardises metadata fields (`strain`, `date`, `country`, `division`, `location`, `data_use`, `clade`)
- Truncates clade names to a configurable number of hierarchy levels (`clade_levels`) for display grouping
- Infers geographic regions from country names
- Deduplicates sequences, preferring local ITpS records

**QC filters** applied by `augur filter`:

| Parameter | Value |
|-----------|-------|
| `qc.genome_quality` | `A`, `B` (grades C and D discarded) |
| `qc.min_coverage` | 0.70 |
| Required columns | `strain`, `date`, `country`, `clade` |

### Clade truncation

RSV follows a deep hierarchical nomenclature (e.g., `A.D.3.3.1`). The pipeline truncates to 3 levels for colour grouping:

| Build | `clade_levels` | Example |
|-------|---------------|---------|
| RSV_A | 3 | `A.D.3.3.1` → `A.D.3` |
| RSV_B | 3 | `B.D.E.1.1` → `B.D.E` |

---

## Stage 3 — Subsampling

Controlled by `config/subsample.yaml`. Strategy: **focal** sequences (ITpS + all Brazil) are always kept in full; **context** sequences (outside Brazil) are subsampled by country × year × clade_truncated.

```yaml
defaults:
  min_date: 2000

samples:
  focal:
    query: "(source == 'ITpS') or (country == 'Brazil')"

  context:
    group_by: [country, year, clade_truncated]
    sequences_per_group: 4
    exclude_where:
      - "country=Brazil"
      - "clade_truncated="
      - "date="
```

Adjust `sequences_per_group` and `min_date` to control dataset size and temporal depth. RSV uses `min_date: 2000` to capture the full modern diversity; flu pipelines use `min_date: 2015`.

---

## Stage 4 — Coordinates and Colours

**Coordinates**: `get_coordinates.py` queries Nominatim (OpenStreetMap) to geocode `country`, `division`, and `location` fields. Results are cached in `config/cache_coordinates.tsv` to avoid redundant API calls. The output `config/latlongs.tsv` is consumed by `augur export`.

**Colours**: `generate_name2hue.py` assigns hues to each unique value in `clade_truncated`, `source`, `data_use`, and geographic columns. `colour_maker.py` produces the final `config/colour_scheme.tsv`.

Colour columns configured in `config.yaml`:

```yaml
colours:
  clade:    "clade_truncated clade"
  geo:      "region country division location"
  source:   "source"
  data_use: "data_use"
```

---

## Stage 5 — Phylogenetic

Run separately after ingest completes:

```bash
snakemake --snakefile phylogenetic/Snakefile --cores 5
```

Steps: `align` (MAFFT) → `mask` → `tree` (IQ-TREE 3 UFBoot) → `refine` (TreeTime) → `ancestral` → `translate` → `traits` → `clades` → `export` → `auspice/results.json`

> **Note:** RSV genomes are ~15 kb (full genome). IQ-TREE can take several hours on large datasets. Masking removes 44 bp at the 5′ end and 155 bp at the 3′ end to eliminate low-quality terminal regions common in RSV assemblies.

### Key phylogenetic parameters (`config.yaml`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `parameters.model` | `MFP` | ModelFinder Plus — auto-selects best substitution model |
| `parameters.ufboot` | `1000` | Ultrafast bootstrap replicates |
| `parameters.root` | `least-squares` | Root method for time-calibrated tree |
| `parameters.coalescent` | `skyline` | Effective population size model in TreeTime |
| `parameters.date_inference` | `marginal` | Marginal date inference for ambiguous dates |
| `parameters.divergence_units` | `mutations` | Branch length units in timetree |
| `parameters.clock_filter_iqd` | `4` | IQD filter for clock outliers |
| `parameters.ancestral_inference` | `joint` | Joint ancestral reconstruction |
| `parameters.mask_5prime` | `44` | Bases masked at 5′ end (RSV-specific) |
| `parameters.mask_3prime` | `155` | Bases masked at 3′ end (RSV-A specific, Use 145 for RSV-B) |
| `options.threads` | `5` | Threads for MAFFT and IQ-TREE |
| `traits.columns` | `country division location clade` | Columns for ancestral trait inference |

---

## Configuration Files

| File | Purpose |
|------|---------|
| `config/config.yaml` | All pipeline parameters (Pathoplexus/NCBI, QC, phylogenetic, colours) |
| `config/subsample.yaml` | Subsampling strategy and group sizes |
| `config/auspice_config.json` | Auspice display settings (colorings, filters, panels) |
| `config/reference.gb` | Full-genome reference sequence in GenBank format |
| `config/clades.tsv` | Clade definitions for `augur clades` |
| `config/keep.txt` | Strains to always include |
| `config/ignore.txt` | Strains to always exclude |
| `data/new_sequences.fasta` | Local ITpS sequences |
| `data/metadata.xlsx` | Local ITpS metadata |

---

## Author

**Thales Bermann** — Instituto Todos pela Saúde (ITpS)
✉️ [thalesbermann@gmail.com](mailto:thalesbermann@gmail.com)

---

## License

This project is licensed under the [MIT License](LICENSE).
