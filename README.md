# Raw2Ready Skills

A collection of **QoderWork Agent Skills** for bioinformatics workflows — from raw GEO data to analysis-ready AnnData datasets.

Each skill is a self-contained module that teaches an AI agent how to perform a specific task, complete with utility scripts and step-by-step instructions.

## Available Skills

| Skill | Description | Dependencies |
|-------|-------------|--------------|
| [**geo-download**](skills/geo-download/) | Download GEO datasets, SOFT metadata files, and expression matrices by accession ID (GSE/GSM). Supports batch download and URL generation. | `requests`, `GEOparse` |
| [**parse-geo-soft**](skills/parse-geo-soft/) | Parse GEO SOFT format files into structured CSV metadata tables. Auto-discovers sample characteristics columns. | `pandas` |
| [**geo-matrix-merge**](skills/geo-matrix-merge/) | Merge GEO single-cell matrix files (10X, counts+cellname, GSM matrices) into a single AnnData `.h5ad` file. Auto-detects format. | `pandas`, `numpy`, `anndata`, `scipy` |
| [**geo-ready**](skills/geo-ready/) | End-to-end pipeline: download a GEO dataset → parse metadata → filter non-human samples → merge matrices → produce analysis-ready `.h5ad`. Orchestrates the three skills above. | All of the above |

## Pipeline Overview

```
                    ┌──────────────────┐
                    │   GEO Accession  │
                    │   (e.g. GSE...)  │
                    └────────┬─────────┘
                             │
              ┌──────────────┼──────────────┐
              │              │              │
              ▼              ▼              ▼
     ┌────────────┐  ┌─────────────┐  ┌────────────┐
     │geo-download│  │parse-geo-   │  │geo-matrix- │
     │  (fetch    │  │soft (extract│  │merge (build│
     │  RAW.tar + │  │ metadata    │  │ AnnData    │
     │  SOFT.gz)  │  │ CSV)        │  │ h5ad)      │
     └─────┬──────┘  └──────┬──────┘  └──────┬─────┘
           │                │                │
           └────────────────┼────────────────┘
                            │
                            ▼
                    ┌───────────────┐
                    │   geo-ready   │
                    │ (all-in-one   │
                    │  pipeline)    │
                    └───────────────┘
```

**geo-ready** calls the other three skills internally. You can use it as a single command, or run each skill individually for more control.

## Quick Start

### One-command pipeline (geo-ready)

```bash
python3 skills/geo-ready/scripts/geo_ready.py GSE213835 -d /path/to/downloads
```

This downloads the dataset, parses metadata, filters non-human samples, merges matrices, and writes `GSE213835_merged.h5ad`.

### Step-by-step

```bash
# 1. Download
python3 skills/geo-download/scripts/geo_download.py soft GSE213835 -o ./GSE213835
curl -L -o ./GSE213835/GSE213835_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213835/suppl/GSE213835_RAW.tar"

# 2. Parse metadata
gunzip -k ./GSE213835/GSE213835_family.soft.gz
python3 skills/parse-geo-soft/scripts/parse_soft.py \
  ./GSE213835/GSE213835_family.soft -o ./GSE213835/metadata.csv

# 3. Extract and merge
mkdir -p ./GSE213835/matrix_files
tar xf ./GSE213835/GSE213835_RAW.tar -C ./GSE213835/matrix_files
python3 skills/geo-matrix-merge/scripts/geo_matrix_merge.py \
  ./GSE213835/matrix_files -o ./GSE213835/output/GSE213835_merged.h5ad
```

## Installation

### Option 1 — Copy skill directories

Download the skill directory you need and place it in your QoderWork skills folder:

| OS | Path |
|----|------|
| macOS / Linux | `~/.qoderworkcn/skills/<skill-name>/` |
| Windows | `%USERPROFILE%\.qoderworkcn\skills\<skill-name>\` |

For example:
```bash
cp -r skills/geo-ready ~/.qoderworkcn/skills/
```

### Option 2 — Clone the whole repo

```bash
git clone git@github.com:Haley202305/Raw2Ready-Skills.git
# Then symlink or copy individual skills to ~/.qoderworkcn/skills/
```

## Skill Structure

Every skill follows this layout:

```
skill-name/
├── SKILL.md          # Agent instructions (required)
└── scripts/          # Utility scripts
    └── tool.py
```

The agent reads `SKILL.md` to learn **what** to do and **when**, then calls scripts in `scripts/` to execute.

## Requirements

Common Python dependencies across all skills:

```bash
pip install pandas numpy anndata scipy requests GEOparse
```

Some skills also require system tools: `curl`, `tar`, `gunzip`.

## Contributing

To add a new skill:

1. Create `skills/your-skill-name/SKILL.md` with YAML frontmatter (`name` + `description`)
2. Add utility scripts under `scripts/`
3. Update this README's skill table
4. Submit a PR

## License

MIT
