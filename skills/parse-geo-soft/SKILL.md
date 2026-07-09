---
name: parse-geo-soft
description: Parse NCBI GEO SOFT format files into structured CSV metadata tables. Use when the user wants to extract sample metadata from GEO datasets, parse .soft files, convert SOFT to CSV, or organize GEO experiment information (sample titles, characteristics, platforms, etc.) into tabular format.
argument-hint: <soft-file-path> [-o output.csv]
---

# Parse GEO SOFT File

Parse any NCBI GEO SOFT format file and extract sample metadata into a CSV table.

## Usage

Run the bundled script:

```bash
python3 ~/.qoder-cn/skills/parse-geo-soft/scripts/parse_soft.py <input.soft> [-o output.csv]
```

- `input.soft` — path to the SOFT file (plain text, not gzipped; decompress first if needed)
- `-o output.csv` — output path (default: `<input_basename>_metadata.csv` in the same directory)

## Output Columns

The script dynamically discovers all `Sample_characteristics_ch1` keys (e.g. tissue, age, genotype) and creates a column for each.

Fixed columns:
- `series_id`, `study_title` — from Series section
- `platform_id`, `platform_title` — from Platform section
- `Sample_geo_accession`, `Sample_title` — sample identifiers
- Dynamic characteristics columns (tissue, age, genotype, etc.)
- `description`, `source_name`, `organism`, `molecule`
- `extract_protocol`, `data_processing`, `instrument_model`
- `library_selection`, `library_source`, `library_strategy`

## Workflow

1. Download the SOFT file from GEO FTP if not already present:
   ```
   curl -L -o <file>.soft.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/GSE<acc>/soft/GSE<acc>_family.soft.gz"
   gunzip <file>.soft.gz
   ```
2. Run the parser script on the decompressed `.soft` file.
3. Verify the output CSV has the expected number of samples and columns.

## Notes

- Handles multi-channel samples, multiple characteristics, and varying SOFT structures.
- If the SOFT file is gzipped (`.soft.gz`), decompress it first with `gunzip`.
- The script outputs a summary line with sample count and column list to stdout.
