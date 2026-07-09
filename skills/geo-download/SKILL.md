---
name: geo-download
description: Download NCBI GEO datasets, SOFT metadata files, and expression matrices. Use when the user wants to download GEO data by accession ID (GSE/GSM), fetch SOFT format metadata files, batch download multiple GEO datasets, generate GEO download URLs, or list previously downloaded GEO data.
argument-hint: <geo-id> [command]
---

# GEO Download Tool

Download NCBI GEO datasets and metadata via command-line subcommands.

## Usage

Run the bundled script with a subcommand:

```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py <command> [args]
```

## Commands

### download — Download GEO data via GEOparse
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py download GSE201425
```
Parses the GEO record and returns title, summary, sample count, and first 5 sample IDs.

### soft — Download SOFT metadata file
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py soft GSE201425 [-o /path/to/output]
```
Downloads the `.soft.gz` (GSE) or `.soft` (GSM/GPL/GDS) metadata file containing full experiment metadata.

### url — Download from a GEO download URL
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py url "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138709&format=file" [-o /path/to/output]
```
Directly downloads data files from a GEO download link.

### batch — Batch download multiple datasets
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py batch GSE12345 GSE67890 [-o /base/dir]
```

### gen-url — Generate a GEO download URL
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py gen-url GSE12345
```

### list — List downloaded GEO data
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py list [-p /path/to/check]
```

## Dependencies

Requires: `requests`, `GEOparse`, `langchain_core`. Install with:
```bash
pip install requests GEOparse langchain-core
```

## Notes

- Default download directory is the system temp directory.
- SOFT files for GSE are gzip-compressed (`.soft.gz`); decompress with `gunzip` before parsing.
- Combine with the `parse-geo-soft` skill to extract metadata from downloaded SOFT files.
