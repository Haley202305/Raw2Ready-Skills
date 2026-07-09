#!/usr/bin/env python3
"""Parse NCBI GEO SOFT files and output a metadata CSV table.

Usage:
    python parse_soft.py <input.soft> [-o output.csv] [--fields field1,field2,...]

The script dynamically extracts all Sample_characteristics fields
(e.g. tissue, age, genotype) as separate columns.
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict


def parse_soft(filepath):
    series_info = OrderedDict()
    platform_info = OrderedDict()
    samples = []
    current_sample = None
    section = None

    with open(filepath, encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n\r")

            if line.startswith("^SERIES"):
                section = "SERIES"
                current_sample = None
                continue
            elif line.startswith("^PLATFORM"):
                section = "PLATFORM"
                current_sample = None
                continue
            elif line.startswith("^SAMPLE"):
                section = "SAMPLE"
                current_sample = {"_characteristics": OrderedDict(), "_relations": []}
                samples.append(current_sample)
                val = line.split(" = ", 1)[1].strip() if " = " in line else ""
                current_sample["geo_accession"] = val
                continue

            if " = " not in line:
                continue

            key, value = line.split(" = ", 1)
            key = key.lstrip("!")

            if section == "SERIES":
                series_info[key] = series_info.get(key, [])
                if isinstance(series_info[key], list):
                    series_info[key].append(value)
                else:
                    series_info[key] = [series_info[key], value]

            elif section == "PLATFORM":
                platform_info[key] = value

            elif section == "SAMPLE" and current_sample is not None:
                if key == "Sample_characteristics_ch1":
                    if ": " in value:
                        ck, cv = value.split(": ", 1)
                        current_sample["_characteristics"][ck.strip()] = cv.strip()
                    else:
                        current_sample["_characteristics"][value] = ""
                elif key == "Sample_relation":
                    current_sample["_relations"].append(value)
                elif key.startswith("Sample_supplementary_file"):
                    current_sample.setdefault("_supplementary_files", []).append(value)
                else:
                    short_key = key.replace("Sample_", "", 1) if key.startswith("Sample_") else key
                    current_sample[short_key] = value

    return series_info, platform_info, samples


def build_rows(series_info, platform_info, samples):
    all_char_keys = []
    for s in samples:
        for k in s["_characteristics"]:
            if k not in all_char_keys:
                all_char_keys.append(k)

    platform_title = platform_info.get("Platform_title", "")
    series_title = _first(series_info.get("Series_title", ""))
    series_id = _first(series_info.get("Series_geo_accession", ""))

    header = ["series_id", "study_title", "platform_id", "platform_title",
              "Sample_geo_accession", "Sample_title"] + all_char_keys + [
              "description", "source_name", "organism", "molecule",
              "extract_protocol", "data_processing", "instrument_model",
              "library_selection", "library_source", "library_strategy"]

    rows = []
    for s in samples:
        row = OrderedDict()
        row["series_id"] = series_id
        row["study_title"] = series_title
        row["platform_id"] = s.get("platform_id", "")
        row["platform_title"] = platform_title
        row["Sample_geo_accession"] = s.get("geo_accession", "")
        row["Sample_title"] = s.get("title", "")

        for ck in all_char_keys:
            row[ck] = s["_characteristics"].get(ck, "")

        row["description"] = s.get("description", "")
        row["source_name"] = s.get("source_name_ch1", "")
        row["organism"] = s.get("organism_ch1", "")
        row["molecule"] = s.get("molecule_ch1", "")
        row["extract_protocol"] = s.get("extract_protocol_ch1", "")
        row["data_processing"] = s.get("data_processing", "")
        row["instrument_model"] = s.get("instrument_model", "")
        row["library_selection"] = s.get("library_selection", "")
        row["library_source"] = s.get("library_source", "")
        row["library_strategy"] = s.get("library_strategy", "")

        rows.append(row)

    return header, rows


def _first(val):
    if isinstance(val, list):
        return val[0] if val else ""
    return val


def main():
    parser = argparse.ArgumentParser(description="Parse GEO SOFT file to CSV")
    parser.add_argument("input", help="Path to the SOFT file (.soft)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output CSV path (default: <input_basename>_metadata.csv)")
    parser.add_argument("--fields", default=None,
                        help="Comma-separated list of extra characteristics keys to include as columns")
    args = parser.parse_args()

    series_info, platform_info, samples = parse_soft(args.input)

    if not samples:
        print("No samples found in the SOFT file.", file=sys.stderr)
        sys.exit(1)

    header, rows = build_rows(series_info, platform_info, samples)

    if args.output is None:
        import os
        base = os.path.splitext(os.path.basename(args.input))[0]
        args.output = f"{base}_metadata.csv"

    with open(args.output, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Parsed {len(samples)} samples -> {args.output}")
    print(f"Columns: {', '.join(header)}")


if __name__ == "__main__":
    main()
