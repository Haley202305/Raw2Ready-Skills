#!/usr/bin/env python3
"""CLI wrapper for GEO download tools."""
import sys
import os
import argparse

sys.path.insert(0, "/public/home/jihongyi/download_data/tools")

from geo_download_tool import (
    download_geo_series,
    download_geo_soft,
    download_geo_from_url,
    download_multiple_geo_datasets,
    get_geo_download_url,
    list_downloaded_geo_data,
)


def main():
    parser = argparse.ArgumentParser(description="GEO Data Download Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    p1 = subparsers.add_parser("download", help="Download GEO series/sample data via GEOparse")
    p1.add_argument("geo_id", help="GEO ID (e.g. GSE201425, GSM12345)")

    p2 = subparsers.add_parser("soft", help="Download GEO SOFT metadata file")
    p2.add_argument("geo_id", help="GEO ID (e.g. GSE201425)")
    p2.add_argument("-o", "--output-dir", default=None, help="Output directory")

    p3 = subparsers.add_parser("url", help="Download from a GEO download URL")
    p3.add_argument("url", help="GEO download URL")
    p3.add_argument("-o", "--output-dir", default=None, help="Output directory")

    p4 = subparsers.add_parser("batch", help="Batch download multiple GEO datasets")
    p4.add_argument("geo_ids", nargs="+", help="GEO IDs (e.g. GSE12345 GSE67890)")
    p4.add_argument("-o", "--output-dir", default=None, help="Base output directory")

    p5 = subparsers.add_parser("gen-url", help="Generate GEO download URL")
    p5.add_argument("geo_id", help="GEO ID")
    p5.add_argument("-f", "--format", default="file", help="Download format (default: file)")

    p6 = subparsers.add_parser("list", help="List downloaded GEO data")
    p6.add_argument("-p", "--path", default=None, help="Base path to check")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "download":
        result = download_geo_series.invoke({"geo_series_id": args.geo_id})
    elif args.command == "soft":
        result = download_geo_soft.invoke({"geo_id": args.geo_id, "output_dir": args.output_dir})
    elif args.command == "url":
        result = download_geo_from_url.invoke({"download_url": args.url, "output_dir": args.output_dir})
    elif args.command == "batch":
        result = download_multiple_geo_datasets.invoke({"geo_ids": args.geo_ids, "output_base_dir": args.output_dir})
    elif args.command == "gen-url":
        result = get_geo_download_url.invoke({"geo_id": args.geo_id, "format_type": args.format})
    elif args.command == "list":
        result = list_downloaded_geo_data.invoke({"base_path": args.path})

    print(result)


if __name__ == "__main__":
    main()
