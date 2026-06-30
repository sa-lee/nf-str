#!/usr/bin/env python3
"""
Repeat Expansion Locus Catalogue Generator
Generates locus catalogues for: ExpansionHunter, LongTR, STRkit,
Atarva, Straglr, and TRGT from simple locus definitions.

Input format — tab-delimited (new format):
    #chrom  start  stop  id  gene  reference_motif_reference_orientation
            pathogenic_motif_reference_orientation  pathogenic_min
            inheritance  disease

Legacy input format (auto-detected):
    chr4:101,106,372-101,106,526    AAAG    PPP3CA

Usage:
    python generate_repeat_catalogs.py --input loci.tsv --outdir catalogues/
    python generate_repeat_catalogs.py --input loci.tsv --prefix STRchive-disease-loci-hg38
"""

import argparse
import json
import os
import re
import subprocess
import sys


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

EXPECTED_NEW_COLS = 10


def parse_locus_string(locus_str: str) -> dict:
    """
    Auto-detect and parse either:
      - New TSV format: chrom  start  stop  id  gene  ref_motif  path_motif  path_min  inheritance  disease
      - Legacy format:  chrom:start-end  motif  [name]

    Coordinates in the new format are assumed to be 0-based half-open (BED).
    Coordinates in the legacy format are assumed to be UCSC 1-based display,
    and are converted to 0-based half-open internally.
    """
    locus_str = locus_str.strip()
    if not locus_str or locus_str.startswith("#"):
        return None

    parts = re.split(r"\t", locus_str)
    if len(parts) == 1:
        parts = re.split(r"\s+", locus_str)

    if re.match(r"^chr[\w.]+$", parts[0]) and len(parts) >= EXPECTED_NEW_COLS:
        return _parse_new_format(parts)

    return _parse_legacy_format(parts, locus_str)


def _parse_new_format(parts: list) -> dict:
    """
    Parse the new 10-column TSV format.
    Coordinates are 0-based half-open (BED) — used as-is.
    Trailing disease description may contain spaces (join remaining fields).
    """
    chrom       = parts[0]
    start       = int(parts[1])
    end         = int(parts[2])
    name        = parts[3]
    gene        = parts[4]
    ref_motif   = parts[5].upper()
    path_motif  = parts[6].upper() if parts[6] not in ("None", ".", "") else None
    path_min    = parts[7] if parts[7] not in ("None", ".", "") else None
    inheritance = parts[8] if len(parts) > 8 else None
    disease     = " ".join(parts[9:]) if len(parts) > 9 else None

    return {
        "chrom":            chrom,
        "start":            start,
        "end":              end,
        "motif":            ref_motif,
        "name":             name,
        "gene":             gene,
        "pathogenic_motif": path_motif,
        "pathogenic_min":   path_min,
        "inheritance":      inheritance,
        "disease":          disease,
        "original_str":     "\t".join(parts),
    }


def _parse_legacy_format(parts: list, original: str) -> dict:
    """
    Parse the legacy 'chrom:start-end  motif  [name]' format.
    Input coordinates are UCSC 1-based display; converted to 0-based half-open.
    """
    if len(parts) < 2:
        raise ValueError(f"Expected at least 2 fields (coords + motif), got: '{original}'")

    coord_str = parts[0].replace(",", "")
    motif     = parts[1].upper()
    name      = parts[2] if len(parts) >= 3 else None

    m = re.match(r"^([\w.]+):(\d+)-(\d+)$", coord_str)
    if not m:
        raise ValueError(f"Cannot parse coordinates: '{coord_str}'. Expected chrom:start-end")

    chrom = m.group(1)
    start = int(m.group(2))
    end   = int(m.group(3))
    start_0based = start - 1

    if name is None:
        name = f"{chrom}_{start_0based}_{end}"

    return {
        "chrom":            chrom,
        "start":            start_0based,
        "end":              end,
        "motif":            motif,
        "name":             name,
        "gene":             None,
        "pathogenic_motif": None,
        "pathogenic_min":   None,
        "inheritance":      None,
        "disease":          None,
        "original_str":     original,
    }


def load_loci_from_file(filepath: str) -> list:
    loci = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            locus = parse_locus_string(line)
            if locus:
                loci.append(locus)
    return loci


def sort_loci(loci: list) -> list:
    """Sort loci by chromosome (natural sort) then start position."""
    def chrom_key(locus):
        chrom = locus["chrom"].lstrip("chr")
        chrom_map = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        try:
            return (int(chrom_map.get(chrom, chrom)), locus["start"])
        except (ValueError, TypeError):
            return (26, locus["start"])
    return sorted(loci, key=chrom_key)


def make_filename(outdir: str, prefix: str, toolname: str, ext: str = "bed") -> str:
    """
    Build an output path like:
        outdir/prefix_toolname.ext   (if prefix given)
        outdir/toolname.ext          (if no prefix)
    """
    basename = f"{prefix}_{toolname}.{ext}" if prefix else f"{toolname}.{ext}"
    return os.path.join(outdir, basename)


# ---------------------------------------------------------------------------
# Catalogue writers
# ---------------------------------------------------------------------------

def write_expansionhunter(loci: list, outdir: str, prefix: str):
    """
    ExpansionHunter variant catalog — JSON.
    Coordinates: 0-based half-open.
    Ref: https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md
    """
    catalog = []
    for locus in loci:
        entry = {
            "LocusId":         locus["name"],
            "LocusStructure":  f"({locus['motif']})*",
            "ReferenceRegion": f"{locus['chrom']}:{locus['start']}-{locus['end']}",
            "VariantType":     "Repeat",
        }
        if locus.get("pathogenic_motif") and locus["pathogenic_motif"] != locus["motif"]:
            entry["OfftargetRegions"] = []
            entry["VariantType"] = "RareRepeat"
        catalog.append(entry)

    outpath = make_filename(outdir, prefix, "expansionhunter", "json")
    with open(outpath, "w") as fh:
        json.dump(catalog, fh, indent=4)
    print(f"  [ExpansionHunter]  -> {outpath}  ({len(loci)} loci)")


def write_trgt(loci: list, outdir: str, prefix: str):
    """
    TRGT repeat definitions — BED4 with structured INFO in column 4.
    Coordinates: 0-based half-open (standard BED).
    Ref: https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md
    """
    outpath = make_filename(outdir, prefix, "trgt")
    with open(outpath, "w") as fh:
        for locus in loci:
            motifs = locus["motif"]
            if locus.get("pathogenic_motif") and locus["pathogenic_motif"] != locus["motif"]:
                motifs = f"{locus['motif']},{locus['pathogenic_motif']}"
            info = f"ID={locus['name']};MOTIFS={motifs};STRUC=<TR>"
            fh.write(f"{locus['chrom']}\t{locus['start']}\t{locus['end']}\t{info}\n")
    print(f"  [TRGT]             -> {outpath}  ({len(loci)} loci)")


def write_longtr(loci: list, outdir: str, prefix: str):
    """
    LongTR region BED file (v1.2+ format).
    Columns: chrom | start(1-based) | end | motif | name
    Ref: https://github.com/gymrek-lab/LongTR
    """
    outpath = make_filename(outdir, prefix, "longtr")
    with open(outpath, "w") as fh:
        for locus in loci:
            start_1based = locus["start"] + 1
            fh.write(
                f"{locus['chrom']}\t{start_1based}\t{locus['end']}\t"
                f"{locus['motif']}\t{locus['name']}\n"
            )
    print(f"  [LongTR]           -> {outpath}  ({len(loci)} loci)")


def write_strkit(loci: list, outdir: str, prefix: str):
    """
    STRkit locus catalog — 5-column BED (chrom, start, end, name, motif).
    Coordinates: 0-based half-open. File must be sorted by position.
    Ref: https://github.com/davidlougheed/strkit/blob/master/docs/caller_catalog.md
    """
    outpath = make_filename(outdir, prefix, "strkit")
    with open(outpath, "w") as fh:
        for locus in sort_loci(loci):
            fh.write(
                f"{locus['chrom']}\t{locus['start']}\t{locus['end']}\t"
                f"{locus['name']}\t{locus['motif']}\n"
            )
    print(f"  [STRkit]           -> {outpath}  ({len(loci)} loci)")


def write_straglr(loci: list, outdir: str, prefix: str):
    """
    Straglr locus catalog — 4-column BED (chrom, start, end, motif).
    Coordinates: 0-based half-open (standard BED).
    Ref: https://github.com/bcgsc/straglr
    """
    outpath = make_filename(outdir, prefix, "straglr")
    with open(outpath, "w") as fh:
        for locus in loci:
            fh.write(
                f"{locus['chrom']}\t{locus['start']}\t{locus['end']}\t-\n") #f"{locus['motif']}\n"

    print(f"  [Straglr]          -> {outpath}  ({len(loci)} loci)")


def write_atarva(loci: list, outdir: str, prefix: str):
    """
    ATaRVa locus catalog — 6-column BED, sorted, bgzipped and tabix indexed.
    Columns: chrom | start | end | motif | motif_length | id
    Coordinates: 0-based half-open (standard BED).
    Ref: https://github.com/SowpatiLab/ATaRVa
    """
    bed_path = make_filename(outdir, prefix, "atarva")
    bgz_path = bed_path + ".gz"
    tbi_path = bgz_path + ".tbi"

    with open(bed_path, "w") as fh:
        for locus in sort_loci(loci):
            motif_len = len(locus["motif"])
            fh.write(
                f"{locus['chrom']}\t{locus['start']}\t{locus['end']}\t"
                f"{locus['motif']}\t{motif_len}\t{locus['name']}\n"
            )

    try:
        subprocess.run(["bgzip", "-f", bed_path], check=True, capture_output=True)
        print(f"  [ATaRVa]           -> {bgz_path}  ({len(loci)} loci, bgzipped)")
    except FileNotFoundError:
        print(
            f"  [ATaRVa] WARNING: bgzip not found. Uncompressed BED written to {bed_path}.\n"
            f"           Install htslib and run manually:\n"
            f"             bgzip -f {bed_path}\n"
            f"             tabix -p bed {bgz_path}"
        )
        return
    except subprocess.CalledProcessError as e:
        print(f"  [ATaRVa] ERROR during bgzip: {e.stderr.decode().strip()}")
        return

    try:
        subprocess.run(["tabix", "-p", "bed", bgz_path], check=True, capture_output=True)
        print(f"  [ATaRVa]           -> {tbi_path}  (tabix index)")
    except FileNotFoundError:
        print(f"  [ATaRVa] WARNING: tabix not found. Run manually: tabix -p bed {bgz_path}")
    except subprocess.CalledProcessError as e:
        print(f"  [ATaRVa] ERROR during tabix: {e.stderr.decode().strip()}")


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def write_summary(loci: list, outdir: str, prefix: str):
    """Write a human-readable TSV summary of all parsed loci."""
    outpath = make_filename(outdir, prefix, "summary", "tsv")
    with open(outpath, "w") as fh:
        fh.write(
            "name\tchrom\tstart_0based\tend\tgene\tmotif\tpathogenic_motif\t"
            "pathogenic_min\tinheritance\tdisease\tmotif_len\tlocus_size_bp\n"
        )
        for locus in loci:
            size = locus["end"] - locus["start"]
            fh.write(
                f"{locus['name']}\t{locus['chrom']}\t{locus['start']}\t{locus['end']}\t"
                f"{locus.get('gene') or ''}\t{locus['motif']}\t"
                f"{locus.get('pathogenic_motif') or ''}\t"
                f"{locus.get('pathogenic_min') or ''}\t"
                f"{locus.get('inheritance') or ''}\t"
                f"{locus.get('disease') or ''}\t"
                f"{len(locus['motif'])}\t{size}\n"
            )
    print(f"  [Summary]          -> {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate repeat expansion locus catalogues for multiple tools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage — toolname-only filenames:
  python generate_repeat_catalogs.py --input loci.tsv --outdir catalogues/

  # With a prefix — e.g. STRchive-disease-loci-hg38_longtr.bed:
  python generate_repeat_catalogs.py \\
      --input loci.tsv \\
      --prefix STRchive-disease-loci-hg38 \\
      --outdir catalogues/

  # Single locus from command line:
  python generate_repeat_catalogs.py \\
      --locus "chrX 147582151 147582229 FMR1_FMR1 FMR1 CGG CGG 55 XL Fragile X syndrome" \\
      --prefix my-loci-hg38

  # Only generate specific tool catalogues:
  python generate_repeat_catalogs.py \\
      --input loci.tsv --tools trgt strkit longtr --prefix STRchive-disease-loci-hg38

Output filenames:
  Without prefix:  longtr.bed, trgt.bed, strkit.bed, straglr.bed,
                   atarva.bed.gz, expansionhunter.json, summary.tsv
  With prefix:     STRchive-disease-loci-hg38_longtr.bed, etc.
        """,
    )
    parser.add_argument(
        "--input", "-i",
        metavar="FILE",
        help="Input file with one locus per line",
    )
    parser.add_argument(
        "--locus", "-l",
        metavar="STR",
        action="append",
        default=[],
        help="Single locus string. Can be specified multiple times.",
    )
    parser.add_argument(
        "--outdir", "-o",
        metavar="DIR",
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--prefix", "-p",
        metavar="STR",
        default="",
        help='Filename prefix, e.g. "STRchive-disease-loci-hg38" -> STRchive-disease-loci-hg38_longtr.bed',
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        choices=["expansionhunter", "trgt", "longtr", "strkit", "straglr", "atarva", "all"],
        default=["all"],
        help="Which tool catalogues to generate (default: all)",
    )
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Skip writing the summary.tsv file",
    )

    args = parser.parse_args()

    # --- Collect loci ---
    loci = []

    if args.input:
        file_loci = load_loci_from_file(args.input)
        loci.extend(file_loci)
        print(f"Loaded {len(file_loci)} loci from {args.input}")

    for locus_str in args.locus:
        locus = parse_locus_string(locus_str)
        if locus:
            loci.append(locus)

    if not loci:
        parser.error("No loci provided. Use --input FILE or one or more --locus arguments.")

    # --- Check for duplicate names ---
    names = [l["name"] for l in loci]
    if len(names) != len(set(names)):
        print("WARNING: Duplicate locus names detected — consider using unique names.")

    # --- Create output directory ---
    os.makedirs(args.outdir, exist_ok=True)

    prefix_display = f"{args.prefix}_" if args.prefix else ""
    print(f"\nGenerating catalogues for {len(loci)} loci -> {args.outdir}/{prefix_display}*.bed\n")

    # --- Determine which tools to run ---
    tools = set(args.tools)
    run_all = "all" in tools

    tool_map = {
        "expansionhunter": write_expansionhunter,
        "trgt":            write_trgt,
        "longtr":          write_longtr,
        "strkit":          write_strkit,
        "straglr":         write_straglr,
        "atarva":          write_atarva,
    }

    for tool_name, writer_fn in tool_map.items():
        if run_all or tool_name in tools:
            writer_fn(loci, args.outdir, args.prefix)

    if not args.no_summary:
        write_summary(loci, args.outdir, args.prefix)

    print(f"\nDone. {len(loci)} loci processed.")


if __name__ == "__main__":
    main()