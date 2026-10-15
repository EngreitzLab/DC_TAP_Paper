#!/usr/bin/env python3

import csv
import math
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = ROOT / "results" / "2026-05-11_SCEPTRE" / "runs"
TRACKS_DIR = ROOT / "results" / "2026-05-11_SCEPTRE" / "tracks"
DATA_DIR = ROOT / "data"

TABLE_S3 = DATA_DIR / "Table_S3_Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements_250908_append.tsv"
OLD_BEDGRAPH = DATA_DIR / "2025-11-26-dctapseq-wtc11-sig-grna-pam-pct-eff-size.bedGraph"

DATE_PREFIX = "2026-05-13"


TARGET_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_target(target):
    match = TARGET_RE.match(target)
    if not match:
        raise ValueError(f"Could not parse target coordinate: {target}")
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def clean_value(value):
    if value is None:
        return ""
    value = value.strip()
    if value.upper() == "NA":
        return ""
    return value


def parse_float(value):
    value = clean_value(value)
    if not value:
        return math.nan
    return float(value)


def parse_int(value):
    value = clean_value(value)
    if not value:
        return None
    return int(float(value))


def base_gene_id(gene_id):
    return clean_value(gene_id).split(".")[0]


def split_guide_ids(guide_ids):
    guide_ids = clean_value(guide_ids)
    if not guide_ids:
        return []
    guides = []
    for guide_id in guide_ids.split(","):
        guide_id = guide_id.strip().strip('"')
        if guide_id:
            guides.append(guide_id.replace("_", "-"))
    return guides


def chrom_key(chrom):
    if chrom.startswith("chr"):
        suffix = chrom[3:]
    else:
        suffix = chrom
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if suffix.isdigit():
        return 0, int(suffix)
    return 1, special.get(suffix, 1000), suffix


def row_pct(row):
    return 100.0 * (float(row["fold_change"]) - 1.0)


def is_true(value):
    return clean_value(value).upper() == "TRUE"


def load_table_s3():
    metadata_by_pair = {}
    guides_by_target = {}

    for row in read_tsv(TABLE_S3):
        if clean_value(row.get("cell_type")) != "WTC11":
            continue

        gene_id = base_gene_id(row.get("gene_id"))
        target = clean_value(row.get("intended_target_name_hg38"))
        if not gene_id or not target:
            continue

        key = (gene_id, target)
        metadata_by_pair.setdefault(key, row)

        guides = split_guide_ids(row.get("guide_ids"))
        if guides and target not in guides_by_target:
            guides_by_target[target] = guides

    return metadata_by_pair, guides_by_target


def load_old_pam_intervals():
    intervals = []
    with OLD_BEDGRAPH.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            chrom, start, end, _value = line.split("\t")[:4]
            if start == "." or end == ".":
                continue
            intervals.append((chrom, int(float(start)), int(float(end))))
    return intervals


def build_guide_coordinate_lookup(guides_by_target):
    old_intervals = load_old_pam_intervals()
    guide_coords = {}
    warnings = []

    for target, guides in sorted(guides_by_target.items()):
        chrom, target_start, target_end = parse_target(target)
        unique_intervals = sorted({
            (interval_chrom, interval_start, interval_end)
            for interval_chrom, interval_start, interval_end in old_intervals
            if (
                interval_chrom == chrom
                and interval_start >= target_start - 1
                and interval_end <= target_end
            )
        }, key=lambda x: (x[1], x[2]))

        if len(unique_intervals) == len(guides):
            for guide_id, interval in zip(guides, unique_intervals):
                guide_coords[guide_id] = interval
        elif unique_intervals:
            warnings.append({
                "target": target,
                "n_guides": str(len(guides)),
                "n_unique_old_pam_intervals": str(len(unique_intervals)),
                "warning": "old_pam_interval_count_does_not_match_guide_count",
            })
        else:
            warnings.append({
                "target": target,
                "n_guides": str(len(guides)),
                "n_unique_old_pam_intervals": "0",
                "warning": "no_old_pam_intervals_for_target",
            })

    return guide_coords, warnings


def target_interval(row):
    chrom, start, end = parse_target(row["grna_target"])
    return chrom, start - 1, end


def arc_interval(metadata, row):
    chrom, target_start, target_end = parse_target(row["grna_target"])
    target_center = (target_start + target_end) / 2.0

    if not metadata:
        return chrom, target_start - 1, target_end

    distance = parse_float(metadata.get("distance_to_abc_canonical_TSS"))
    tss_start = parse_int(metadata.get("startTSS_hg38"))
    arc_chrom = clean_value(metadata.get("chrTSS_hg38")) or chrom

    if math.isnan(distance) or tss_start is None:
        return arc_chrom, target_start - 1, target_end

    if distance == 0:
        direction = 0
    elif tss_start < target_center:
        direction = -1
    else:
        direction = 1

    abc_tss = target_center + direction * distance
    start = max(0, int(math.floor(abc_tss)) - 251)
    return arc_chrom, start, start + 501


def pair_name(metadata, row):
    if metadata:
        name = clean_value(metadata.get("element_gene_pair_identifier_hg19"))
        if not name:
            name = clean_value(metadata.get("element_gene_pair_identifier_hg38"))
        if name:
            return f"{name}:."
        gene_symbol = clean_value(metadata.get("gene_symbol"))
        if gene_symbol:
            return f"{gene_symbol}|{row['grna_target']}:."
    return f"{base_gene_id(row['response_id'])}|{row['grna_target']}:."


def color_for_pct(pct):
    return "255,0,0" if pct < 0 else "0,255,0"


def write_bedgraph(rows, path):
    rows = sorted(rows, key=lambda row: (chrom_key(row[0]), row[1], row[2], row[3]))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_bedpe(rows, path):
    rows = sorted(rows, key=lambda row: (chrom_key(row[0]), row[1], row[2], chrom_key(row[3]), row[4], row[5], row[6]))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chr1", "x1", "x2", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2", "color"])
        writer.writerows(rows)


def collapse_by_abs_pct(rows, key_fields):
    best = {}
    for row in rows:
        key = tuple(row[field] for field in key_fields)
        pct = row["_pct_change_effect_size"]
        if key not in best or abs(pct) > abs(best[key]["_pct_change_effect_size"]):
            best[key] = row
    return list(best.values())


def make_bedgraph_rows(rows, guide_coords):
    out = []
    fallback_rows = []
    for row in rows:
        guide_id = row.get("grna_id")
        if guide_id and guide_id in guide_coords:
            chrom, start, end = guide_coords[guide_id]
        else:
            chrom, start, end = target_interval(row)
            fallback_rows.append({
                "run": row["_run"],
                "response_id": row["response_id"],
                "grna_id": row.get("grna_id", ""),
                "grna_target": row["grna_target"],
                "pct_change_effect_size": row["_pct_change_effect_size"],
            })
        out.append([chrom, start, end, row["_pct_change_effect_size"]])
    return out, fallback_rows


def make_bedpe_rows(rows, metadata_by_pair):
    out = []
    for row in rows:
        key = (base_gene_id(row["response_id"]), row["grna_target"])
        metadata = metadata_by_pair.get(key)
        chrom1, start1, end1 = target_interval(row)
        chrom2, start2, end2 = arc_interval(metadata, row)
        pct = row["_pct_change_effect_size"]
        out.append([
            chrom1,
            start1,
            end1,
            chrom2,
            start2,
            end2,
            pair_name(metadata, row),
            1000,
            "+",
            "+",
            color_for_pct(pct),
        ])
    return out


def read_run_results(path):
    run = path.parent.name
    rows = []
    for row in read_tsv(path):
        if not is_true(row.get("significant")):
            continue
        if "pass_qc" in row and not is_true(row.get("pass_qc")):
            continue
        row = dict(row)
        row["_run"] = run
        row["_pct_change_effect_size"] = row_pct(row)
        rows.append(row)
    return rows


def write_manifest(rows):
    path = TRACKS_DIR / "track_manifest.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "file",
                "track_type",
                "source",
                "n_rows",
                "n_target_interval_fallback_rows",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_warnings(warnings):
    path = TRACKS_DIR / "guide_coordinate_mapping_warnings.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["target", "n_guides", "n_unique_old_pam_intervals", "warning"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(warnings)


def write_fallback_rows(rows):
    path = TRACKS_DIR / "bedgraph_target_interval_fallback_rows.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run", "response_id", "grna_id", "grna_target", "pct_change_effect_size"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    TRACKS_DIR.mkdir(parents=True, exist_ok=True)

    metadata_by_pair, guides_by_target = load_table_s3()
    guide_coords, warnings = build_guide_coordinate_lookup(guides_by_target)
    write_warnings(warnings)

    manifest = []
    all_singleton = []
    all_union = []
    all_fallback_rows = []

    for result_path in sorted(RUNS_DIR.glob("*/sceptre_discovery_results.tsv")):
        run = result_path.parent.name
        rows = read_run_results(result_path)
        if not rows:
            continue

        is_singleton = "grna_id" in rows[0]
        if is_singleton:
            all_singleton.extend(rows)
            bedgraph_rows, fallback_rows = make_bedgraph_rows(rows, guide_coords)
            all_fallback_rows.extend(fallback_rows)
            bedgraph_path = TRACKS_DIR / f"{run}_sig_grna_pam_pct_eff_size.bedGraph"
            write_bedgraph(bedgraph_rows, bedgraph_path)
            manifest.append({
                "file": str(bedgraph_path.relative_to(ROOT)),
                "track_type": "bedGraph",
                "source": run,
                "n_rows": str(len(bedgraph_rows)),
                "n_target_interval_fallback_rows": str(len(fallback_rows)),
            })
        else:
            all_union.extend(rows)

        arc_input = collapse_by_abs_pct(rows, ["response_id", "grna_target"])
        bedpe_rows = make_bedpe_rows(arc_input, metadata_by_pair)
        bedpe_path = TRACKS_DIR / f"{run}_signif_arcs.bedpe"
        write_bedpe(bedpe_rows, bedpe_path)
        manifest.append({
            "file": str(bedpe_path.relative_to(ROOT)),
            "track_type": "bedpe",
            "source": run,
            "n_rows": str(len(bedpe_rows)),
            "n_target_interval_fallback_rows": "0",
        })

    if all_singleton:
        singleton_best = collapse_by_abs_pct(all_singleton, ["response_id", "grna_target", "grna_id"])
        bedgraph_rows, fallback_rows = make_bedgraph_rows(singleton_best, guide_coords)
        for bedgraph_path in [
            TRACKS_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-sig-grna-pam-pct-eff-size.bedGraph",
            DATA_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-sig-grna-pam-pct-eff-size.bedGraph",
        ]:
            write_bedgraph(bedgraph_rows, bedgraph_path)
        manifest.append({
            "file": str((TRACKS_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-sig-grna-pam-pct-eff-size.bedGraph").relative_to(ROOT)),
            "track_type": "bedGraph",
            "source": "singleton_any_best_abs_effect_per_gene_target_guide",
            "n_rows": str(len(bedgraph_rows)),
            "n_target_interval_fallback_rows": str(len(fallback_rows)),
        })

        singleton_arc_input = collapse_by_abs_pct(all_singleton, ["response_id", "grna_target"])
        singleton_bedpe_rows = make_bedpe_rows(singleton_arc_input, metadata_by_pair)
        for bedpe_path in [
            TRACKS_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-signif-arcs.bedpe",
            DATA_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-signif-arcs.bedpe",
        ]:
            write_bedpe(singleton_bedpe_rows, bedpe_path)
        manifest.append({
            "file": str((TRACKS_DIR / f"{DATE_PREFIX}-sceptre-singleton-any-signif-arcs.bedpe").relative_to(ROOT)),
            "track_type": "bedpe",
            "source": "singleton_any_best_abs_effect_per_gene_target",
            "n_rows": str(len(singleton_bedpe_rows)),
            "n_target_interval_fallback_rows": "0",
        })

    if all_union:
        union_arc_input = collapse_by_abs_pct(all_union, ["response_id", "grna_target"])
        union_bedpe_rows = make_bedpe_rows(union_arc_input, metadata_by_pair)
        for bedpe_path in [
            TRACKS_DIR / f"{DATE_PREFIX}-sceptre-union-any-signif-arcs.bedpe",
            DATA_DIR / f"{DATE_PREFIX}-sceptre-union-any-signif-arcs.bedpe",
        ]:
            write_bedpe(union_bedpe_rows, bedpe_path)
        manifest.append({
            "file": str((TRACKS_DIR / f"{DATE_PREFIX}-sceptre-union-any-signif-arcs.bedpe").relative_to(ROOT)),
            "track_type": "bedpe",
            "source": "union_any_best_abs_effect_per_gene_target",
            "n_rows": str(len(union_bedpe_rows)),
            "n_target_interval_fallback_rows": "0",
        })

    write_manifest(manifest)
    write_fallback_rows(all_fallback_rows)


if __name__ == "__main__":
    main()
