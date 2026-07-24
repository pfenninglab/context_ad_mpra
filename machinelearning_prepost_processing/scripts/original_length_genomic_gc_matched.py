#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline for ORIGINAL-LENGTH positive peaks + genomic GC-matched negatives.

Key changes relative to the old script:
  1) Positive peaks are NOT expanded around summit.
  2) Train/validation/test splits are made directly from the original peak coordinates.
  3) Negative sequences are generated from SHUFFLED GENOMIC LOCATIONS that preserve
     chromosome and interval length, then matched to each positive by GC content.
  4) Overlaps with all_possible peaks are removed using bedtools shuffle -excl.

Outputs:
  POS:
    {bed_out_dir}/{name}.original_length.bed
    {bed_out_dir}/{name}.narrowPeak.train.bed (+ validation/test)
    {fasta_out_dir}/{name}.original_length.fasta
    {fasta_out_dir}/{name}.narrowPeak.train.bed.fasta (+ validation/test)

  NEG (real genomic coordinates):
    {neg_out_dir}/{name}.genomicGCmatch.original_length.fasta
    {neg_out_dir}/{name}.genomicGCmatch.original_length.train.fasta (+ validation/test)
    {bed_out_dir}/{name}.narrowPeak.train.Neg.bed (+ validation/test)
    {bed_out_dir}/{name}.narrowPeak.train.withNeg.bed (+ validation/test)
    {bed_out_dir}/{name}.narrowPeak.train.Neg_Part1.bed ... Neg_Part4.bed

  POSNEG FASTA:
    {posneg_out_dir}/{name}.narrowPeak.train.posneg.fasta (+ validation/test)

Requirements:
  - Python 3.7+
  - pandas, numpy
  - bedtools in PATH

Notes:
  - The genomic negatives are chosen greedily from an oversampled bedtools shuffle pool.
  - If too few negatives are found, increase oversample_factor.
  - exclude_padding_bp defaults to 0 because positives are no longer expanded.
"""

import os
import gzip
import glob
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


# =========================================================
# helpers
# =========================================================
def require_file(p: str, label: str = ""):
    pth = Path(p)
    if not pth.is_file():
        raise FileNotFoundError("{} missing file: {}".format(label, pth))


def require_dir(p: str, label: str = ""):
    pth = Path(p)
    if not pth.is_dir():
        raise FileNotFoundError("{} missing dir: {}".format(label, pth))


def ensure_dir(p: str):
    Path(p).mkdir(parents=True, exist_ok=True)


def which_or_raise(exe: str):
    if shutil.which(exe) is None:
        raise EnvironmentError("Missing executable in PATH: {}".format(exe))


def run(cmd: str):
    print("[RUN] {}".format(cmd))
    subprocess.run(cmd, shell=True, check=True)


def gunzip_to_temp(gz_path: str, out_path: str):
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def maybe_decompress(path_in: str, work_dir: str) -> str:
    p = Path(path_in)
    if p.suffix == ".gz":
        out_path = str(Path(work_dir) / p.name.replace(".gz", ""))
        gunzip_to_temp(str(p), out_path)
        return out_path
    return str(p)


def write_fasta_from_df(df: pd.DataFrame, fasta_path: str, header_col: str = "header", seq_col: str = "sequence"):
    ensure_dir(str(Path(fasta_path).parent))
    with open(fasta_path, "w") as fout:
        for _, row in df.iterrows():
            fout.write(">{header}\n{seq}\n".format(header=row[header_col], seq=row[seq_col]))


def gc_fraction(seq: str) -> float:
    seq = str(seq).upper()
    valid = [b for b in seq if b in set(["A", "C", "G", "T"])]
    if len(valid) == 0:
        return np.nan
    gc = sum([1 for b in valid if b in set(["G", "C"])])
    return float(gc) / float(len(valid))


def read_bedtools_tab(tab_file: str) -> pd.DataFrame:
    rows = []
    with open(tab_file, "r") as fin:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            name, seq = line.split("\t", 1)
            rows.append((name, seq))
    return pd.DataFrame(rows, columns=["name", "sequence"])


def bed_to_fasta_df(bed_file: str, genome_fa: str, work_dir: str, prefix: str) -> pd.DataFrame:
    tab_out = str(Path(work_dir) / (prefix + ".bedtools.tab"))
    run("bedtools getfasta -fi {genome} -bed {bed} -name -tab > {out}".format(
        genome=genome_fa,
        bed=bed_file,
        out=tab_out,
    ))
    return read_bedtools_tab(tab_out)


def write_bed(df: pd.DataFrame, bed_path: str, cols: List[str]):
    ensure_dir(str(Path(bed_path).parent))
    df.loc[:, cols].to_csv(bed_path, sep="\t", header=False, index=False)


def split_by_chrom(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    val = df[df["chrom"] == "chr4"].copy()
    test = df[df["chrom"].isin(["chr8", "chr9"])].copy()
    train = df[~df["chrom"].isin(["chr4", "chr8", "chr9"])].copy()
    return train, val, test


# =========================================================
# positive peaks (original length)
# =========================================================
def load_original_narrowpeak(narrowpeak_in: str) -> pd.DataFrame:
    df = pd.read_csv(narrowpeak_in, sep="\t", header=None)
    if df.shape[1] < 10:
        raise ValueError("{} has {} cols; expected >=10 for narrowPeak.".format(narrowpeak_in, df.shape[1]))

    df = df.iloc[:, :10].copy()
    df.columns = [
        "chrom", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"
    ]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["peak_id"] = ["pos_{:07d}".format(i) for i in range(df.shape[0])]
    return df


def write_positive_outputs(df_pos: pd.DataFrame, cfg: Dict[str, str]) -> Dict[str, str]:
    name = cfg["name"]
    genome_fa = cfg["genome_fa"]
    work_dir = cfg["work_dir"]
    bed_out_dir = cfg["bed_out_dir"]
    fasta_out_dir = cfg["fasta_out_dir"]

    ensure_dir(bed_out_dir)
    ensure_dir(fasta_out_dir)

    pos_full_bed = str(Path(bed_out_dir) / (name + ".original_length.bed"))
    write_bed(df_pos, pos_full_bed, [
        "chrom", "start", "end", "peak_id", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"
    ])

    train_df, val_df, test_df = split_by_chrom(df_pos)

    pos_train_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.train.bed"))
    pos_val_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.validation.bed"))
    pos_test_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.test.bed"))

    for dfi, out_bed in [(train_df, pos_train_bed), (val_df, pos_val_bed), (test_df, pos_test_bed)]:
        write_bed(dfi, out_bed, [
            "chrom", "start", "end", "peak_id", "score", "strand",
            "signalValue", "pValue", "qValue", "summit"
        ])

    pos_full_seq = bed_to_fasta_df(pos_full_bed, genome_fa, work_dir, name + ".pos_full")
    pos_full = df_pos.merge(pos_full_seq, left_on="peak_id", right_on="name", how="left")
    pos_full["header"] = pos_full["peak_id"] + "|" + pos_full["chrom"] + ":" + pos_full["start"].astype(str) + "-" + pos_full["end"].astype(str)
    pos_full["gc"] = pos_full["sequence"].map(gc_fraction)

    pos_full_fa = str(Path(fasta_out_dir) / (name + ".original_length.fasta"))
    write_fasta_from_df(pos_full, pos_full_fa)

    out = {
        "pos_full_bed": pos_full_bed,
        "pos_full_fa": pos_full_fa,
        "pos_train_bed": pos_train_bed,
        "pos_val_bed": pos_val_bed,
        "pos_test_bed": pos_test_bed,
        "pos_train_fa": str(Path(fasta_out_dir) / (Path(pos_train_bed).name + ".fasta")),
        "pos_val_fa": str(Path(fasta_out_dir) / (Path(pos_val_bed).name + ".fasta")),
        "pos_test_fa": str(Path(fasta_out_dir) / (Path(pos_test_bed).name + ".fasta")),
    }

    for subset_df, out_fa in [(train_df, out["pos_train_fa"]), (val_df, out["pos_val_fa"]), (test_df, out["pos_test_fa"])]:
        sub = pos_full[pos_full["peak_id"].isin(subset_df["peak_id"])].copy()
        write_fasta_from_df(sub, out_fa)

    return out, pos_full, train_df, val_df, test_df


# =========================================================
# exclusion bed for no-overlap genomic negatives
# =========================================================
def build_exclusion_bed(all_possible_glob: str, work_dir: str, prefix: str, padding_bp: int = 0, merge_distance: int = 200) -> str:
    files = sorted(glob.glob(all_possible_glob))
    if len(files) == 0:
        raise FileNotFoundError("No files matched: {}".format(all_possible_glob))

    cated = str(Path(work_dir) / (prefix + ".all_possible_peaks.narrowPeak"))
    merged = str(Path(work_dir) / (prefix + ".all_possible_peaks.merged.bed"))
    expanded = str(Path(work_dir) / (prefix + ".all_possible_peaks.merged.pad{}.bed".format(int(padding_bp))))

    with open(cated, "wb") as out:
        for fp in files:
            p = Path(fp)
            if p.suffix == ".gz":
                with gzip.open(fp, "rb") as f_in:
                    shutil.copyfileobj(f_in, out)
            else:
                with open(fp, "rb") as f_in:
                    shutil.copyfileobj(f_in, out)

    run("bedtools sort -i {cated} | bedtools merge -d {d} -i - > {merged}".format(
        cated=cated,
        d=int(merge_distance),
        merged=merged,
    ))

    df = pd.read_csv(merged, sep="\t", header=None, names=["chrom", "start", "end"])
    df["start"] = (df["start"].astype(int) - int(padding_bp)).clip(lower=0)
    df["end"] = df["end"].astype(int) + int(padding_bp)
    df[["chrom", "start", "end"]].to_csv(expanded, sep="\t", header=False, index=False)
    return expanded


# =========================================================
# genomic GC-matched negatives
# =========================================================
def build_oversampled_template_bed(df_pos: pd.DataFrame, out_bed: str, oversample_factor: int):
    rows = []
    for _, row in df_pos.iterrows():
        for rep in range(int(oversample_factor)):
            rows.append([
                row["chrom"], int(row["start"]), int(row["end"]),
                "{}__rep{:03d}".format(row["peak_id"], rep), 0, "+"
            ])
    out = pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "score", "strand"])
    out.to_csv(out_bed, sep="\t", header=False, index=False)


def parse_shuffle_candidates(shuffled_bed: str, seq_tab: str, df_pos_full: pd.DataFrame) -> pd.DataFrame:
    bed_cols = ["chrom", "start", "end", "candidate_name", "score", "strand"]
    df_bed = pd.read_csv(shuffled_bed, sep="\t", header=None, names=bed_cols)
    df_seq = read_bedtools_tab(seq_tab)
    df = df_bed.merge(df_seq, left_on="candidate_name", right_on="name", how="left")

    df["peak_id"] = df["candidate_name"].str.replace(r"__rep\d+$", "", regex=True)
    df["gc"] = df["sequence"].map(gc_fraction)

    pos_gc = df_pos_full[["peak_id", "sequence", "gc"]].copy()
    pos_gc.columns = ["peak_id", "pos_sequence", "pos_gc"]
    df = df.merge(pos_gc, on="peak_id", how="left")
    df["gc_abs_diff"] = (df["gc"] - df["pos_gc"]).abs()
    return df


def greedy_select_best_candidates(df_candidates: pd.DataFrame) -> pd.DataFrame:
    df = df_candidates.sort_values(["gc_abs_diff", "peak_id", "chrom", "start", "end"]).copy()
    selected_rows = []
    used_peak_ids = set()
    used_coords = set()
    used_sequences = set()

    positive_sequences = set(df["pos_sequence"].dropna().astype(str).tolist())

    for _, row in df.iterrows():
        peak_id = row["peak_id"]
        coord = (row["chrom"], int(row["start"]), int(row["end"]))
        seq = str(row["sequence"])
        if peak_id in used_peak_ids:
            continue
        if coord in used_coords:
            continue
        if seq in used_sequences:
            continue
        if seq in positive_sequences:
            continue
        selected_rows.append(row)
        used_peak_ids.add(peak_id)
        used_coords.add(coord)
        used_sequences.add(seq)

    if len(selected_rows) == 0:
        return pd.DataFrame(columns=df.columns)
    return pd.DataFrame(selected_rows)


def make_negative_label_bed(df_neg: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    out[0] = df_neg["chrom"].astype(str)
    out[1] = df_neg["start"].astype(int)
    out[2] = df_neg["end"].astype(int)
    out[3] = df_neg["peak_id"].astype(str) + "_NEG"
    out[4] = 0
    out[5] = "."
    out[6] = 0
    out[7] = "NA"
    out[8] = "NA"
    out[9] = "NA"
    return out


def split_negative_df(df_neg: pd.DataFrame):
    val = df_neg[df_neg["chrom"] == "chr4"].copy()
    test = df_neg[df_neg["chrom"].isin(["chr8", "chr9"])].copy()
    train = df_neg[~df_neg["chrom"].isin(["chr4", "chr8", "chr9"])].copy()
    return train, val, test


def write_neg_parts(train_neg_bed: str):
    df = pd.read_csv(train_neg_bed, sep="\t", header=None)
    df = df.sample(frac=1.0, random_state=0).reset_index(drop=True)
    split_size = len(df) // 4
    parts = [
        df.iloc[:split_size],
        df.iloc[split_size:2 * split_size],
        df.iloc[2 * split_size:3 * split_size],
        df.iloc[3 * split_size:],
    ]
    for i, part in enumerate(parts, start=1):
        part.to_csv(train_neg_bed[:-4] + ".Neg_Part{}.bed".format(i), sep="\t", header=None, index=False)


def cat_fastas(f1: str, f2: str, out_fa: str):
    ensure_dir(str(Path(out_fa).parent))
    with open(out_fa, "w") as fout:
        for fp in [f1, f2]:
            with open(fp, "r") as fin:
                shutil.copyfileobj(fin, fout)


def generate_genomic_gc_matched_negatives(df_pos_full: pd.DataFrame, cfg: Dict[str, str]) -> pd.DataFrame:
    work_dir = cfg["work_dir"]
    genome_fa = cfg["genome_fa"]
    name = cfg["name"]
    oversample_factor = int(cfg.get("oversample_factor", 30))
    seed = int(cfg.get("shuffle_seed", 19951124))

    exclusion_bed = build_exclusion_bed(
        all_possible_glob=cfg["all_possible_glob"],
        work_dir=work_dir,
        prefix=name,
        padding_bp=int(cfg.get("exclude_padding_bp", 0)),
        merge_distance=int(cfg.get("merge_distance", 200)),
    )

    template_bed = str(Path(work_dir) / (name + ".shuffle_template.bed"))
    shuffled_bed = str(Path(work_dir) / (name + ".shuffle_candidates.bed"))
    shuffled_tab = str(Path(work_dir) / (name + ".shuffle_candidates.tab"))

    build_oversampled_template_bed(df_pos_full, template_bed, oversample_factor)

    shuffle_cmd = (
        "bedtools shuffle "
        "-i {inp} "
        "-g {genome_sizes} "
        "-chrom "
        "-seed {seed} "
        "-excl {excl} "
        "> {out}"
    ).format(
        inp=template_bed,
        genome_sizes=cfg["genome_sizes"],
        seed=seed,
        excl=exclusion_bed,
        out=shuffled_bed,
    )
    run(shuffle_cmd)

    run("bedtools getfasta -fi {genome} -bed {bed} -name -tab > {out}".format(
        genome=genome_fa,
        bed=shuffled_bed,
        out=shuffled_tab,
    ))

    candidates = parse_shuffle_candidates(shuffled_bed, shuffled_tab, df_pos_full)
    best = greedy_select_best_candidates(candidates)

    n_pos = df_pos_full["peak_id"].nunique()
    n_neg = best["peak_id"].nunique()
    print("[INFO] matched negatives: {}/{}".format(n_neg, n_pos))
    if n_neg < n_pos:
        missing = sorted(list(set(df_pos_full["peak_id"]) - set(best["peak_id"])))
        miss_out = str(Path(work_dir) / (name + ".missing_negatives.txt"))
        with open(miss_out, "w") as fout:
            for x in missing:
                fout.write(str(x) + "\n")
        print("[WARN] {} positives did not get a unique genomic negative. Increase oversample_factor. Missing IDs saved to {}".format(n_pos - n_neg, miss_out))

    best = best.merge(
        df_pos_full[["peak_id", "chrom", "start", "end"]].rename(columns={
            "chrom": "pos_chrom",
            "start": "pos_start",
            "end": "pos_end",
        }),
        on="peak_id",
        how="left",
    )
    best["header"] = best["peak_id"] + "|NEG|" + best["chrom"] + ":" + best["start"].astype(str) + "-" + best["end"].astype(str)
    return best


# =========================================================
# main
# =========================================================
def validate_config(cfg: Dict[str, str]):
    require_file(cfg["narrowpeak_gz"], "CONFIG")
    require_file(cfg["genome_fa"], "CONFIG")
    require_file(cfg["genome_sizes"], "CONFIG")
    require_dir(cfg["neg_out_dir"], "CONFIG")

    ensure_dir(cfg["work_dir"])
    ensure_dir(cfg["bed_out_dir"])
    ensure_dir(cfg["fasta_out_dir"])
    ensure_dir(cfg["posneg_out_dir"])

    which_or_raise("bedtools")


def main(cfg: Dict[str, str]):
    validate_config(cfg)

    narrowpeak_txt = maybe_decompress(cfg["narrowpeak_gz"], cfg["work_dir"])
    df_pos = load_original_narrowpeak(narrowpeak_txt)

    paths, df_pos_full, train_pos, val_pos, test_pos = write_positive_outputs(df_pos, cfg)
    df_neg = generate_genomic_gc_matched_negatives(df_pos_full, cfg)

    name = cfg["name"]
    neg_out_dir = cfg["neg_out_dir"]
    bed_out_dir = cfg["bed_out_dir"]
    posneg_out_dir = cfg["posneg_out_dir"]

    full_neg_fa = str(Path(neg_out_dir) / (name + ".genomicGCmatch.original_length.fasta"))
    write_fasta_from_df(df_neg, full_neg_fa)

    train_neg, val_neg, test_neg = split_negative_df(df_neg)

    train_neg_fa = str(Path(neg_out_dir) / (name + ".genomicGCmatch.original_length.train.fasta"))
    val_neg_fa = str(Path(neg_out_dir) / (name + ".genomicGCmatch.original_length.validation.fasta"))
    test_neg_fa = str(Path(neg_out_dir) / (name + ".genomicGCmatch.original_length.test.fasta"))
    write_fasta_from_df(train_neg, train_neg_fa)
    write_fasta_from_df(val_neg, val_neg_fa)
    write_fasta_from_df(test_neg, test_neg_fa)

    train_neg_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.train.Neg.bed"))
    val_neg_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.validation.Neg.bed"))
    test_neg_bed = str(Path(bed_out_dir) / (name + ".narrowPeak.test.Neg.bed"))

    make_negative_label_bed(train_neg).to_csv(train_neg_bed, sep="\t", header=False, index=False)
    make_negative_label_bed(val_neg).to_csv(val_neg_bed, sep="\t", header=False, index=False)
    make_negative_label_bed(test_neg).to_csv(test_neg_bed, sep="\t", header=False, index=False)

    for pos_bed, neg_bed in [
        (paths["pos_train_bed"], train_neg_bed),
        (paths["pos_val_bed"], val_neg_bed),
        (paths["pos_test_bed"], test_neg_bed),
    ]:
        pos_df = pd.read_csv(pos_bed, sep="\t", header=None)
        neg_df = pd.read_csv(neg_bed, sep="\t", header=None)
        withneg_path = pos_bed[:-4] + ".withNeg.bed"
        pd.concat([pos_df, neg_df], axis=0).to_csv(withneg_path, sep="\t", header=False, index=False)

    write_neg_parts(train_neg_bed)

    train_posneg = str(Path(posneg_out_dir) / (name + ".narrowPeak.train.posneg.fasta"))
    val_posneg = str(Path(posneg_out_dir) / (name + ".narrowPeak.validation.posneg.fasta"))
    test_posneg = str(Path(posneg_out_dir) / (name + ".narrowPeak.test.posneg.fasta"))

    cat_fastas(paths["pos_train_fa"], train_neg_fa, train_posneg)
    cat_fastas(paths["pos_val_fa"], val_neg_fa, val_posneg)
    cat_fastas(paths["pos_test_fa"], test_neg_fa, test_posneg)

    summary_path = str(Path(cfg["work_dir"]) / (name + ".genomicGCmatch.summary.tsv"))
    df_neg.loc[:, ["peak_id", "pos_chrom", "pos_start", "pos_end", "chrom", "start", "end", "pos_gc", "gc", "gc_abs_diff"]].to_csv(
        summary_path, sep="\t", index=False
    )

    print("\n[DONE]")
    print("POS full fasta:", paths["pos_full_fa"])
    print("NEG full fasta:", full_neg_fa)
    print("Train/val/test posneg:", train_posneg, val_posneg, test_posneg)
    print("GC match summary:", summary_path)


# =========================================================
# EDIT CONFIG ONLY
# =========================================================
CONFIG = {
    "name": "fullard_DLPFC_neurons.optimal_peak",
    "narrowpeak_gz": "../peak_files/fullard_DLPFC/fullard_DLPFC_neurons.optimal_peak.narrowPeak.gz",
    "genome_fa": "../../genome/hg38.fa",
    "genome_sizes": "../../genome/hg38.chrom.sizes",

    "all_possible_glob": "../peak_files/fullard_DLPFC/*.narrowPeak.gz",
    "merge_distance": 200,
    "exclude_padding_bp": 0,

    "oversample_factor": 30,
    "shuffle_seed": 19951124,

    "work_dir": "./_work_fullard_DLPFC_original_length_genomic_neg",
    "bed_out_dir": "../original_length_peaks_hg38",
    "fasta_out_dir": "../original_length_peaks_fasta_hg38",
    "neg_out_dir": "/media/zihengc/T7/THP1_machinelearning/background_negatives_original_length_hg38",
    "posneg_out_dir": "../positive_negative_original_length_hg38",
}


if __name__ == "__main__":
    main(CONFIG)
