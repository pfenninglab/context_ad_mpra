#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

import pandas as pd


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for path in [start, *start.parents]:
        if (path / "src" / "mpra").is_dir():
            return path
    raise FileNotFoundError("Could not find repository root. Expected src/mpra.")


def paired_fastq_files(fastq_dir: Path) -> dict[str, dict[str, Path]]:
    r1_suffix = "_R1_001.fastq.gz"
    r2_suffix = "_R2_001.fastq.gz"
    paired: dict[str, dict[str, Path]] = {}

    for path in sorted(fastq_dir.glob("*.fastq.gz")):
        name = path.name
        if name.startswith(".") or name.startswith("._"):
            continue
        if name.endswith(r1_suffix):
            sample = name[: -len(r1_suffix)]
            paired.setdefault(sample, {})["R1"] = path
        elif name.endswith(r2_suffix):
            sample = name[: -len(r2_suffix)]
            paired.setdefault(sample, {})["R2"] = path

    return paired


def fastq_sequences(path: Path):
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            yield seq


def find_barcode_and_r1_flank(read1: str, flank_bp: int) -> tuple[str | None, str | None]:
    for offset in range(39, 42):
        if read1[offset : offset + 5] in {"TCTAA", "TCTGA"}:
            barcode = read1[offset + 5 : offset + 21]
            r1_flank = read1[offset + 33 : offset + 33 + flank_bp]
            return barcode, r1_flank
    return None, None


def build_r2_r1_key(read2: str, barcode: str, r1_flank: str, flank_bp: int) -> str | None:
    for offset in (0, -1, 1):
        if read2[46 + offset : 51 + offset] == "TGCTA":
            r2_flank = read2[offset + 51 : offset + 51 + flank_bp]
            return f"{barcode}_{r2_flank}_{r1_flank}"
    return None


def count_sample(
    sample: str,
    r1_path: Path,
    r2_path: Path,
    valid_keys: set[str],
    flank_bp: int,
    progress_every: int,
) -> tuple[Counter, int]:
    counts: Counter[str] = Counter()
    total = 0

    for read1, read2 in zip(fastq_sequences(r1_path), fastq_sequences(r2_path)):
        total += 1
        barcode, r1_flank = find_barcode_and_r1_flank(read1, flank_bp)
        if barcode:
            key = build_r2_r1_key(read2, barcode, r1_flank, flank_bp)
            if key in valid_keys:
                counts[key] += 1

        if progress_every and total % progress_every == 0:
            print(
                f"[PROGRESS] {sample}: processed {total / 1_000_000:.1f}M reads, "
                f"mapped {sum(counts.values()) / 1_000_000:.1f}M",
                flush=True,
            )

    return counts, total


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream Neuron FASTQs into a 10bp R2+R1 count table."
    )
    parser.add_argument("--base-dir", default=None)
    parser.add_argument("--fastq-dir", default="mpra_fastq")
    parser.add_argument(
        "--library",
        default="indexing/MPRA3_Contributor_20231108_unique_GeneName_BarcodeEnhancerPair.csv",
    )
    parser.add_argument("--out", default="outputs/read_counts_R1R2/Neuron_DNA_RNA.csv")
    parser.add_argument("--flank-bp", default=10, type=int)
    parser.add_argument("--progress-every", default=5_000_000, type=int)
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve() if args.base_dir else find_repo_root(Path(__file__).parent)
    fastq_dir = (base_dir / args.fastq_dir).resolve()
    library_path = (base_dir / args.library).resolve()
    out_path = (base_dir / args.out).resolve()

    match_col = f"Barcode_RC_Enhancer{args.flank_bp}_R2_R1"
    df = pd.read_csv(library_path)
    if match_col not in df.columns:
        raise KeyError(f"Library is missing mapping column: {match_col}")

    valid_keys = set(df[match_col].dropna().astype(str))
    paired = paired_fastq_files(fastq_dir)
    if not paired:
        raise FileNotFoundError(f"No paired FASTQ files found in {fastq_dir}")

    print(f"[INFO] base_dir: {base_dir}", flush=True)
    print(f"[INFO] fastq_dir: {fastq_dir}", flush=True)
    print(f"[INFO] library: {library_path}", flush=True)
    print(f"[INFO] output: {out_path}", flush=True)
    print(f"[INFO] mapping column: {match_col}", flush=True)
    print(f"[INFO] samples: {sorted(paired)}", flush=True)

    for sample in sorted(paired):
        files = paired[sample]
        if "R1" not in files or "R2" not in files:
            print(f"[SKIP] {sample}: missing R1 or R2", flush=True)
            continue

        print(f"[START] {sample}", flush=True)
        counts, total = count_sample(
            sample=sample,
            r1_path=files["R1"],
            r2_path=files["R2"],
            valid_keys=valid_keys,
            flank_bp=args.flank_bp,
            progress_every=args.progress_every,
        )
        df[sample] = df[match_col].map(counts)
        mapped = int(sum(counts.values()))
        frac = round(mapped / total, 4) if total else 0
        print(
            f"[OK] {sample}: mapped={mapped} total={total} fraction={frac}",
            flush=True,
        )

        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=False)

    # Put count columns first, mirroring the previous Neuron_DNA_RNA.csv layout.
    count_cols = [sample for sample in sorted(paired) if sample in df.columns]
    metadata_cols = [col for col in df.columns if col not in count_cols]
    df = df[count_cols + metadata_cols]
    df.to_csv(out_path, index=False)
    print(f"[DONE] wrote {out_path} shape={df.shape}", flush=True)


if __name__ == "__main__":
    main()
