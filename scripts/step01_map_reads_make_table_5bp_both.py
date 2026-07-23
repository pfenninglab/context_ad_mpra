from pathlib import Path
import sys


def find_repo_root(start: Path) -> Path:
    start = start.resolve()

    for path in [start, *start.parents]:
        if (path / "src" / "mpra").is_dir():
            return path

    raise FileNotFoundError("Could not find repository root. Expected to find 'src/mpra'.")


def main():
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    sys.path.insert(0, str(repo_root / "src"))

    from mpra.map_reads import MPRA

    library = repo_root / "indexing" / "MPRA3_Contributor_20231108_unique_GeneName_BarcodeEnhancerPair.csv"
    fastq_dir = repo_root / "mpra_fastq"

    mpra = MPRA(
        test_library=str(library),
        flank_bp=5,
        match_mode="both",
    )

    print("Reference read, start building count tables")
    print(f"Using flank_bp = {mpra.flank_bp}")
    print(f"Using match_mode = {mpra.match_mode}")
    print(f"Using mapping column = {mpra.match_col}")

    mpra.loop_through_samples_and_build_count_tables(dir_folder=str(fastq_dir))


if __name__ == "__main__":
    main()
