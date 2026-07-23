from pathlib import Path
import sys


def find_repo_root(start: Path) -> Path:
    """
    Find repository root by searching upward for src/mpra.
    This makes the script portable across macOS, Linux, and Windows.
    """
    start = start.resolve()

    for path in [start, *start.parents]:
        if (path / "src" / "mpra").is_dir():
            return path

    raise FileNotFoundError(
        "Could not find repository root. Expected to find 'src/mpra' "
        "in the current directory or one of its parent directories."
    )


def main():
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    src_dir = repo_root / "src"
    sys.path.insert(0, str(src_dir))

    from mpra.map_reads import MPRA

    library = repo_root / "indexing" / "MPRA3_Contributor_20231108_unique_GeneName_BarcodeEnhancerPair.csv"
    fastq_dir = repo_root / "mpra_fastq"

    if not library.exists():
        raise FileNotFoundError(f"Library file not found: {library}")

    if not fastq_dir.exists():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")

    mpra = MPRA(
        test_library=str(library),
        flank_bp=5,
    )

    print("Reference read, start building count tables")
    print(f"Repo root: {repo_root}")
    print(f"Library: {library}")
    print(f"FASTQ directory: {fastq_dir}")
    print(f"Using flank_bp = {mpra.flank_bp}")
    print(f"Using mapping column = {mpra.match_col}")

    mpra.loop_through_samples_and_build_count_tables(
        dir_folder=str(fastq_dir)
    )


if __name__ == "__main__":
    main()
