"""Data loading and intactness filtering for HIV sequence analysis.

This module implements progressive intactness filtering based on CFEIntact defect categories.
"""

from typing import Literal
from functools import cache
import csv
from pathlib import Path


# CFEIntact defect codes categorized by what they depend on:
STRUCTURAL_DEFECTS = {
    "APOBECHypermutation",
    "NonHIV",
    "Scramble",
    "InternalInversion",
    "UnknownNucleotide",
    "MissingORF",
    "LongDeletion",
}

DISTANCE_DEFECTS = STRUCTURAL_DEFECTS | {
    "Deletion",
    "Insertion",
    "MutatedStartCodon",
    "MutatedStopCodon",
    "InternalStop",
}

INDEL_DEFECTS = STRUCTURAL_DEFECTS | {
    "Deletion",
    "Insertion",
}

OTHER_DEFECTS = {
    "PackagingSignalDeletion",
    "MajorSpliceDonorSiteMutated",
}

def get_errors(ret: dict[str, list[dict]], path: Path) -> None:
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            qseqid = row["qseqid"]
            if qseqid not in ret:
                ret[qseqid] = []
            ret[qseqid].append(row)


@cache
def errors(path: Path) -> dict[str, list[dict]]:
    ret: dict[str, list[dict]] = {}
    get_errors(ret=ret, path=path)
    return ret


def is_structurally_intact(qseqid: str, path: Path) -> bool:
    """Check only structural defects (for size analysis)."""
    errs = errors(path).get(qseqid)
    if errs is None:
        return True
    codes = [err["code"] for err in errs]
    return not any(code in STRUCTURAL_DEFECTS for code in codes)


def is_distance_intact(qseqid: str, path: Path) -> bool:
    """Check structural + size defects (for distance analysis)."""
    errs = errors(path).get(qseqid)
    if errs is None:
        return True
    codes = [err["code"] for err in errs]
    return not any(code in DISTANCE_DEFECTS for code in codes)


def is_indel_intact(qseqid: str, path: Path) -> bool:
    """Check all CFEIntact defects (for indel analysis)."""
    errs = errors(path).get(qseqid)
    if errs is None:
        return True
    codes = [err["code"] for err in errs]
    return not any(code in INDEL_DEFECTS for code in codes)


def get_joined_it(source: Literal["los-alamos/plasma", "cfeintact/plasma", "cfeintact/all"]):
    if source == "los-alamos/plasma":
        path = "output/individual-plasma/joined.csv"
        defects_path = None
    elif source == "cfeintact/plasma":
        path = "output/fullgenomes-plasma/regions.csv"
        defects_path = Path("output/fullgenomes-plasma/defects.csv")
    elif source == "cfeintact/all":
        path = "output/fullgenomes-all/regions.csv"
        defects_path = Path("output/fullgenomes-all/defects.csv")
    else:
        raise ValueError(f"Invalid choice for source: {source!r}.")

    with open(path, "r") as f:
        r = csv.DictReader(f)
        for row in r:
            row["distance"] = float(row["distance"])
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])

            key = row["qseqid"]

            if defects_path is None:
                aminos = row["aminoacids"]
                has_stop = "*" in aminos[10:-10]
                row["size_structural_intact"] = not has_stop
                row["distance_intact"] = not has_stop
                row["indel_intact"] = not has_stop
            else:
                row["size_structural_intact"] = is_structurally_intact(key, defects_path)
                row["distance_intact"] = is_distance_intact(key, defects_path)
                row["indel_intact"] = is_indel_intact(key, defects_path)

            yield row


@cache
def get_joined(source: Literal["los-alamos/plasma", "cfeintact/all", "cfeintact/plasma"]):
    """Get cached joined data for a source.

    Args:
        source: Data source identifier

    Returns:
        List of sequence dictionaries with intactness annotations
    """
    return list(get_joined_it(source))
