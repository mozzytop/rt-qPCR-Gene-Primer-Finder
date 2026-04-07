"""
Utility functions for the Gene Primer Lookup Tool.

Provides helpers for:
  - DNA character filtering
  - Reverse complement generation
  - Formatted sequence block rendering
"""

import re
from typing import Optional


# Standard IUPAC complement map
_COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def filter_dna(sequence: str) -> str:
    """Strip everything except A, T, C, G (case-insensitive) from *sequence*."""
    return re.sub(r"[^ATCGatcg]", "", sequence)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA string (preserves case)."""
    return sequence.translate(_COMPLEMENT)[::-1]


def format_origin_block(sequence: str, line_width: int = 60, block_size: int = 10) -> str:
    """
    Format a DNA sequence in NCBI GenBank ORIGIN style.

    Example output line:
        1 atgttcaact cgatgacccc accaccaatc agtagctatg gcgagccctg ctgtctccgg

    Parameters
    ----------
    sequence : str
        Raw (already-filtered) lowercase DNA string.
    line_width : int
        Number of bases per display line (default 60).
    block_size : int
        Number of bases per space-separated block (default 10).

    Returns
    -------
    str
        The fully formatted ORIGIN block WITHOUT the leading "ORIGIN" header
        or trailing "//".
    """
    seq = sequence.lower()
    lines: list[str] = []
    for i in range(0, len(seq), line_width):
        chunk = seq[i : i + line_width]
        blocks = [chunk[j : j + block_size] for j in range(0, len(chunk), block_size)]
        # Right-justify the position number to 9 characters (GenBank style)
        lines.append(f"{i + 1:>9} {' '.join(blocks)}")
    return "\n".join(lines)


def format_filtered_dna(sequence: str, line_width: int = 60) -> str:
    """
    Format filtered DNA as continuous lowercase lines of *line_width* bases.

    No numbering, no spaces — just raw FASTA-body style.
    """
    seq = sequence.lower()
    return "\n".join(seq[i : i + line_width] for i in range(0, len(seq), line_width))


def format_translation(translation: str, line_width: int = 60) -> str:
    """Wrap a protein translation string to *line_width* characters per line."""
    return "\n".join(translation[i : i + line_width] for i in range(0, len(translation), line_width))


def find_primer_in_sequence(primer: str, sequence: str) -> Optional[int]:
    """
    Search for *primer* inside *sequence* (case-insensitive).

    Returns the 0-based start position or ``None`` if not found.
    """
    primer_upper = primer.upper()
    seq_upper = sequence.upper()
    idx = seq_upper.find(primer_upper)
    return idx if idx != -1 else None


def find_primer_on_either_strand(primer: str, sequence: str) -> tuple[Optional[int], str]:
    """
    Try to locate *primer* on the sense strand first, then on the antisense
    (reverse complement of the sequence).

    Returns
    -------
    (position, strand)
        *position* is 0-based index or None;  *strand* is ``"sense"`` or
        ``"antisense"`` or ``""`` if not found.
    """
    idx = find_primer_in_sequence(primer, sequence)
    if idx is not None:
        return idx, "sense"
    # Check antisense: reverse-complement of the CDS, then search
    rc_seq = reverse_complement(sequence)
    idx = find_primer_in_sequence(primer, rc_seq)
    if idx is not None:
        return idx, "antisense"
    return None, ""


def extract_primers_from_text(text: str) -> list[str]:
    """
    Use regex to pull out primer-like sequences enclosed in 5'-...-3' notation
    from free text (abstract / full-text snippet).

    Returns a list of uppercase DNA strings (the primer bodies only).
    """
    # Pattern: 5' followed by optional dash/space, then DNA bases, then
    # optional dash/space and 3'
    pattern = r"5['\u2019]\s*-?\s*([ACGTacgt]{15,40})\s*-?\s*3['\u2019]"
    return [m.upper() for m in re.findall(pattern, text)]
