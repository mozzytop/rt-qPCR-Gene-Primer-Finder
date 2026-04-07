"""
Utility functions for the Gene Primer Lookup Tool.

Provides helpers for:
  - DNA character filtering
  - Reverse complement generation
  - Formatted sequence block rendering
  - Robust primer extraction with direction detection
"""

import re
from typing import Optional
from dataclasses import dataclass


# Standard IUPAC complement map
_COMPLEMENT = str.maketrans("ATCGatcgUu", "TAGCtagcAa")


def filter_dna(sequence: str) -> str:
    """Strip everything except A, T, C, G, U (case-insensitive) from *sequence*."""
    return re.sub(r"[^ATCGUatcgu]", "", sequence)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA string (preserves case)."""
    return sequence.translate(_COMPLEMENT)[::-1]


def format_origin_block(sequence: str, line_width: int = 60, block_size: int = 10) -> str:
    """
    Format a DNA sequence in NCBI GenBank ORIGIN style.

    Example output line:
        1 atgttcaact cgatgacccc accaccaatc agtagctatg gcgagccctg ctgtctccgg
    """
    seq = sequence.lower()
    lines: list[str] = []
    for i in range(0, len(seq), line_width):
        chunk = seq[i : i + line_width]
        blocks = [chunk[j : j + block_size] for j in range(0, len(chunk), block_size)]
        lines.append(f"{i + 1:>9} {' '.join(blocks)}")
    return "\n".join(lines)


def format_filtered_dna(sequence: str, line_width: int = 60) -> str:
    """Format filtered DNA as continuous lowercase lines of *line_width* bases."""
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
    primer_upper = primer.upper().replace("U", "T")
    seq_upper = sequence.upper().replace("U", "T")
    idx = seq_upper.find(primer_upper)
    return idx if idx != -1 else None


def find_primer_on_either_strand(primer: str, sequence: str) -> tuple[Optional[int], str]:
    """
    Try to locate *primer* on the sense strand first, then on the antisense
    (reverse complement of the sequence).

    Returns (position, strand) where strand is "sense", "antisense", or "".
    """
    idx = find_primer_in_sequence(primer, sequence)
    if idx is not None:
        return idx, "sense"
    rc_seq = reverse_complement(sequence)
    idx = find_primer_in_sequence(primer, rc_seq)
    if idx is not None:
        return idx, "antisense"
    return None, ""


# ── Robust primer extraction ─────────────────────────────────────────


@dataclass
class ExtractedPrimer:
    """A primer sequence extracted from text with optional direction label."""
    sequence: str           # uppercase DNA only
    direction: str = ""     # "forward", "reverse", or "" if unknown
    context: str = ""       # surrounding text snippet for debugging


def extract_primers_from_text(text: str) -> list[str]:
    """
    Basic extraction: pull out primer-like sequences in 5'-...-3' notation.
    Returns a list of uppercase DNA strings (primer bodies only).
    """
    pattern = r"5['\u2019\u0027\u2032]\s*-?\s*([ACGTUacgtu]{15,40})\s*-?\s*3['\u2019\u0027\u2032]"
    return [m.upper().replace("U", "T") for m in re.findall(pattern, text)]


def extract_primers_with_direction(text: str) -> list[ExtractedPrimer]:
    """
    Advanced extraction: scan full text for primer sequences with direction
    context clues.

    Searches for:
      1. Explicit 5'-SEQUENCE-3' notation
      2. Bare DNA strings of 18-30 bp near keywords like "forward", "reverse",
         "sense", "antisense", "primer"

    For each hit, inspects the surrounding ~200 characters for direction
    keywords to assign "forward" or "reverse".

    Returns a deduplicated list of ``ExtractedPrimer``.
    """
    results: list[ExtractedPrimer] = []
    seen: set[str] = set()

    # ── Pattern 1: Explicit 5'-...-3' notation ────────────────────────
    pattern_explicit = (
        r"5['\u2019\u0027\u2032]\s*-?\s*"
        r"([ACGTUacgtu]{15,40})"
        r"\s*-?\s*3['\u2019\u0027\u2032]"
    )
    for m in re.finditer(pattern_explicit, text):
        seq = m.group(1).upper().replace("U", "T")
        if seq not in seen:
            seen.add(seq)
            start = max(0, m.start() - 200)
            end = min(len(text), m.end() + 200)
            context = text[start:end]
            direction = _infer_direction(context, m.start() - start)
            results.append(ExtractedPrimer(sequence=seq, direction=direction, context=context))

    # ── Pattern 2: Bare DNA strings near primer keywords ──────────────
    # Look for runs of 18-30 DNA bases that appear near "primer", "forward",
    # "reverse", "sense", "antisense", etc.
    primer_keyword_pattern = re.compile(
        r"(?:primer|forward|reverse|sense|antisense|fwd|rev|F\s*:|R\s*:)",
        re.IGNORECASE,
    )
    bare_dna_pattern = re.compile(r"\b([ACGTUacgtu]{18,35})\b")

    for m in bare_dna_pattern.finditer(text):
        seq = m.group(1).upper().replace("U", "T")
        if seq in seen:
            continue
        # Check if any primer keyword appears within 300 chars
        start_window = max(0, m.start() - 300)
        end_window = min(len(text), m.end() + 300)
        window = text[start_window:end_window]
        if primer_keyword_pattern.search(window):
            seen.add(seq)
            context = window
            direction = _infer_direction(context, m.start() - start_window)
            results.append(ExtractedPrimer(sequence=seq, direction=direction, context=context))

    return results


def _infer_direction(context: str, seq_position: int) -> str:
    """
    Inspect text surrounding a primer match to guess Forward vs Reverse.

    Looks at the 250 characters *before* the sequence position for:
      - "forward", "fwd", "sense" (not "antisense"), "F:" → forward
      - "reverse", "rev", "antisense", "R:" → reverse

    The closest keyword wins if both are found.
    """
    before = context[:seq_position].lower()

    fwd_keywords = [
        (before.rfind("forward"), "forward"),
        (before.rfind(" fwd"), "forward"),
        (before.rfind("f:"), "forward"),
        (before.rfind("f :"), "forward"),
    ]
    rev_keywords = [
        (before.rfind("reverse"), "reverse"),
        (before.rfind(" rev"), "reverse"),
        (before.rfind("antisense"), "reverse"),
        (before.rfind("r:"), "reverse"),
        (before.rfind("r :"), "reverse"),
    ]

    # Handle "sense" carefully — only count if NOT preceded by "anti"
    sense_pos = before.rfind("sense")
    if sense_pos != -1:
        prefix_start = max(0, sense_pos - 5)
        if "anti" in before[prefix_start:sense_pos]:
            rev_keywords.append((sense_pos, "reverse"))
        else:
            fwd_keywords.append((sense_pos, "forward"))

    best_fwd = max((pos for pos, _ in fwd_keywords if pos != -1), default=-1)
    best_rev = max((pos for pos, _ in rev_keywords if pos != -1), default=-1)

    if best_fwd == -1 and best_rev == -1:
        return ""
    if best_fwd > best_rev:
        return "forward"
    if best_rev > best_fwd:
        return "reverse"
    return ""


def verify_primer_pair(
    forward: str,
    reverse: str,
    cds_sequence: str,
) -> dict:
    """
    Verify that a forward and reverse primer pair maps to the CDS.

    The forward primer should match the sense strand.
    The reverse primer should match the antisense strand (i.e. its reverse
    complement is found in the sense strand).

    Returns a dict with verification results.
    """
    fwd_clean = forward.upper().replace("U", "T")
    rev_clean = reverse.upper().replace("U", "T")
    cds_upper = cds_sequence.upper().replace("U", "T")
    rc_cds = reverse_complement(cds_upper)

    fwd_sense = find_primer_in_sequence(fwd_clean, cds_upper)
    fwd_anti = find_primer_in_sequence(fwd_clean, rc_cds)

    rev_sense = find_primer_in_sequence(rev_clean, cds_upper)
    rev_anti = find_primer_in_sequence(rev_clean, rc_cds)

    # Forward should be on sense, reverse on antisense (= its RC on sense)
    rev_rc = reverse_complement(rev_clean)
    rev_rc_on_sense = find_primer_in_sequence(rev_rc, cds_upper)

    return {
        "forward_seq": fwd_clean,
        "reverse_seq": rev_clean,
        "reverse_complement": rev_rc,
        "fwd_sense_pos": fwd_sense,
        "fwd_antisense_pos": fwd_anti,
        "rev_sense_pos": rev_sense,
        "rev_antisense_pos": rev_anti,
        "rev_rc_sense_pos": rev_rc_on_sense,
        "fwd_maps": fwd_sense is not None or fwd_anti is not None,
        "rev_maps": rev_sense is not None or rev_anti is not None,
        "both_map": (fwd_sense is not None or fwd_anti is not None) and (rev_sense is not None or rev_anti is not None),
    }
