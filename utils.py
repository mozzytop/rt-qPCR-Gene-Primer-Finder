"""
Utility functions for the Gene Primer Lookup Tool.

Provides helpers for:
  - DNA character filtering
  - Reverse complement generation
  - Formatted sequence block rendering
  - Robust primer extraction with direction detection
"""

from __future__ import annotations

import re
from io import BytesIO
from typing import Optional
from dataclasses import dataclass
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


# Standard IUPAC complement map
_COMPLEMENT = str.maketrans("ATCGatcgUu", "TAGCtagcAa")


def filter_dna(sequence: str) -> str:
    """Strip everything except A, T, C, G, U (case-insensitive) from *sequence*."""
    return re.sub(r"[^ATCGUatcgu]", "", sequence)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA string (preserves case)."""
    return sequence.translate(_COMPLEMENT)[::-1]


def complement_only(sequence: str) -> str:
    """Return the complement of a DNA string WITHOUT reversing (preserves case)."""
    return sequence.translate(_COMPLEMENT)


def reverse_only(sequence: str) -> str:
    """Return the reversed sequence WITHOUT complementing (preserves case)."""
    return sequence[::-1]


def _normalize_dna(sequence: str) -> str:
    """Normalize DNA/RNA input to uppercase DNA characters."""
    return filter_dna(sequence).upper().replace("U", "T")


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
    primer_upper = _normalize_dna(primer)
    seq_upper = _normalize_dna(sequence)
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
    fwd_clean = _normalize_dna(forward)
    rev_clean = _normalize_dna(reverse)
    cds_upper = _normalize_dna(cds_sequence)
    rc_cds = reverse_complement(cds_upper)

    fwd_sense = find_primer_in_sequence(fwd_clean, cds_upper)
    fwd_anti = find_primer_in_sequence(fwd_clean, rc_cds)

    rev_sense = find_primer_in_sequence(rev_clean, cds_upper)
    rev_anti = find_primer_in_sequence(rev_clean, rc_cds)

    # Forward should be on sense, reverse on antisense (= its RC on sense)
    rev_rc = reverse_complement(rev_clean)
    rev_rc_on_sense = find_primer_in_sequence(rev_rc, cds_upper)
    forward_verified = fwd_sense is not None
    reverse_verified = rev_rc_on_sense is not None
    reverse_antisense_pos = find_primer_in_sequence(rev_clean, rc_cds)

    return {
        "forward_seq": fwd_clean,
        "reverse_seq": rev_clean,
        "reverse_complement": rev_rc,
        "fwd_sense_pos": fwd_sense,
        "fwd_antisense_pos": fwd_anti,
        "rev_sense_pos": rev_sense,
        "rev_antisense_pos": rev_anti,
        "rev_rc_sense_pos": rev_rc_on_sense,
        "reverse_antisense_binding_pos": reverse_antisense_pos,
        "fwd_maps": forward_verified,
        "rev_maps": reverse_verified,
        "both_map": forward_verified and reverse_verified,
    }


def build_report_pdf(report_text: str) -> bytes:
    """
    Build a simple monospaced PDF document from plain text using only stdlib.
    """
    lines = report_text.splitlines() or [""]
    lines_per_page = 52
    pages = [lines[i : i + lines_per_page] for i in range(0, len(lines), lines_per_page)] or [[""]]

    objects: list[bytes] = []
    page_ids: list[int] = []
    content_ids: list[int] = []

    def _add_object(data: str | bytes) -> int:
        payload = data.encode("latin-1", errors="replace") if isinstance(data, str) else data
        objects.append(payload)
        return len(objects)

    def _pdf_escape(text: str) -> str:
        return (
            text.replace("\\", "\\\\")
            .replace("(", "\\(")
            .replace(")", "\\)")
        )

    catalog_id = _add_object("")
    pages_id = _add_object("")
    font_id = _add_object("<< /Type /Font /Subtype /Type1 /BaseFont /Courier >>")

    for page_lines in pages:
        stream_lines = ["BT", "/F1 10 Tf", "72 770 Td", "12 TL"]
        for idx, line in enumerate(page_lines):
            prefix = "" if idx == 0 else "T* "
            stream_lines.append(f"{prefix}({_pdf_escape(line)}) Tj")
        stream_lines.append("ET")
        stream = "\n".join(stream_lines).encode("latin-1", errors="replace")
        content_id = _add_object(
            b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream"
        )
        page_id = _add_object(
            f"<< /Type /Page /Parent {pages_id} 0 R /MediaBox [0 0 612 792] "
            f"/Resources << /Font << /F1 {font_id} 0 R >> >> /Contents {content_id} 0 R >>"
        )
        content_ids.append(content_id)
        page_ids.append(page_id)

    objects[catalog_id - 1] = f"<< /Type /Catalog /Pages {pages_id} 0 R >>".encode("latin-1")
    kids = " ".join(f"{page_id} 0 R" for page_id in page_ids)
    objects[pages_id - 1] = f"<< /Type /Pages /Count {len(page_ids)} /Kids [{kids}] >>".encode("latin-1")

    output = BytesIO()
    output.write(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")

    offsets = [0]
    for index, obj in enumerate(objects, start=1):
        offsets.append(output.tell())
        output.write(f"{index} 0 obj\n".encode("ascii"))
        output.write(obj)
        output.write(b"\nendobj\n")

    xref_offset = output.tell()
    output.write(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    output.write(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        output.write(f"{offset:010d} 00000 n \n".encode("ascii"))

    output.write(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root {catalog_id} 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF"
        ).encode("ascii")
    )
    return output.getvalue()


def build_report_docx(report_text: str) -> bytes:
    """
    Build a minimal DOCX document from plain text using only stdlib.
    """
    body_parts = []
    for line in report_text.splitlines() or [""]:
        if line:
            run_parts = []
            for part in re.split(r"(\s+)", line):
                if not part:
                    continue
                preserve_attr = ' xml:space="preserve"' if part != part.strip() else ""
                run_parts.append(f"<w:r><w:t{preserve_attr}>{escape(part)}</w:t></w:r>")
            runs = "".join(run_parts)
        else:
            runs = "<w:r />"
        body_parts.append(f"<w:p>{runs}</w:p>")

    document_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas" '
        'xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" '
        'xmlns:o="urn:schemas-microsoft-com:office:office" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math" '
        'xmlns:v="urn:schemas-microsoft-com:vml" '
        'xmlns:wp14="http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing" '
        'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" '
        'xmlns:w10="urn:schemas-microsoft-com:office:word" '
        'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
        'xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml" '
        'xmlns:wpg="http://schemas.microsoft.com/office/word/2010/wordprocessingGroup" '
        'xmlns:wpi="http://schemas.microsoft.com/office/word/2010/wordprocessingInk" '
        'xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml" '
        'xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape" '
        'mc:Ignorable="w14 wp14">'
        '<w:body>'
        f'{"".join(body_parts)}'
        '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" '
        'w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/>'
        '<w:cols w:space="708"/><w:docGrid w:linePitch="360"/></w:sectPr>'
        '</w:body></w:document>'
    )

    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        '</Types>'
    )
    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="word/document.xml"/>'
        '</Relationships>'
    )
    doc_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships" />'
    )

    output = BytesIO()
    with ZipFile(output, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", rels_xml)
        zf.writestr("word/document.xml", document_xml)
        zf.writestr("word/_rels/document.xml.rels", doc_rels_xml)
    return output.getvalue()
