"""GenBank cleaner: automate the mundane annotation cleanup, never touch the sequence.

Ground rules:
- the sequence and topology must be identical before and after, or the file fails loudly;
- every dropped or merged feature leaves a trace (a /note on the survivor or a COMMENT line);
- the original file is never overwritten: callers write a cleaned copy plus a change report;
- running the cleaner twice gives the same result (idempotent).

Rules (ids are stable; DEFAULT_ON is the safe set):
  exact_duplicates  same key, location and all qualifiers -> keep one
  merge_same        same key + strand + span + label (qualifiers differ) -> one feature, qualifiers united
  generic_labels    same key + strand + span, labels generic ("Insert", "Insert Sequence", unlabelled
                    misc_feature...) -> keep one, other labels noted
  unlabelled_dupe   unlabelled feature with exactly the span of a labelled one -> drop
  unescape_labels   'LjUBI1\\promoter\\(no\\BsaI)' -> 'LjUBI1 promoter (no BsaI)'
  map_keys          non-INSDC keys (Vector, Site, Exon, Misc._feature, Position, seq, Tag) -> INSDC key,
                    with /note="original key: X"; resistance-gene keys become CDS only if a clean ORF
  drop_vntifkey     Vector NTI internal qualifier
  synonyms          normalise common labels (table below); the old label goes into /note
  map_unknown_keys  any key still not an INSDC (or legacy promoter/terminator...) key: case-folded
                    ('Gene' -> gene), a few with a real target ('mutation' -> variation, 'Ori' ->
                    rep_origin), everything else -> misc_feature; /note="original key: X", and an
                    unlabelled feature takes X as its label. Importers drop or rename keys they
                    do not know (CLC, Vector NTI, ApE and hand-made types: 'Module', 'Seq-region')
  strip_quotes      CLC exports a second, quoted copy of a label (/label=""35S promoter"") ->
                    quotes stripped, the duplicate dropped
  english_dates     LOCUS dates in Danish ('14-okt-2008') -> '14-OCT-2008'; Biopython would
                    otherwise write 01-JAN-1980
Opt-in (lose structure or need judgement):
  drop_colours      ApE colour qualifiers
  fragments_to_comment   'Fragment (...)' assembly-history features -> COMMENT, then dropped
  vector_to_comment      whole-plasmid 'Vector' features -> COMMENT, then dropped
  clc_primer_hits   CLC auto primer annotations '<name> (mismatches: N)', N>0, whose name is not a
                    lab oligo name/ID -> dropped
  name_unlabelled   give unlabelled features an obvious placeholder label
Always reported, never fixed: features past the sequence end, CDS problems, origin-wrapping features
on a linear record, "no BsaI/BpiI" label claims contradicted by the sequence.
"""

from __future__ import annotations

import io
import re
import warnings
from collections import defaultdict
from dataclasses import dataclass, field

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from ..formats import decode
from ..sequtil import revcomp

DEFAULT_ON = [
    "exact_duplicates",
    "merge_same",
    "generic_labels",
    "unlabelled_dupe",
    "unescape_labels",
    "map_keys",
    "drop_vntifkey",
    "synonyms",
    "map_unknown_keys",
    "strip_quotes",
    "english_dates",
]
OPT_IN = [
    "drop_colours",
    "fragments_to_comment",
    "vector_to_comment",
    "clc_primer_hits",
    "name_unlabelled",
]
ALL_RULES = DEFAULT_ON + OPT_IN

LABEL_KEYS = ("label", "gene", "product", "standard_name", "note")
GENERIC_LABELS = {
    "insert",
    "insert sequence",
    "region",
    "feature",
    "new feature",
    "misc_feature",
    "misc feature",
    "fragment",
    "sequence",
    "selection",
    "unnamed",
}
KEY_MAP = {
    "Misc._feature": "misc_feature",
    "Exon": "exon",
    "Vector": "misc_feature",
    "Site": "misc_feature",
    "Position": "misc_feature",
    "seq": "misc_feature",
    "Tag": "misc_feature",
    "MSC": "misc_feature",
    "Intron": "intron",
    "Promoter": "promoter",
    "Terminator": "terminator",
    "Primer": "primer_bind",
    "primer": "primer_bind",
}
INSDC_KEYS = set(
    """assembly_gap C_region CDS centromere D-loop D_segment exon gap gene iDNA intron J_segment
    mat_peptide misc_binding misc_difference misc_feature misc_recomb misc_RNA misc_structure
    mobile_element modified_base mRNA ncRNA N_region old_sequence operon oriT polyA_site
    precursor_RNA prim_transcript primer_bind propeptide protein_bind regulatory repeat_region
    rep_origin rRNA S_region sig_peptide source stem_loop STS telomere tmRNA transit_peptide tRNA
    unsure V_region V_segment variation 3'UTR 5'UTR""".split()
)
# retired from INSDC in 2014, still read by SnapGene, Benchling and Geneious: kept as they are
LEGACY_KEYS = set(
    """promoter terminator enhancer polyA_signal RBS -10_signal -35_signal misc_signal attenuator
    CAAT_signal TATA_signal GC_signal LTR repeat_unit satellite conflict""".split()
)
KNOWN_KEYS = INSDC_KEYS | LEGACY_KEYS
_KNOWN_FOLDED = {k.lower(): k for k in KNOWN_KEYS}
# lower-cased non-standard key -> INSDC key; anything else unknown becomes misc_feature
UNKNOWN_KEY_MAP = {
    "mutation": "variation",
    "snp": "variation",
    "difference": "misc_difference",
    "ori": "rep_origin",
    "origin": "rep_origin",
    "replication_ori": "rep_origin",
    "replication_origin": "rep_origin",
    "primer": "primer_bind",
    "primer_binding": "primer_bind",
    "primer_binding_": "primer_bind",
    "primer_binding_site": "primer_bind",
    "exons": "exon",
    "introns": "intron",
    "polya": "polyA_signal",
    "polya_signal": "polyA_signal",
}
MONTHS = "JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split()
# Danish abbreviations as CLC writes them on a Danish-locale machine, plus full English names
MONTH_ALIASES = {"MAJ": "MAY", "OKT": "OCT"} | {
    m: m[:3]
    for m in "JANUARY FEBRUARY MARCH APRIL JUNE JULY AUGUST SEPTEMBER OCTOBER NOVEMBER DECEMBER".split()
}
DATE = re.compile(r"^(\d{1,2})-([A-Za-z]+)\.?-(\d{4})$")
RESISTANCE_KEYS = {
    "AmpR",
    "KanR",
    "spec/strepR",
    "Spec/StrepR",
    "CmR",
    "TetR",
    "HygR",
    "BarR",
}
# Editable synonym table: lower-cased variant -> canonical label. Labels only; keys are never changed.
SYNONYMS = {
    "amp r": "AmpR",
    "ampr": "AmpR",
    "amp": "AmpR",
    "bla": "AmpR",
    "ampicillin": "AmpR",
    "amp resistance": "AmpR",
    "kan": "KanR",
    "kanr": "KanR",
    "kan r": "KanR",
    "kanamycin": "KanR",
    "spec/strepr": "SpecR/StrepR",
    "spec/strep": "SpecR/StrepR",
    "specr": "SpecR/StrepR",
    "aada": "SpecR/StrepR",
    "35s ter": "35S terminator",
    "t35s": "35S terminator",
    "35s term": "35S terminator",
    "35s terminator": "35S terminator",
    "ori": "ori",
    "origin": "ori",
    "p(lac)": "lac promoter",
    "p lac prom": "lac promoter",
    "lac promoter": "lac promoter",
    "plac": "lac promoter",
}
COLOUR_QUALS = ("ApEinfo_fwdcolor", "ApEinfo_revcolor", "ApEinfo_graphicformat")
CLC_PRIMER = re.compile(r"^(?P<name>.+?) \(mismatches: (?P<n>\d+)\)$")
FRAGMENT = re.compile(r"^Fragment \((?P<what>.*)\)$")


@dataclass
class Change:
    rule: str
    action: (
        str  # dropped | merged | relabelled | rekeyed | edited | commented | flagged
    )
    feature: str
    detail: str = ""


@dataclass
class CleanResult:
    record: object  # Bio.SeqRecord
    changes: list[Change] = field(default_factory=list)
    flags: list[Change] = field(default_factory=list)
    n_before: int = 0
    n_after: int = 0

    def genbank(self) -> str:
        # some exports carry no molecule type, and Biopython refuses to write without one
        self.record.annotations.setdefault("molecule_type", "DNA")
        out = io.StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            SeqIO.write(self.record, out, "genbank")
        return out.getvalue()


def label_of(f: SeqFeature) -> str | None:
    for k in LABEL_KEYS:
        vals = [
            v
            for v in f.qualifiers.get(k, [])
            if not str(v).startswith(("original key:", "cleaner:"))
        ]
        if vals:
            return str(vals[0])
    return None


def _desc(f: SeqFeature) -> str:
    return f"{f.type} '{label_of(f) or '(unlabelled)'}' {f.location}"


def _span(f: SeqFeature) -> tuple:
    return (str(f.location), f.location.strand or 0)


def _nspan(f: SeqFeature) -> tuple:
    """Numeric span: ignores partial-end markers (<, >), which tools add and drop inconsistently."""
    return tuple((int(p.start), int(p.end)) for p in f.location.parts) + (
        f.location.strand or 0,
    )


def _set_label(f: SeqFeature, new: str) -> None:
    key = next(
        (
            k
            for k in ("label", "gene", "product", "standard_name")
            if f.qualifiers.get(k)
        ),
        "label",
    )
    f.qualifiers[key] = [new]


def _add_note(f: SeqFeature, note: str) -> None:
    notes = f.qualifiers.setdefault("note", [])
    if note not in notes:
        notes.append(note)


def _is_generic(f: SeqFeature) -> bool:
    lab = label_of(f)
    return (
        lab is None
        or lab.strip().lower() in GENERIC_LABELS
        or lab.strip().lower() == f.type.lower()
    )


def _clean_orf(seq: str) -> bool:
    if len(seq) % 3 or len(seq) < 90:
        return False
    prot = str(Seq(seq).translate())
    return prot[0] == "M" and prot.endswith("*") and "*" not in prot[:-1]


def _feature_seq(f: SeqFeature, rec) -> str:
    return str(f.location.extract(rec.seq)).upper()


def _comment(rec, line: str) -> None:
    c = rec.annotations.get("comment", "")
    if line not in c:
        rec.annotations["comment"] = (c + "\n" if c else "") + line


def clean_record(
    rec, rules: list[str] | None = None, lab_oligo_names: set[str] | None = None
) -> CleanResult:
    rules = DEFAULT_ON if rules is None else rules
    unknown = set(rules) - set(ALL_RULES)
    if unknown:
        raise ValueError(f"unknown rule(s): {', '.join(sorted(unknown))}")
    seq_before, topo_before = str(rec.seq), rec.annotations.get("topology")
    res = CleanResult(rec, n_before=len(rec.features))
    feats = list(rec.features)
    L = len(rec.seq)

    def drop(f, rule, detail=""):
        feats.remove(f)
        res.changes.append(Change(rule, "dropped", _desc(f), detail))

    # a location Biopython cannot read (CLC writes 'Site -1^1') cannot be compared or written back
    for f in [f for f in feats if f.location is None]:
        drop(f, "check", "location unreadable, cannot be written to GenBank")

    # 1. label hygiene first, so later comparisons see clean labels
    for f in feats:
        if "strip_quotes" in rules:
            for k in LABEL_KEYS:
                vals = f.qualifiers.get(k)
                if not vals or not any(
                    isinstance(v, str) and len(v) > 1 and v[0] == v[-1] == '"'
                    for v in vals
                ):
                    continue
                new = list(
                    dict.fromkeys(
                        v[1:-1].strip()
                        if isinstance(v, str) and len(v) > 1 and v[0] == v[-1] == '"'
                        else v
                        for v in vals
                    )
                )
                f.qualifiers[k] = new
                res.changes.append(
                    Change("strip_quotes", "edited", _desc(f), f"/{k} {vals} -> {new}")
                )
        lab = label_of(f)
        if "unescape_labels" in rules and lab and "\\" in lab:
            new = re.sub(r"\s+", " ", lab.replace("\\", " ")).strip()
            _set_label(f, new)
            res.changes.append(
                Change("unescape_labels", "relabelled", _desc(f), f"'{lab}' -> '{new}'")
            )
        if "map_keys" in rules and f.type in RESISTANCE_KEYS | set(KEY_MAP):
            old = f.type
            if old in RESISTANCE_KEYS:
                new_key = "CDS" if _clean_orf(_feature_seq(f, rec)) else "misc_feature"
            else:
                new_key = KEY_MAP[old]
            if not label_of(f) and old.lower() != new_key.lower():
                _set_label(
                    f, old
                )  # an unlabelled 'Vector'/'Site' keeps what it was as its label
            vnti = f.qualifiers.get("vntifkey")
            f.type = new_key
            _add_note(f, f"original key: {old}")
            res.changes.append(
                Change(
                    "map_keys",
                    "rekeyed",
                    _desc(f),
                    f"{old} -> {new_key}" + (f" (vntifkey {vnti[0]})" if vnti else ""),
                )
            )
        if "map_unknown_keys" in rules and f.type not in KNOWN_KEYS:
            old = f.type
            new_key = _KNOWN_FOLDED.get(old.lower()) or UNKNOWN_KEY_MAP.get(
                old.lower(), "misc_feature"
            )
            if not label_of(f):
                _set_label(f, old)
            f.type = new_key
            _add_note(f, f"original key: {old}")
            res.changes.append(
                Change("map_unknown_keys", "rekeyed", _desc(f), f"{old} -> {new_key}")
            )
        if "drop_vntifkey" in rules and "vntifkey" in f.qualifiers:
            del f.qualifiers["vntifkey"]
            res.changes.append(
                Change("drop_vntifkey", "edited", _desc(f), "removed /vntifkey")
            )
        if "drop_colours" in rules and any(q in f.qualifiers for q in COLOUR_QUALS):
            for q in COLOUR_QUALS:
                f.qualifiers.pop(q, None)
            res.changes.append(
                Change("drop_colours", "edited", _desc(f), "removed ApE colours")
            )
        lab = label_of(f)
        if (
            "synonyms" in rules
            and lab
            and lab.strip().lower() in SYNONYMS
            and SYNONYMS[lab.strip().lower()] != lab
        ):
            new = SYNONYMS[lab.strip().lower()]
            _set_label(f, new)
            _add_note(f, f"cleaner: label was '{lab}'")
            res.changes.append(
                Change("synonyms", "relabelled", _desc(f), f"'{lab}' -> '{new}'")
            )

    # 2. provenance features become COMMENT lines (opt-in)
    for f in list(feats):
        lab = label_of(f) or ""
        if "fragments_to_comment" in rules and (m := FRAGMENT.match(lab)):
            _comment(rec, f"assembled from: {m.group('what')} [{f.location}]")
            drop(f, "fragments_to_comment", "moved to COMMENT")
        elif (
            "vector_to_comment" in rules
            and (
                f.type in ("Vector",)
                or "original key: Vector" in f.qualifiers.get("note", [])
            )
            and len(f.location) >= 0.5 * L
        ):
            _comment(rec, f"backbone: {lab or 'vector'} [{f.location}]")
            drop(f, "vector_to_comment", "moved to COMMENT")
        elif (
            "clc_primer_hits" in rules
            and (m := CLC_PRIMER.match(lab))
            and int(m.group("n")) > 0
        ):
            name = m.group("name")
            if lab_oligo_names is None or name not in lab_oligo_names:
                drop(
                    f,
                    "clc_primer_hits",
                    f"CLC primer hit with {m.group('n')} mismatches",
                )

    # 3. duplicates and stacked features
    if "exact_duplicates" in rules:
        seen = {}
        for f in list(feats):
            key = (
                f.type,
                _span(f),
                tuple(sorted((k, tuple(v)) for k, v in f.qualifiers.items())),
            )
            if key in seen:
                drop(f, "exact_duplicates", "identical to another feature")
            else:
                seen[key] = f
    groups: dict[tuple, list[SeqFeature]] = defaultdict(list)
    for f in feats:
        groups[(f.type, _nspan(f))].append(f)
    for (ftype, _), group in groups.items():
        if len(group) < 2:
            continue
        keep = group[0]
        for f in group[1:]:
            same_label = (label_of(f) or "").strip().lower() == (
                label_of(keep) or ""
            ).strip().lower()
            if "merge_same" in rules and same_label:
                for k, v in f.qualifiers.items():
                    if k not in keep.qualifiers:
                        keep.qualifiers[k] = v
                    elif (
                        k not in ("translation", "codon_start", "product")
                        and v != keep.qualifiers[k]
                    ):
                        keep.qualifiers[k] = list(dict.fromkeys(keep.qualifiers[k] + v))
                _add_note(keep, "cleaner: merged a duplicate annotation")
                drop(f, "merge_same", f"merged into {_desc(keep)}")
            elif "generic_labels" in rules and _is_generic(f) and _is_generic(keep):
                if label_of(f) and label_of(f) != label_of(keep):
                    _add_note(keep, f"cleaner: also labelled '{label_of(f)}'")
                drop(
                    f,
                    "generic_labels",
                    f"stacked generic annotation, kept {_desc(keep)}",
                )
            elif "generic_labels" in rules and _is_generic(f) and not _is_generic(keep):
                drop(f, "generic_labels", f"generic annotation under {_desc(keep)}")
            elif "generic_labels" in rules and _is_generic(keep) and not _is_generic(f):
                feats.remove(keep)
                res.changes.append(
                    Change(
                        "generic_labels",
                        "dropped",
                        _desc(keep),
                        f"generic annotation under {_desc(f)}",
                    )
                )
                keep = f
    if (
        "generic_labels" in rules
    ):  # generic misc_feature exactly under a named feature of another key
        named = {}
        for f in feats:
            if not _is_generic(f):
                named.setdefault(_nspan(f), f)
        for f in list(feats):
            if (
                f.type == "misc_feature"
                and _is_generic(f)
                and label_of(f)
                and _nspan(f) in named
            ):
                other = named[_nspan(f)]
                _add_note(other, f"cleaner: also annotated as '{label_of(f)}'")
                drop(f, "generic_labels", f"generic annotation under {_desc(other)}")
    if "unlabelled_dupe" in rules:
        labelled_spans = {_span(f) for f in feats if label_of(f)}
        for f in list(feats):
            if not label_of(f) and _span(f) in labelled_spans:
                drop(f, "unlabelled_dupe", "a labelled feature has exactly this span")
    if "name_unlabelled" in rules:
        for f in feats:
            if not label_of(f) and f.type != "source":
                new = f"unlabelled_{f.type}_{int(f.location.start) + 1}-{int(f.location.end)}"
                _set_label(f, new)
                res.changes.append(
                    Change("name_unlabelled", "relabelled", _desc(f), new)
                )

    # 4. report-only checks
    generic = sorted(
        (f for f in feats if _is_generic(f) and f.type != "source"),
        key=lambda f: int(f.location.start),
    )
    for i, f in enumerate(generic):
        for g in generic[i + 1 : i + 6]:
            if g is f or _nspan(g) == _nspan(f):
                continue
            if (
                abs(int(g.location.start) - int(f.location.start)) <= 5
                and abs(int(g.location.end) - int(f.location.end)) <= 5
            ):
                res.flags.append(
                    Change(
                        "check",
                        "flagged",
                        _desc(g),
                        f"near-duplicate (±5 bp) of {_desc(f)}; not merged",
                    )
                )

    circular = rec.annotations.get("topology") == "circular"
    for f in feats:
        if f.type == "source":
            continue
        if int(f.location.end) > L:
            res.flags.append(
                Change(
                    "check", "flagged", _desc(f), f"runs past the sequence end ({L} bp)"
                )
            )
        if (
            not circular
            and len(f.location.parts) > 1
            and int(f.location.parts[0].start) > int(f.location.parts[-1].start)
            and (f.location.strand or 1) == 1
        ):
            res.flags.append(
                Change(
                    "check",
                    "flagged",
                    _desc(f),
                    "wraps the origin but the record is linear",
                )
            )
        if f.type == "CDS":
            s = _feature_seq(f, rec)
            cs = int(f.qualifiers.get("codon_start", ["1"])[0]) - 1
            body = s[cs:]
            if len(body) % 3:
                res.flags.append(
                    Change(
                        "check",
                        "flagged",
                        _desc(f),
                        f"CDS length {len(body)} is not a multiple of 3",
                    )
                )
            elif body:
                prot = str(Seq(body).translate())
                if "*" in prot[:-1]:
                    res.flags.append(
                        Change(
                            "check",
                            "flagged",
                            _desc(f),
                            f"internal stop codon at aa {prot.index('*') + 1}",
                        )
                    )
                tr = f.qualifiers.get("translation", [None])[0]
                if tr and tr.rstrip("*") != prot.rstrip("*"):
                    res.flags.append(
                        Change(
                            "check",
                            "flagged",
                            _desc(f),
                            "/translation disagrees with the sequence",
                        )
                    )
        lab = (label_of(f) or "").lower()
        for enz, site in (("BsaI", "GGTCTC"), ("BpiI", "GAAGAC")):
            if f"no {enz.lower()}" in lab or (
                enz == "BpiI" and "or bpii" in lab and "no " in lab
            ):
                s = _feature_seq(f, rec)
                if site in s or revcomp(site) in s:
                    res.flags.append(
                        Change(
                            "check",
                            "flagged",
                            _desc(f),
                            f"label says no {enz} site, but the feature contains one",
                        )
                    )
    for f in feats:
        if not label_of(f) and f.type != "source" and "name_unlabelled" not in rules:
            res.flags.append(Change("check", "flagged", _desc(f), "unlabelled"))

    if "english_dates" in rules and isinstance(rec.annotations.get("date"), str):
        old = rec.annotations["date"]
        if m := DATE.match(old.strip()):
            mon = m.group(2).upper()
            mon = MONTH_ALIASES.get(mon, mon[:3])
            if mon in MONTHS:
                new = f"{int(m.group(1)):02d}-{mon}-{m.group(3)}"
                if new != old:
                    rec.annotations["date"] = new
                    res.changes.append(
                        Change("english_dates", "edited", "LOCUS", f"{old} -> {new}")
                    )

    rec.features = feats
    res.n_after = len(feats)
    if str(rec.seq) != seq_before or rec.annotations.get("topology") != topo_before:
        raise RuntimeError(
            "cleaner changed the sequence or topology - refusing to continue"
        )
    return res


def read_genbank(path) -> object:
    with open(path, "rb") as fh:
        text, _ = decode(fh.read())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return next(SeqIO.parse(io.StringIO(text), "genbank"))
