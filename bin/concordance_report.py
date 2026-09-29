#!/usr/bin/env python3
"""
concordance_report.py - turn somalier output into sample-identity verdicts.

somalier reports what the data says (relatedness, IBS0, depth, contamination).
This reduces that to a decision against the identity structure declared in the
samplesheet: which pairs were supposed to be the same individual, which were
supposed to be different, and which of those expectations the data violates.

Decision statistic
------------------
A pair is IDENTICAL when relatedness >= identity_min AND the IBS0 rate is at or
below identity_max_ibs0_rate.

Both conditions matter. IBS0 counts sites where one sample is hom-ref and the
other hom-alt, which is genetically impossible within one individual, so it is
what separates self (relatedness ~1, IBS0 ~0) from parent-child (relatedness
~0.5, IBS0 ~0) and from full sibs (relatedness ~0.5, IBS0 > 0). It is also the
statistic that survives tumour LOH: loss of heterozygosity turns hets into homs
and depresses het-based concordance for a genuine same-patient pair, but it
cannot manufacture a hom-ref/hom-alt transition.

Raw genotype concordance percentage is deliberately not the decision statistic.
It is uncalibrated and both MAF- and depth-dependent, so a threshold tuned on
one panel or one coverage regime does not transfer.

Stdlib only, Python 3.8+.
"""

import argparse
import csv
import html
import json
import os
import sys
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# observed relationship classes
# ---------------------------------------------------------------------------
IDENTICAL     = "IDENTICAL"
PARENT_CHILD  = "PARENT_CHILD"
SIBLING       = "SIBLING"
SECOND_DEGREE = "SECOND_DEGREE"
UNRELATED     = "UNRELATED"
INCONCLUSIVE  = "INCONCLUSIVE"

RELATED_CLASSES = (PARENT_CHILD, SIBLING, SECOND_DEGREE)

# expectation derived from the samplesheet
EXP_SAME_SUBJECT = "same_subject"
EXP_SAME_FAMILY  = "same_family"
EXP_UNRELATED    = "unrelated"
EXP_UNKNOWN      = "unknown"

PASS, WARN, FAIL = "pass", "warn", "fail"


# ---------------------------------------------------------------------------
# parsing helpers
# ---------------------------------------------------------------------------
def read_tsv(path):
    """Read a somalier TSV. The header line may be prefixed with '#'."""
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path, "r", newline="") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if header is None:
                header = [c.strip() for c in line.lstrip("#").split("\t")]
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts[: len(header)])))
    return rows


def as_float(value, default=0.0):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return default
    if f != f or f in (float("inf"), float("-inf")):  # NaN / inf
        return default
    return f


def as_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def pick(row, *names, **kwargs):
    """First present, non-empty column from a list of aliases."""
    default = kwargs.get("default", "")
    for n in names:
        if n in row and row[n] != "":
            return row[n]
    return default


def normalise_sex(value):
    v = (value or "").strip().lower()
    if v in ("m", "male", "1"):
        return "male"
    if v in ("f", "female", "2"):
        return "female"
    return "unknown"


# ---------------------------------------------------------------------------
# samplesheet -> expectations
# ---------------------------------------------------------------------------
def read_samplesheet(path):
    """sample -> {subject, family, sex, batch}. Keyed on the `sample` column."""
    declared = {}
    with open(path, "r", newline="") as fh:
        for row in csv.DictReader(fh):
            row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
            sample = row.get("sample")
            if not sample:
                continue
            subject = row.get("subject") or sample
            declared[sample] = {
                "sample":  sample,
                "subject": subject,
                # a row with no family is its own family: relatedness to any
                # other subject is then unexpected and gets flagged
                "family":  row.get("family") or subject,
                "sex":     normalise_sex(row.get("sex")),
                "batch":   row.get("batch") or "default",
            }
    return declared


def expectation_for(a, b, declared):
    da, db = declared.get(a), declared.get(b)
    if da is None or db is None:
        return EXP_UNKNOWN
    if da["subject"] == db["subject"]:
        return EXP_SAME_SUBJECT
    if da["family"] == db["family"]:
        return EXP_SAME_FAMILY
    return EXP_UNRELATED


# ---------------------------------------------------------------------------
# classification and verdicts
# ---------------------------------------------------------------------------
def classify_pair(relatedness, ibs0_rate, n_sites, t):
    if n_sites < t["min_sites"]:
        return INCONCLUSIVE
    if relatedness >= t["identity_min_relatedness"] and ibs0_rate <= t["identity_max_ibs0_rate"]:
        return IDENTICAL
    if relatedness >= t["first_degree_min_relatedness"]:
        # both first-degree; IBS0 splits parent-child (~0) from full sibs (>0)
        return PARENT_CHILD if ibs0_rate <= t["identity_max_ibs0_rate"] else SIBLING
    if relatedness >= t["second_degree_min_relatedness"]:
        return SECOND_DEGREE
    return UNRELATED


def verdict_for(expectation, observed):
    """(verdict, severity, note)"""
    if observed == INCONCLUSIVE:
        return ("INCONCLUSIVE", WARN,
                "too few usable sites to decide; check depth and that the sites file "
                "matches the reference contig naming")

    if expectation == EXP_SAME_SUBJECT:
        if observed == IDENTICAL:
            return ("CONFIRMED_MATCH", PASS, "")
        return ("SAMPLE_SWAP", FAIL,
                "declared as the same subject but the genotypes disagree")

    if expectation == EXP_SAME_FAMILY:
        if observed == IDENTICAL:
            return ("UNEXPECTED_DUPLICATE", FAIL,
                    "declared as different subjects in one family but the genotypes are identical")
        if observed in RELATED_CLASSES:
            return ("CONFIRMED_RELATED", PASS, "")
        return ("MISSING_RELATEDNESS", WARN,
                "declared as the same family but no relatedness detected")

    if expectation == EXP_UNRELATED:
        if observed == IDENTICAL:
            return ("UNEXPECTED_DUPLICATE", FAIL,
                    "declared as different subjects but the genotypes are identical")
        if observed in RELATED_CLASSES:
            return ("UNEXPECTED_RELATEDNESS", WARN,
                    "declared as unrelated but relatedness detected; check for cross-contamination "
                    "before assuming a pedigree error")
        return ("CONFIRMED_DISTINCT", PASS, "")

    # one or both sides come from the sketch registry, so there is no declared
    # expectation. A match against history is the signal worth surfacing.
    if observed == IDENTICAL:
        return ("REGISTRY_MATCH", WARN,
                "matches a sketch from the registry that is not in this samplesheet")
    if observed in RELATED_CLASSES:
        return ("REGISTRY_RELATED", WARN, "related to a registry sketch")
    return ("NOT_REPORTED", PASS, "")


# ---------------------------------------------------------------------------
# per-sample QC
# ---------------------------------------------------------------------------
def infer_sex(row):
    """
    Prefer somalier's own inference (the `sex` column, written when relate is
    run with --infer: 1=male, 2=female). Fall back to X heterozygosity and
    relative Y depth when that column is unset.
    """
    declared_by_somalier = (pick(row, "sex") or "").strip()
    if declared_by_somalier == "1":
        return "male", "somalier --infer"
    if declared_by_somalier == "2":
        return "female", "somalier --infer"

    x_hom_ref = as_int(pick(row, "X_hom_ref"))
    x_het     = as_int(pick(row, "X_het"))
    x_hom_alt = as_int(pick(row, "X_hom_alt"))
    x_total   = x_hom_ref + x_het + x_hom_alt
    if x_total < 10:
        return "unknown", "too few X sites"

    x_het_rate = x_het / float(x_total)
    gt_depth   = as_float(pick(row, "gt_depth_mean"))
    y_depth    = as_float(pick(row, "Y_depth_mean"))
    y_ratio    = (y_depth / gt_depth) if gt_depth > 0 else 0.0

    if x_het_rate < 0.10 and y_ratio > 0.10:
        return "male", "X het %.3f, Y/auto depth %.2f" % (x_het_rate, y_ratio)
    if x_het_rate > 0.20 and y_ratio < 0.05:
        return "female", "X het %.3f, Y/auto depth %.2f" % (x_het_rate, y_ratio)
    return "ambiguous", "X het %.3f, Y/auto depth %.2f" % (x_het_rate, y_ratio)


def build_sample_table(sample_rows, declared, contam_samples, contam_pairs, ancestry, t):
    # top contaminator per sample, from the anchored pairwise MLE
    worst_contaminator = {}
    for row in contam_pairs:
        sample = pick(row, "sample_name")
        anchor = pick(row, "anchor_sample")
        mle    = as_float(pick(row, "contamination_mle"))
        if not sample:
            continue
        current = worst_contaminator.get(sample)
        if current is None or mle > current[1]:
            worst_contaminator[sample] = (anchor, mle)

    charr_by_sample = {}
    for row in contam_samples:
        name = pick(row, "sample_name")
        if name:
            charr_by_sample[name] = as_float(pick(row, "contamination_charr"))

    ancestry_by_sample = {}
    for row in ancestry:
        name = pick(row, "#sample_id", "sample_id", "sample")
        if name:
            ancestry_by_sample[name] = pick(row, "predicted_ancestry", "given_ancestry", default="")

    out = []
    for row in sample_rows:
        sid = pick(row, "sample_id")
        if not sid:
            continue
        decl = declared.get(sid)

        n_sites = (as_int(pick(row, "n_hom_ref")) + as_int(pick(row, "n_het"))
                   + as_int(pick(row, "n_hom_alt")))
        depth   = as_float(pick(row, "gt_depth_mean"))

        inferred, sex_evidence = infer_sex(row)
        declared_sex = decl["sex"] if decl else "unknown"

        # charr from the contamination subcommand if present, else from relate
        charr = charr_by_sample.get(sid)
        if charr is None:
            charr = as_float(pick(row, "contamination_charr"), default=-1.0)

        notes, status = [], PASS

        def escalate(level):
            order = {PASS: 0, WARN: 1, FAIL: 2}
            return level if order[level] > order[status] else status

        if n_sites < t["min_sites"]:
            notes.append("only %d usable sites (< %d)" % (n_sites, t["min_sites"]))
            status = escalate(FAIL)
        if depth < t["min_depth"]:
            notes.append("mean genotype depth %.1f (< %.1f)" % (depth, t["min_depth"]))
            status = escalate(WARN)

        if declared_sex != "unknown" and inferred in ("male", "female") and inferred != declared_sex:
            notes.append("sex mismatch: declared %s, inferred %s (%s)"
                         % (declared_sex, inferred, sex_evidence))
            status = escalate(FAIL)
        elif inferred == "ambiguous":
            notes.append("sex ambiguous (%s)" % sex_evidence)
            status = escalate(WARN)

        if charr >= 0:
            if charr >= t["contamination_fail"]:
                notes.append("contamination %.3f (>= %.3f)" % (charr, t["contamination_fail"]))
                status = escalate(FAIL)
            elif charr >= t["contamination_warn"]:
                notes.append("contamination %.3f (>= %.3f)" % (charr, t["contamination_warn"]))
                status = escalate(WARN)

        # The anchored pairwise MLE only means something once there is
        # contamination to attribute. Below the warn threshold it is noise, and
        # a value pinned at 1.0 is the optimiser hitting its bound rather than a
        # real estimate - in both cases reporting a named contaminator would
        # read as a finding when it is not one.
        top_src, top_mle = worst_contaminator.get(sid, ("", 0.0))
        if charr < t["contamination_warn"]:
            top_src, top_mle = "", 0.0
        elif top_src and top_mle >= 0.99:
            notes.append("pairwise contamination estimate saturated at %.2f against %s; "
                         "treat the attribution as unreliable" % (top_mle, top_src))
            top_src = top_src + " (saturated)"

        out.append({
            "sample_id":            sid,
            "subject":              decl["subject"] if decl else "",
            "family":               decl["family"] if decl else "",
            "batch":                decl["batch"] if decl else "",
            "in_samplesheet":       "yes" if decl else "no (registry)",
            "n_sites":              n_sites,
            "gt_depth_mean":        round(depth, 2),
            "declared_sex":         declared_sex,
            "inferred_sex":         inferred,
            "sex_evidence":         sex_evidence,
            "contamination_charr":  round(charr, 4) if charr >= 0 else "",
            "top_contaminator":     top_src,
            "top_contamination_mle": round(top_mle, 4) if top_src else "",
            "ancestry":             ancestry_by_sample.get(sid, ""),
            "status":               status,
            "notes":                "; ".join(notes),
        })
    return out


# ---------------------------------------------------------------------------
# pair analysis
# ---------------------------------------------------------------------------
def build_pair_table(pair_rows, declared, t):
    results = []
    for row in pair_rows:
        a = pick(row, "sample_a")
        b = pick(row, "sample_b")
        if not a or not b:
            continue

        relatedness = as_float(pick(row, "relatedness"))
        ibs0        = as_int(pick(row, "ibs0"))
        ibs2        = as_int(pick(row, "ibs2"))
        n_sites     = as_int(pick(row, "n"))
        # 0.3.2 renamed homozygous_concordance -> concordance
        conc        = as_float(pick(row, "concordance", "homozygous_concordance"), default=-1.0)
        ibs0_rate   = (ibs0 / float(n_sites)) if n_sites > 0 else 1.0

        expectation = expectation_for(a, b, declared)
        observed    = classify_pair(relatedness, ibs0_rate, n_sites, t)
        verdict, severity, note = verdict_for(expectation, observed)

        results.append({
            "sample_a":         a,
            "sample_b":         b,
            "subject_a":        declared.get(a, {}).get("subject", ""),
            "subject_b":        declared.get(b, {}).get("subject", ""),
            "expectation":      expectation,
            "relatedness":      round(relatedness, 4),
            "ibs0":             ibs0,
            "ibs0_rate":        round(ibs0_rate, 6),
            "ibs2":             ibs2,
            "concordance":      round(conc, 4) if conc >= 0 else "",
            "shared_hets":      as_int(pick(row, "shared_hets")),
            "shared_hom_alts":  as_int(pick(row, "shared_hom_alts")),
            "n_sites":          n_sites,
            "observed_class":   observed,
            "verdict":          verdict,
            "severity":         severity,
            "note":             note,
        })
    return results


def is_reportable(pair):
    """
    Bound the output. Every declared pair is reported; registry pairs only when
    they are not plainly unrelated, otherwise a large registry drowns the table
    in N-squared uninteresting rows.
    """
    if pair["expectation"] != EXP_UNKNOWN:
        return True
    return pair["verdict"] != "NOT_REPORTED"


# ---------------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------------
PAIR_COLUMNS = [
    "sample_a", "sample_b", "subject_a", "subject_b", "expectation",
    "relatedness", "ibs0", "ibs0_rate", "ibs2", "concordance",
    "shared_hets", "shared_hom_alts", "n_sites",
    "observed_class", "verdict", "severity", "note",
]

SAMPLE_COLUMNS = [
    "sample_id", "subject", "family", "batch", "in_samplesheet",
    "n_sites", "gt_depth_mean", "declared_sex", "inferred_sex", "sex_evidence",
    "contamination_charr", "top_contaminator", "top_contamination_mle",
    "ancestry", "status", "notes",
]


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_multiqc(path, columns, rows, mqc_id, section, description):
    with open(path, "w", newline="") as fh:
        fh.write("# id: '%s'\n" % mqc_id)
        fh.write("# section_name: '%s'\n" % section)
        fh.write("# description: '%s'\n" % description)
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
CSS = """
:root {
  --bg:#ffffff; --fg:#1a1d21; --muted:#5b6672; --line:#e3e8ee; --panel:#f7f9fc;
  --pass:#1a7f52; --pass-bg:#e7f5ee; --warn:#8a6100; --warn-bg:#fdf3dd;
  --fail:#a4262c; --fail-bg:#fdeaea; --accent:#2b5cb8;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg:#15181c; --fg:#e6e9ed; --muted:#9aa5b1; --line:#2b3138; --panel:#1c2026;
    --pass:#4cc38a; --pass-bg:#16281f; --warn:#e2b35c; --warn-bg:#2b2418;
    --fail:#f2777a; --fail-bg:#2e1b1c; --accent:#7aa2f7;
  }
}
* { box-sizing:border-box; }
body { margin:0; padding:2rem 1.5rem 4rem; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; }
.wrap { max-width:1240px; margin:0 auto; }
h1 { font-size:1.5rem; margin:0 0 .25rem; letter-spacing:-.01em; }
h2 { font-size:1.05rem; margin:2.5rem 0 .75rem; padding-bottom:.4rem;
  border-bottom:1px solid var(--line); letter-spacing:-.01em; }
.sub { color:var(--muted); margin:0 0 1.5rem; font-size:.875rem; }
.cards { display:flex; flex-wrap:wrap; gap:.75rem; margin:1.25rem 0; }
.card { flex:1 1 140px; border:1px solid var(--line); border-radius:8px;
  padding:.85rem 1rem; background:var(--panel); }
.card .n { font-size:1.6rem; font-weight:650; line-height:1.1; font-variant-numeric:tabular-nums; }
.card .l { color:var(--muted); font-size:.75rem; text-transform:uppercase;
  letter-spacing:.06em; margin-top:.2rem; }
.card.pass .n { color:var(--pass); } .card.warn .n { color:var(--warn); }
.card.fail .n { color:var(--fail); }
.scroll { overflow-x:auto; border:1px solid var(--line); border-radius:8px; }
table { border-collapse:collapse; width:100%; font-size:.8125rem; }
th, td { text-align:left; padding:.45rem .65rem; border-bottom:1px solid var(--line);
  white-space:nowrap; }
th { background:var(--panel); font-weight:600; font-size:.7rem; text-transform:uppercase;
  letter-spacing:.05em; color:var(--muted); position:sticky; top:0; }
tr:last-child td { border-bottom:none; }
td.num { text-align:right; font-variant-numeric:tabular-nums; }
td.wrap-note { white-space:normal; min-width:260px; color:var(--muted); }
.tag { display:inline-block; padding:.1rem .45rem; border-radius:4px; font-size:.7rem;
  font-weight:650; letter-spacing:.02em; }
.tag.pass { color:var(--pass); background:var(--pass-bg); }
.tag.warn { color:var(--warn); background:var(--warn-bg); }
.tag.fail { color:var(--fail); background:var(--fail-bg); }
.empty { padding:1.1rem; color:var(--muted); background:var(--panel);
  border:1px solid var(--line); border-radius:8px; }
details { margin-top:.75rem; }
summary { cursor:pointer; color:var(--accent); font-weight:600; font-size:.8125rem;
  padding:.35rem 0; }
.matrix td.cell { text-align:center; font-variant-numeric:tabular-nums; font-size:.7rem;
  min-width:52px; }
.matrix th.rowhead { position:sticky; left:0; background:var(--panel); z-index:1; }
footer { margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line);
  color:var(--muted); font-size:.75rem; }
code { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.9em;
  background:var(--panel); padding:.05rem .3rem; border-radius:3px; }
"""


def esc(v):
    return html.escape(str(v), quote=True)


def tag(severity):
    return '<span class="tag %s">%s</span>' % (severity, esc(severity.upper()))


def matrix_cell_style(rel):
    """Shade by relatedness: white -> accent. Pure inline, no JS."""
    if rel is None:
        return "", ""
    r = max(0.0, min(1.0, rel))
    alpha = 0.06 + 0.78 * r
    return ("background: rgba(43,92,184,%.3f);" % alpha), ("%.2f" % rel)


def render_table(columns, rows, numeric=(), notes=(), severity_col=None):
    out = ['<div class="scroll"><table><thead><tr>']
    for c in columns:
        out.append("<th>%s</th>" % esc(c.replace("_", " ")))
    out.append("</tr></thead><tbody>")
    for r in rows:
        out.append("<tr>")
        for c in columns:
            v = r.get(c, "")
            if c == severity_col:
                out.append("<td>%s</td>" % tag(str(v)))
            elif c in notes:
                out.append('<td class="wrap-note">%s</td>' % esc(v))
            elif c in numeric:
                out.append('<td class="num">%s</td>' % esc(v))
            else:
                out.append("<td>%s</td>" % esc(v))
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def render_matrix(pairs, declared, limit=60):
    names = sorted(declared.keys())
    if len(names) < 2:
        return ""
    truncated = len(names) > limit
    names = names[:limit]
    index = {n: i for i, n in enumerate(names)}

    grid = [[None] * len(names) for _ in names]
    for p in pairs:
        a, b = p["sample_a"], p["sample_b"]
        if a in index and b in index:
            grid[index[a]][index[b]] = p["relatedness"]
            grid[index[b]][index[a]] = p["relatedness"]

    out = ['<div class="scroll"><table class="matrix"><thead><tr><th class="rowhead"></th>']
    for n in names:
        out.append("<th>%s</th>" % esc(n))
    out.append("</tr></thead><tbody>")
    for i, n in enumerate(names):
        out.append('<tr><th class="rowhead">%s</th>' % esc(n))
        for j in range(len(names)):
            if i == j:
                out.append('<td class="cell" style="background:rgba(43,92,184,0.9);'
                           'color:#fff">1.00</td>')
            else:
                style, label = matrix_cell_style(grid[i][j])
                out.append('<td class="cell" style="%s">%s</td>' % (style, label))
        out.append("</tr>")
    out.append("</tbody></table></div>")
    if truncated:
        out.append('<p class="sub">Matrix truncated to the first %d samples; '
                   "the full pair table above is complete.</p>" % limit)
    return "".join(out)


def render_html(run_name, pairs, samples, summary, t):
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    action_pairs = [p for p in pairs if p["severity"] in (FAIL, WARN)]
    action_pairs.sort(key=lambda p: (0 if p["severity"] == FAIL else 1, -p["relatedness"]))
    action_samples = [s for s in samples if s["status"] in (FAIL, WARN)]
    action_samples.sort(key=lambda s: (0 if s["status"] == FAIL else 1, s["sample_id"]))

    pair_cols = ["sample_a", "sample_b", "expectation", "observed_class", "relatedness",
                 "ibs0_rate", "n_sites", "verdict", "severity", "note"]
    sample_cols = ["sample_id", "subject", "n_sites", "gt_depth_mean", "declared_sex",
                   "inferred_sex", "contamination_charr", "top_contaminator",
                   "status", "notes"]
    numeric = {"relatedness", "ibs0", "ibs0_rate", "ibs2", "n_sites", "concordance",
               "shared_hets", "shared_hom_alts", "gt_depth_mean",
               "contamination_charr", "top_contamination_mle"}

    if action_pairs:
        action_block = render_table(pair_cols, action_pairs, numeric=numeric,
                                    notes={"note"}, severity_col="severity")
    else:
        action_block = ('<div class="empty">Every evaluated pair agrees with the identity '
                        "structure declared in the samplesheet.</div>")

    if action_samples:
        sample_action_block = render_table(sample_cols, action_samples, numeric=numeric,
                                           notes={"notes"}, severity_col="status")
    else:
        sample_action_block = '<div class="empty">No per-sample QC flags.</div>'

    matrix = render_matrix(pairs, {s["sample_id"]: s for s in samples
                                   if s["in_samplesheet"] == "yes"})
    matrix_block = matrix if matrix else '<div class="empty">Not enough samples to plot.</div>'

    parts = []
    parts.append("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    parts.append("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    parts.append("<title>Sample concordance - %s</title>" % esc(run_name))
    parts.append("<style>%s</style></head><body><div class='wrap'>" % CSS)

    parts.append("<h1>Sample concordance</h1>")
    parts.append("<p class='sub'>%s &middot; generated %s &middot; somalier sketches, "
                 "all-vs-all</p>" % (esc(run_name), esc(generated)))

    parts.append("<div class='cards'>")
    for label, key, cls in [
        ("Samples", "n_samples", ""),
        ("Pairs evaluated", "n_pairs_evaluated", ""),
        ("Pass", "n_pass", "pass"),
        ("Warn", "n_warn", "warn"),
        ("Fail", "n_fail", "fail"),
    ]:
        parts.append("<div class='card %s'><div class='n'>%s</div>"
                     "<div class='l'>%s</div></div>" % (cls, summary[key], esc(label)))
    parts.append("</div>")

    parts.append("<h2>Action required &mdash; pairs</h2>")
    parts.append(action_block)

    parts.append("<h2>Action required &mdash; samples</h2>")
    parts.append(sample_action_block)

    parts.append("<h2>Relatedness matrix</h2>")
    parts.append("<p class='sub'>Declared samples only. Shading is relatedness; "
                 "the diagonal is self.</p>")
    parts.append(matrix_block)

    parts.append("<h2>All samples</h2>")
    parts.append(render_table(sample_cols, sorted(samples, key=lambda s: s["sample_id"]),
                              numeric=numeric, notes={"notes"}, severity_col="status"))

    parts.append("<h2>All evaluated pairs</h2>")
    parts.append("<details><summary>Show %d pairs</summary>%s</details>" % (
        len(pairs),
        render_table(pair_cols, sorted(pairs, key=lambda p: -p["relatedness"]),
                     numeric=numeric, notes={"note"}, severity_col="severity")))

    parts.append("<footer>")
    parts.append("<p><strong>Decision rule.</strong> A pair is called identical when "
                 "relatedness &ge; <code>%s</code> <em>and</em> the IBS0 rate is "
                 "&le; <code>%s</code>. IBS0 counts sites where one sample is hom-ref and "
                 "the other hom-alt &mdash; impossible within one individual &mdash; so it is "
                 "what separates self from parent-child, and it survives tumour LOH where "
                 "het-based concordance does not.</p>"
                 % (t["identity_min_relatedness"], t["identity_max_ibs0_rate"]))
    parts.append("<p>Thresholds: first-degree &ge; <code>%s</code>, second-degree &ge; "
                 "<code>%s</code>, min sites <code>%s</code>, min depth <code>%s</code>, "
                 "contamination warn/fail <code>%s</code>/<code>%s</code>.</p>"
                 % (t["first_degree_min_relatedness"], t["second_degree_min_relatedness"],
                    t["min_sites"], t["min_depth"],
                    t["contamination_warn"], t["contamination_fail"]))
    parts.append("<p>Relatedness, IBS0 and contamination estimates are produced by "
                 "somalier (Pedersen &amp; Quinlan, <em>Genome Medicine</em> 2020, "
                 "doi:10.1186/s13073-020-00761-2). This page only applies decision "
                 "thresholds to them.</p>")
    parts.append("</footer></div></body></html>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="Sample-identity verdicts from somalier output")
    ap.add_argument("--pairs", required=True, help="somalier relate *.pairs.tsv")
    ap.add_argument("--samples", required=True, help="somalier relate *.samples.tsv")
    ap.add_argument("--samplesheet", required=True, help="pipeline samplesheet CSV")
    ap.add_argument("--contamination-pairs", default=None)
    ap.add_argument("--contamination-samples", default=None)
    ap.add_argument("--ancestry", default=None)
    ap.add_argument("--run-name", default="cohort")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--min-sites", type=int, default=1000)
    ap.add_argument("--min-depth", type=float, default=5.0)
    ap.add_argument("--identity-min-relatedness", type=float, default=0.85)
    ap.add_argument("--identity-max-ibs0-rate", type=float, default=0.005)
    ap.add_argument("--first-degree-min-relatedness", type=float, default=0.35)
    ap.add_argument("--second-degree-min-relatedness", type=float, default=0.15)
    ap.add_argument("--contamination-warn", type=float, default=0.02)
    ap.add_argument("--contamination-fail", type=float, default=0.05)
    ap.add_argument("--fail-on-discordance", action="store_true",
                    help="exit 1 when any pair or sample verdict is FAIL")
    args = ap.parse_args(argv)

    t = {
        "min_sites":                     args.min_sites,
        "min_depth":                     args.min_depth,
        "identity_min_relatedness":      args.identity_min_relatedness,
        "identity_max_ibs0_rate":        args.identity_max_ibs0_rate,
        "first_degree_min_relatedness":  args.first_degree_min_relatedness,
        "second_degree_min_relatedness": args.second_degree_min_relatedness,
        "contamination_warn":            args.contamination_warn,
        "contamination_fail":            args.contamination_fail,
    }

    declared        = read_samplesheet(args.samplesheet)
    pair_rows       = read_tsv(args.pairs)
    sample_rows     = read_tsv(args.samples)
    contam_pairs    = read_tsv(args.contamination_pairs)
    contam_samples  = read_tsv(args.contamination_samples)
    ancestry        = read_tsv(args.ancestry)

    all_pairs = build_pair_table(pair_rows, declared, t)
    pairs     = [p for p in all_pairs if is_reportable(p)]
    samples   = build_sample_table(sample_rows, declared, contam_samples, contam_pairs,
                                   ancestry, t)

    n_fail = sum(1 for p in pairs if p["severity"] == FAIL) + \
             sum(1 for s in samples if s["status"] == FAIL)
    n_warn = sum(1 for p in pairs if p["severity"] == WARN) + \
             sum(1 for s in samples if s["status"] == WARN)
    n_pass = sum(1 for p in pairs if p["severity"] == PASS) + \
             sum(1 for s in samples if s["status"] == PASS)

    # samples declared in the sheet that somalier never produced a sketch for
    seen = {s["sample_id"] for s in samples}
    missing = sorted(set(declared.keys()) - seen)

    summary = {
        "run_name":            args.run_name,
        "generated_utc":       datetime.now(timezone.utc).isoformat(),
        "n_samples":           len(samples),
        "n_samples_declared":  len(declared),
        "n_samples_registry":  sum(1 for s in samples if s["in_samplesheet"] != "yes"),
        "n_pairs_total":       len(all_pairs),
        "n_pairs_evaluated":   len(pairs),
        "n_pass":              n_pass,
        "n_warn":              n_warn,
        "n_fail":              n_fail,
        "missing_samples":     missing,
        "thresholds":          t,
        "verdict_counts":      _count_by(pairs, "verdict"),
        "sample_status_counts": _count_by(samples, "status"),
        "failures": [
            {k: p[k] for k in ("sample_a", "sample_b", "verdict", "relatedness",
                               "ibs0_rate", "n_sites", "note")}
            for p in pairs if p["severity"] == FAIL
        ],
        "sample_failures": [
            {k: s[k] for k in ("sample_id", "subject", "status", "notes")}
            for s in samples if s["status"] == FAIL
        ],
    }

    outdir = args.outdir
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    def out(name):
        return os.path.join(outdir, name)

    write_tsv(out("concordance.pairs.tsv"), PAIR_COLUMNS, pairs)
    write_tsv(out("concordance.samples.tsv"), SAMPLE_COLUMNS, samples)
    with open(out("concordance.summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(out("concordance_report.html"), "w") as fh:
        fh.write(render_html(args.run_name, pairs, samples, summary, t))

    write_multiqc(out("concordance_pairs_mqc.tsv"),
                  ["sample_a", "sample_b", "expectation", "observed_class",
                   "relatedness", "ibs0_rate", "n_sites", "verdict"],
                  [p for p in pairs if p["severity"] in (FAIL, WARN)] or pairs[:50],
                  "concordance_pairs", "Sample concordance - pairs",
                  "Pairwise identity verdicts against the declared samplesheet structure.")
    write_multiqc(out("concordance_samples_mqc.tsv"),
                  ["sample_id", "subject", "n_sites", "gt_depth_mean", "declared_sex",
                   "inferred_sex", "contamination_charr", "status"],
                  samples, "concordance_samples", "Sample concordance - samples",
                  "Per-sample identity QC: site count, depth, sex check, contamination.")

    # stdout summary - this is what shows up in .nextflow.log and CI output
    print("[concordance] %d samples (%d declared, %d from registry), %d pairs evaluated"
          % (summary["n_samples"], summary["n_samples_declared"],
             summary["n_samples_registry"], summary["n_pairs_evaluated"]))
    if missing:
        print("[concordance] WARNING: declared but no sketch produced: %s" % ", ".join(missing))
    for p in pairs:
        if p["severity"] == FAIL:
            print("[concordance] FAIL %s <-> %s : %s (relatedness %.3f, IBS0 rate %.4f)"
                  % (p["sample_a"], p["sample_b"], p["verdict"],
                     p["relatedness"], p["ibs0_rate"]))
    for s in samples:
        if s["status"] == FAIL:
            print("[concordance] FAIL %s : %s" % (s["sample_id"], s["notes"]))
    print("[concordance] pass=%d warn=%d fail=%d" % (n_pass, n_warn, n_fail))

    if args.fail_on_discordance and n_fail > 0:
        print("[concordance] --fail-on-discordance set and %d failures found" % n_fail,
              file=sys.stderr)
        return 1
    return 0


def _count_by(rows, key):
    counts = {}
    for r in rows:
        counts[r[key]] = counts.get(r[key], 0) + 1
    return counts


if __name__ == "__main__":
    sys.exit(main())
