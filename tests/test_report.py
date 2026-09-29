#!/usr/bin/env python3
"""
Regression test for bin/concordance_report.py.

Builds synthetic somalier output that exercises every verdict path, runs the
report, and asserts the verdicts. No somalier, no containers, no data needed:

    python3 tests/test_report.py
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "..", "bin", "concordance_report.py")

PAIR_HEADER = ("#sample_a\tsample_b\trelatedness\tibs0\tibs2\tconcordance\thets_a\thets_b\t"
               "hets_ab\tshared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\t"
               "x_ibs0\tx_ibs2\texpected_relatedness")

SAMPLE_HEADER = ("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t"
                 "original_pedigree_sex\tgt_depth_mean\tgt_depth_sd\tdepth_mean\tdepth_sd\t"
                 "ab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\t"
                 "contamination_charr\tcontamination_charr_n_sites\t"
                 "X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n")

N = 17000

# sample, subject, family, sex declared in the samplesheet
SAMPLESHEET = [
    ("S1", "P1", "F1", "female"),
    ("S2", "P1", "F1", "female"),   # true duplicate of S1
    ("S3", "P2", "F1", "male"),     # P1's child, same family
    ("S4", "P3", "F3", "male"),
    ("S5", "P4", "F4", "female"),   # data says male -> sex mismatch
    ("S6", "P5", "F5", "male"),     # S6/S7 declared same subject...
    ("S7", "P5", "F5", "male"),     # ...but genotypes disagree -> swap
    ("S8", "P6", "F6", "male"),     # identical to S4 -> unexpected duplicate
    ("S9", "P7", "F7", "male"),     # contaminated
]

# (a, b, relatedness, ibs0, n)
PAIRS = [
    ("S1", "S2", 0.992, 4,    N),      # CONFIRMED_MATCH
    ("S1", "S3", 0.486, 6,    N),      # CONFIRMED_RELATED (parent-child, same family)
    ("S2", "S3", 0.478, 8,    N),      # CONFIRMED_RELATED
    ("S6", "S7", 0.013, 1850, N),      # SAMPLE_SWAP
    ("S4", "S8", 0.981, 3,    N),      # UNEXPECTED_DUPLICATE (different families)
    ("S1", "S4", 0.008, 1920, N),      # CONFIRMED_DISTINCT
    ("S1", "S5", 0.004, 6,    40),     # INCONCLUSIVE (too few sites)
    ("S3", "S9", 0.210, 320,  N),      # UNEXPECTED_RELATEDNESS (2nd degree, diff family)
    ("S1", "R1", 0.988, 5,    N),      # REGISTRY_MATCH (R1 not in the samplesheet)
    ("S2", "S4", 0.006, 1890, N),
    ("S5", "S9", 0.002, 1955, N),
]

# sample -> (somalier_sex_code, gt_depth, n_hom_ref, n_het, n_hom_alt, charr,
#            X_hom_ref, X_het, X_hom_alt, Y_depth)
SAMPLES = {
    "S1": ("2", 32.0, 9000, 5500, 2500, 0.001, 300, 140, 120, 0.1),
    "S2": ("2", 30.0, 9010, 5480, 2510, 0.002, 302, 138, 118, 0.1),
    "S3": ("1", 31.0, 9100, 5400, 2500, 0.001, 520,   4,  40, 15.0),
    "S4": ("1", 28.0, 9050, 5450, 2500, 0.003, 515,   6,  41, 14.2),
    "S5": ("1", 29.0, 9080, 5420, 2500, 0.001, 518,   5,  39, 14.8),   # declared female
    "S6": ("1", 27.0, 9020, 5470, 2510, 0.002, 512,   7,  42, 13.9),
    "S7": ("1", 26.0, 9040, 5460, 2500, 0.001, 514,   6,  40, 14.1),
    "S8": ("1", 25.0, 9060, 5440, 2500, 0.002, 516,   5,  41, 14.0),
    "S9": ("1", 30.0, 9030, 5465, 2505, 0.071, 513,   6,  41, 14.3),   # contaminated
    "R1": ("2", 33.0, 9005, 5495, 2500, 0.001, 301, 139, 119, 0.1),    # registry sketch
    # a declared sample with far too few sites -> per-sample FAIL
    "S10": ("2", 2.0,   200,   90,   40, 0.001,  10,   4,   3, 0.05),
}

EXPECTED_PAIR_VERDICTS = {
    ("S1", "S2"): "CONFIRMED_MATCH",
    ("S1", "S3"): "CONFIRMED_RELATED",
    ("S2", "S3"): "CONFIRMED_RELATED",
    ("S6", "S7"): "SAMPLE_SWAP",
    ("S4", "S8"): "UNEXPECTED_DUPLICATE",
    ("S1", "S4"): "CONFIRMED_DISTINCT",
    ("S1", "S5"): "INCONCLUSIVE",
    ("S3", "S9"): "UNEXPECTED_RELATEDNESS",
    ("S1", "R1"): "REGISTRY_MATCH",
}


def write_fixtures(d):
    sheet = os.path.join(d, "samplesheet.csv")
    with open(sheet, "w") as fh:
        fh.write("sample,subject,family,alignment,index,sex,batch\n")
        for s, subj, fam, sex in SAMPLESHEET:
            fh.write("%s,%s,%s,/dev/null/%s.cram,,%s,run1\n" % (s, subj, fam, s, sex))
        # declared but its sketch comes back nearly empty
        fh.write("S10,P8,F8,/dev/null/S10.cram,,female,run1\n")

    pairs = os.path.join(d, "somalier.pairs.tsv")
    with open(pairs, "w") as fh:
        fh.write(PAIR_HEADER + "\n")
        for a, b, rel, ibs0, n in PAIRS:
            ibs2 = n - ibs0
            fh.write("\t".join([a, b, "%.3f" % rel, str(ibs0), str(ibs2), "0.950",
                                "5500", "5480", "5400", "5200", "2500", "2510",
                                "2400", str(n), "0", "500", "-1"]) + "\n")

    samples = os.path.join(d, "somalier.samples.tsv")
    with open(samples, "w") as fh:
        fh.write(SAMPLE_HEADER + "\n")
        for sid, (sex, depth, hr, het, ha, charr, xhr, xhet, xha, ydepth) in SAMPLES.items():
            fh.write("\t".join([
                "fam_" + sid, sid, "-9", "-9", sex, "-9", "-9",
                "%.1f" % depth, "3.0", "%.1f" % depth, "3.0", "0.50", "0.05",
                str(hr), str(het), str(ha), "0", "0.01",
                "%.4f" % charr, "12000",
                "%.1f" % (depth / 2.0), "560", str(xhr), str(xhet), str(xha),
                "%.1f" % ydepth, "300",
            ]) + "\n")

    contam = os.path.join(d, "somalier.contamination.pairs.tsv")
    with open(contam, "w") as fh:
        fh.write("#sample_name\tanchor_sample\tn_sites_usable\tcontamination_mle\n")
        fh.write("S9\tS3\t12000\t0.0710\n")
        fh.write("S9\tS4\t12000\t0.0120\n")

    contam_s = os.path.join(d, "somalier.contamination.samples.tsv")
    with open(contam_s, "w") as fh:
        fh.write("#sample_name\tn_sites_usable\tcontamination_charr\n")
        for sid, vals in SAMPLES.items():
            fh.write("%s\t12000\t%.4f\n" % (sid, vals[5]))

    return sheet, pairs, samples, contam, contam_s


def main():
    failures = []
    with tempfile.TemporaryDirectory() as d:
        sheet, pairs, samples, contam, contam_s = write_fixtures(d)
        out = os.path.join(d, "out")
        cmd = [sys.executable, REPORT,
               "--pairs", pairs, "--samples", samples, "--samplesheet", sheet,
               "--contamination-pairs", contam, "--contamination-samples", contam_s,
               "--run-name", "unit-test", "--outdir", out,
               "--min-sites", "1000", "--min-depth", "5"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        print(proc.stdout)
        if proc.returncode != 0:
            print(proc.stderr)
            failures.append("report exited %d" % proc.returncode)

        # ---- pair verdicts ----
        got = {}
        with open(os.path.join(out, "concordance.pairs.tsv")) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                got[(row["sample_a"], row["sample_b"])] = row["verdict"]

        for key, want in EXPECTED_PAIR_VERDICTS.items():
            have = got.get(key)
            status = "ok " if have == want else "FAIL"
            print("%s %-12s expected %-24s got %s" % (status, "%s/%s" % key, want, have))
            if have != want:
                failures.append("pair %s: expected %s, got %s" % (key, want, have))

        # registry pairs that are plainly unrelated must not be reported
        if ("S2", "S4") in got and got[("S2", "S4")] == "NOT_REPORTED":
            failures.append("NOT_REPORTED rows should be filtered out")

        # ---- sample verdicts ----
        sample_status, sample_notes = {}, {}
        with open(os.path.join(out, "concordance.samples.tsv")) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                sample_status[row["sample_id"]] = row["status"]
                sample_notes[row["sample_id"]] = row["notes"]

        checks = [
            ("S5", "fail", "sex mismatch"),
            ("S9", "fail", "contamination"),
            ("S10", "fail", "usable sites"),
            ("S1", "pass", ""),
            ("S3", "pass", ""),
        ]
        for sid, want_status, want_note in checks:
            have_status = sample_status.get(sid)
            have_note = sample_notes.get(sid, "")
            ok = have_status == want_status and want_note in have_note
            print("%s %-4s expected %-5s got %-5s  %s" % (
                "ok " if ok else "FAIL", sid, want_status, have_status, have_note))
            if not ok:
                failures.append("sample %s: expected %s/%r, got %s/%r"
                                % (sid, want_status, want_note, have_status, have_note))

        # R1 should appear, flagged as registry
        if sample_status.get("R1") is None:
            failures.append("registry sample R1 missing from the sample table")

        # ---- summary + html sanity ----
        with open(os.path.join(out, "concordance.summary.json")) as fh:
            summary = json.load(fh)
        if summary["n_fail"] < 2:
            failures.append("summary n_fail too low: %s" % summary["n_fail"])
        if summary["verdict_counts"].get("SAMPLE_SWAP") != 1:
            failures.append("summary should record exactly one SAMPLE_SWAP")

        html = open(os.path.join(out, "concordance_report.html")).read()
        for needle in ("<title>", "Relatedness matrix", "SAMPLE_SWAP", "</html>"):
            if needle not in html:
                failures.append("HTML missing %r" % needle)
        if "<script" in html.lower():
            failures.append("HTML should be script-free")

        for name in ("concordance_pairs_mqc.tsv", "concordance_samples_mqc.tsv"):
            p = os.path.join(out, name)
            if not os.path.exists(p):
                failures.append("missing %s" % name)
            elif not open(p).readline().startswith("# id:"):
                failures.append("%s missing MultiQC header" % name)

        # ---- exit code contract ----
        rc = subprocess.run(cmd + ["--fail-on-discordance"],
                            capture_output=True, text=True).returncode
        print("%s --fail-on-discordance exit code %d (expected 1)"
              % ("ok " if rc == 1 else "FAIL", rc))
        if rc != 1:
            failures.append("--fail-on-discordance should exit 1, got %d" % rc)

    print()
    if failures:
        print("FAILED (%d):" % len(failures))
        for f in failures:
            print("  - " + f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
