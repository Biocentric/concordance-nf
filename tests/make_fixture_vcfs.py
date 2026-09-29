#!/usr/bin/env python3
"""
Generate a synthetic cohort with a known identity structure, for end-to-end
testing against the real somalier binary.

    python3 tests/make_fixture_vcfs.py /path/to/outdir

Produces a reference FASTA, a somalier sites VCF, one single-sample VCF per
sample, and a samplesheet whose declared structure the pipeline should
reproduce exactly:

    s_a / s_b   same subject, same genotypes      -> CONFIRMED_MATCH
    s_c / s_e   different subjects, same genotypes -> UNEXPECTED_DUPLICATE
    s_f / s_g   same subject, different genotypes  -> SAMPLE_SWAP
    everything else                                -> CONFIRMED_DISTINCT
"""
import gzip
import os
import random
import subprocess
import sys

N_SITES = 4000
# 20 bp apart purely to keep the committed reference small - somalier treats
# sites as independent, so spacing has no effect on any statistic here
SPACING = 20
CHROM = "chr1"
DEPTH = 30
# per-read sequencing error. Realistic Q30 WGS is ~0.1%; anything near 1%
# pushes hom sites into somalier's "middling allele balance" band, which
# inflates het counts, depresses IBS0 and lifts the unrelated relatedness
# baseline off zero. It also drives CHARR, which estimates contamination from
# exactly these minor-allele reads.
READ_ERROR = 0.001

# sample -> (subject, genotype-seed, sex)
# samples sharing a seed get identical genotypes
COHORT = [
    ("s_a", "P1", 101, "female"),
    ("s_b", "P1", 101, "female"),   # true duplicate of s_a
    ("s_c", "P2", 202, "male"),
    ("s_d", "P3", 303, "male"),
    ("s_e", "P4", 202, "male"),     # same genotypes as s_c, declared as P4
    ("s_f", "P5", 404, "female"),   # s_f and s_g declared as one subject...
    ("s_g", "P5", 505, "female"),   # ...but drawn independently
]


def genotypes(seed, error_rate=0.002):
    """Hardy-Weinberg draws at MAF 0.5, with a little genotyping error."""
    rng = random.Random(seed)
    err = random.Random(seed * 7919 + 13)
    out = []
    for _ in range(N_SITES):
        a1 = rng.random() < 0.5
        a2 = rng.random() < 0.5
        alt = int(a1) + int(a2)
        if err.random() < error_rate:
            alt = err.choice([g for g in (0, 1, 2) if g != alt])
        out.append(alt)
    return out


def write_fasta(path):
    length = N_SITES * SPACING + SPACING
    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    with open(path, "w") as fh:
        fh.write(">%s\n" % CHROM)
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    with open(path + ".fai", "w") as fh:
        fh.write("%s\t%d\t%d\t60\t61\n" % (CHROM, length, len(CHROM) + 2))
    return seq


def positions():
    return [(i + 1) * SPACING for i in range(N_SITES)]


def write_sites(path, seq):
    lines = ["##fileformat=VCFv4.2",
             '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
             "##contig=<ID=%s,length=%d>" % (CHROM, len(seq)),
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for i, pos in enumerate(positions()):
        ref = seq[pos - 1]
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
        lines.append("%s\t%d\ts%d\t%s\t%s\t100\tPASS\tAF=0.5" % (CHROM, pos, i, ref, alt))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_sample_vcf(path, sample, gts, seq):
    rng = random.Random(hash(sample) % 100000)
    lines = ["##fileformat=VCFv4.2",
             "##contig=<ID=%s,length=%d>" % (CHROM, len(seq)),
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
             '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample]
    for i, pos in enumerate(positions()):
        ref = seq[pos - 1]
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
        n_alt = gts[i]
        gt = {0: "0/0", 1: "0/1", 2: "1/1"}[n_alt]
        dp = max(8, int(rng.gauss(DEPTH, 4)))
        if n_alt == 0:
            ad_alt = sum(1 for _ in range(dp) if rng.random() < READ_ERROR)
        elif n_alt == 2:
            ad_alt = dp - sum(1 for _ in range(dp) if rng.random() < READ_ERROR)
        else:
            ad_alt = sum(1 for _ in range(dp) if rng.random() < 0.5)
        ad_ref = dp - ad_alt
        lines.append("%s\t%d\ts%d\t%s\t%s\t100\tPASS\t.\tGT:AD:DP\t%s:%d,%d:%d"
                     % (CHROM, pos, i, ref, alt, gt, ad_ref, ad_alt, dp))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def bgzip_and_index(path):
    """Prefer real bgzip/tabix; fall back to gzip when htslib is absent."""
    try:
        subprocess.run(["bgzip", "-f", path], check=True, capture_output=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", path + ".gz"], check=True,
                       capture_output=True)
        return path + ".gz"
    except (FileNotFoundError, subprocess.CalledProcessError):
        with open(path, "rb") as src, gzip.open(path + ".gz", "wb") as dst:
            dst.write(src.read())
        os.remove(path)
        return path + ".gz"


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "fixture"
    os.makedirs(outdir, exist_ok=True)

    fasta = os.path.join(outdir, "ref.fasta")
    seq = write_fasta(fasta)

    sites = os.path.join(outdir, "sites.vcf")
    write_sites(sites, seq)
    sites_gz = bgzip_and_index(sites)

    rows = []
    for sample, subject, seed, sex in COHORT:
        gts = genotypes(seed)
        vcf = os.path.join(outdir, sample + ".vcf")
        write_sample_vcf(vcf, sample, gts, seq)
        vcf_gz = bgzip_and_index(vcf)
        rows.append((sample, subject, os.path.abspath(vcf_gz), sex))

    sheet = os.path.join(outdir, "samplesheet.csv")
    with open(sheet, "w") as fh:
        fh.write("sample,subject,family,alignment,index,sex,batch\n")
        for sample, subject, path, sex in rows:
            fh.write("%s,%s,%s,%s,,%s,fixture\n" % (sample, subject, subject, path, sex))

    print("fixture written to", os.path.abspath(outdir))
    print("  reference   ", fasta)
    print("  sites       ", sites_gz)
    print("  samplesheet ", sheet)
    print("  samples     ", ", ".join(r[0] for r in rows))


if __name__ == "__main__":
    main()
