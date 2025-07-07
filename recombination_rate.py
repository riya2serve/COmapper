import pandas as pd

def parse_vcf(vcf_path):
    """Return dictionary: {(chrom, pos): genotype}"""
    data = {}
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            chrom, pos, ref, alt, fmt, gt_field = parts[0], int(parts[1]), parts[3], parts[4], parts[8], parts[9]
            if "GT" not in fmt: continue
            genotype = gt_field.split(":")[0]
            data[(chrom, pos)] = (ref, alt, genotype)
    return data

def estimate_recombination_rate_v2(parent1_vcf, parent2_vcf, f1_vcf):
    p1 = parse_vcf(parent1_vcf)
    p2 = parse_vcf(parent2_vcf)
    f1 = parse_vcf(f1_vcf)

    shared_sites = sorted(set(p1.keys()) & set(p2.keys()) & set(f1.keys()), key=lambda x: (x[0], x[1]))
    informative_sites = 0
    prev_origin = None
    switches = 0

    for site in shared_sites:
        p1_gt = p1[site][2]
        p2_gt = p2[site][2]
        f1_gt = f1[site][2]

        # Use homozygous parent markers only
        if p1_gt == "1/1" and p2_gt == "0/0":
            origin = "P1" if f1_gt == "1/0" or f1_gt == "0/1" else None
        elif p2_gt == "1/1" and p1_gt == "0/0":
            origin = "P2" if f1_gt == "1/0" or f1_gt == "0/1" else None
        else:
            origin = None

        if origin:
            informative_sites += 1
            if prev_origin and origin != prev_origin:
                switches += 1
            prev_origin = origin

    return switches, informative_sites
