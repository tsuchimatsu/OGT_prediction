import pandas as pd
import warnings
from Bio import SeqIO
import os
import re
import numpy as np

warnings.simplefilter("ignore")


def calc_genome(df, genome_path):
    # calculate genome size
    genome_size = sum(len(record.seq) for record in SeqIO.parse(genome_path, "fasta"))
    df["genome_size"] = genome_size
    return df


def calc_rrna(df, rrna_path):
    # calculate base fraction of rRNA
    target = ["A", "T", "G", "C"]
    rrna_types = {"5S": "5S_rRNA", "16S": "16S_rRNA", "23S": "23S_rRNA"}

    for rrna, prefix in rrna_types.items():
        for nuc in target:
            df[f"{prefix}_{nuc}"] = 0

    total_lengths = {rrna: 0 for rrna in rrna_types}

    for record in SeqIO.parse(rrna_path, "fasta"):
        for rrna, prefix in rrna_types.items():
            if record.description.startswith(rrna):
                seq = record.seq
                for nuc in target:
                    df[f"{prefix}_{nuc}"] += seq.count(nuc)
                total_lengths[rrna] += len(seq)

    for rrna, prefix in rrna_types.items():
        if total_lengths[rrna] > 0:
            for nuc in target:
                df[f"{prefix}_{nuc}"] /= total_lengths[rrna]

    return df


def calc_rrna_mfe(df, rrna_li_path):
    # calculate MFE of rRNA
    rrna_types = {"5": "5S", "16": "16S", "23": "23S"}
    mfe_counts = {key: 0 for key in rrna_types.values()}
    mfe_totals = {key: 0.0 for key in rrna_types.values()}

    for rrna in rrna_types.values():
        df[f"{rrna}_mfe_len"] = 0

    for record in SeqIO.parse(rrna_li_path, "fasta"):
        desc_part = record.description
        seq = record.seq

        p = next((i for i in range(100) if seq[-i] == "("), -1)
        if p == -1:
            continue
        mfe = float(str(seq[-p + 1 : -2]))

        q = p + 1
        while True:
            if seq[-q] != "(" and seq[-q] != "." and seq[-q] != ")":
                break
            q += 1
        structure_len = len(str(seq[-q + 1 : -p]))

        for key, rrna in rrna_types.items():
            if desc_part.startswith(key):
                mfe_totals[rrna] += mfe / structure_len
                mfe_counts[rrna] += 1
                break

    for rrna in rrna_types.values():
        if mfe_counts[rrna] > 0:
            df[f"{rrna}_mfe_len"] = mfe_totals[rrna] / mfe_counts[rrna]

    return df


def calc_trna(df, trna_path):
    # calculate base fraction of tRNA
    target = ["A", "T", "G", "C"]

    for nuc in target:
        df[f"tRNA_{nuc}"] = 0

    total_length = 0

    for record in SeqIO.parse(trna_path, "fasta"):
        if "Undet" in record.description.split(" ")[3]:
            continue

        seq = record.seq
        for nuc in target:
            df[f"tRNA_{nuc}"] += seq.count(nuc)

        total_length += len(seq)

    if total_length > 0:
        df[[f"tRNA_{nuc}" for nuc in target]] /= total_length

    return df


def calc_trna_mfe(df, trna_li_path):
    # calculate MFE of tRNA
    df["trna_mfe_len"] = 0
    num = 0

    for record in SeqIO.parse(trna_li_path, "fasta"):
        seq = record.seq

        p = next((i for i in range(100) if seq[-i] == "("), -1)
        if p == -1:
            continue

        mfe = float(str(seq[-p + 1 : -2]))

        q = p + 1
        while True:
            if seq[-q] != "(" and seq[-q] != "." and seq[-q] != ")":
                break
            q += 1
        structure_len = len(str(seq[-q + 1 : -p]))

        if structure_len > 0:
            df["trna_mfe_len"] += mfe / structure_len

        num += 1

    if num > 0:
        df["trna_mfe_len"] /= num

    return df


def calc_protein(df, cds_faa_path):
    # calculate amino acids fraction
    target = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]
    target2 = []
    for amino1 in target:
        for amino2 in target:
            target2.append(amino1 + amino2)

    for aa in target + target2:
        df["Pro_" + aa] = 0

    len_pro = 0
    count_pro = 0

    for record in SeqIO.parse(cds_faa_path, "fasta"):
        if record.description.split(";")[1][-2:] != "00":
            continue

        seq = record.seq
        len_pro += len(seq)
        count_pro += 1

        for aa in target:
            df[f"Pro_{aa}"] += seq.count(aa)
        for aa in target2:
            df["Pro_" + aa] += seq.count_overlap(aa)

    cols = [f"Pro_{x}" for x in target]
    df.loc[:, cols] /= df[cols].sum(axis=1).values[:, None]

    cols2 = [f"Pro_{x}" for x in target2]
    df.loc[:, cols2] /= df[cols2].sum(axis=1).values[:, None]

    df["Pro_mean_length"] = len_pro / count_pro

    Pro_polar_uncharged = sum(df[f"Pro_{x}"] for x in ["S", "T", "N", "Q"])
    Pro_polar_charged = sum(df[f"Pro_{x}"] for x in ["D", "E", "K", "R", "H"])

    df["Pro_polar_charged"] = Pro_polar_uncharged / Pro_polar_charged
    Pro_polar_hydrophobic = sum(
        df[f"Pro_{x}"] for x in ["A", "V", "I", "L", "M", "F", "Y", "W"]
    )

    df["Pro_polar_hydrophobic"] = Pro_polar_uncharged / Pro_polar_hydrophobic
    df["Pro_LK/Q"] = (df["Pro_L"] + df["Pro_K"]) / df["Pro_Q"]
    df["Pro_EK/QH"] = (df["Pro_E"] + df["Pro_K"]) / (df["Pro_Q"] + df["Pro_H"])

    return df


def calc_codon(df, cds_fna_path):
    # calculate codon fraction
    DNA2Protein = {
        "TTT": "F",
        "TCT": "S",
        "TAT": "Y",
        "TGT": "C",
        "TTC": "F",
        "TCC": "S",
        "TAC": "Y",
        "TGC": "C",
        "TTA": "L",
        "TCA": "S",
        "TTG": "L",
        "TCG": "S",
        "TGG": "W",
        "CTT": "L",
        "CCT": "P",
        "CAT": "H",
        "CGT": "R",
        "CTC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CTG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "ATT": "I",
        "ACT": "T",
        "AAT": "N",
        "AGT": "S",
        "ATC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "ATA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "ATG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GTT": "V",
        "GCT": "A",
        "GAT": "D",
        "GGT": "G",
        "GTC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GTA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GTG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }

    target1 = ["A", "T", "G", "C"]
    target2 = [
        "AA",
        "AT",
        "AG",
        "AC",
        "TA",
        "TT",
        "TG",
        "TC",
        "GA",
        "GT",
        "GG",
        "GC",
        "CA",
        "CT",
        "CG",
        "CC",
    ]
    target3 = []

    for nuc1 in target1:
        for nuc2 in target2:
            target3.append(nuc1 + nuc2)

    dict = {}

    for nuc in target3:
        dict["ORF_" + nuc] = 0

    for record in SeqIO.parse(cds_fna_path, "fasta"):
        desc_part = record.description
        if desc_part.split(";")[1][-2:] != "00":
            continue
        seq = record.seq

        codon_l = re.split("(...)", str(seq))[1::2]
        for codon in codon_l:
            if "N" in codon:
                continue
            dict["ORF_" + codon] += 1

    target3.remove("TAA")
    target3.remove("TAG")
    target3.remove("TGA")
    target3.remove("ATG")
    target3.remove("TGG")

    for codon in target3:
        same_codon = []
        for c in target3:
            if DNA2Protein[codon] == DNA2Protein[c]:
                same_codon.append(c)
        sum = 0
        for c in same_codon:
            sum += dict["ORF_" + c]
        df[DNA2Protein[codon] + "_ORF_" + codon] = dict["ORF_" + codon] / sum

    return df


def calc_fraction(df, cds_fna_path, cds_faa_path):
    # calculate of fraction around start codon

    min_length = 600  # remove cds less than 600 nt
    calc_length = 60

    # guanine fraction
    target = ["A", "T", "G", "C"]
    d = {}
    index = 0
    for nuc in target:
        d[nuc] = index
        index += 1

    count = np.zeros((len(target), calc_length))

    for record in SeqIO.parse(cds_fna_path, "fasta"):
        id_part = record.id
        desc_part = record.description
        if desc_part.split(";")[1][-2:] != "00":
            continue
        seq = record.seq

        if len(seq) < min_length:
            continue

        for i in range(calc_length):
            if seq[i] in target:
                count[d[seq[i]], i] += 1

    count /= count.sum(axis=0)

    nuc = "G"
    window_size = 9
    l = count[d[nuc]]
    l = np.convolve(count[d[nuc]], np.ones(window_size) / window_size, mode="valid")
    df["G_fraction_normalized"] = np.log2(l[5] / l[49])

    # amino acids
    min_length = 200
    calc_length = 20
    target = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]

    d = {}
    index = 0
    for aa in target:
        d[aa] = index
        index += 1

    count = np.zeros((len(target), calc_length))

    for record in SeqIO.parse(cds_faa_path, "fasta"):
        id_part = record.id
        desc_part = record.description
        if desc_part.split(";")[1][-2:] != "00":
            continue
        seq = record.seq

        if len(seq) < min_length:
            continue

        for i in range(calc_length):
            if seq[i] in target:
                count[d[seq[i]], i] += 1

    count /= count.sum(axis=0)

    for aa in target:
        df["Pro_" + aa + "_frac_normed"] = np.log2(
            np.max([count[d[aa]][1:5].mean(), 0.00001]) / df["Pro_" + aa] # avoid 0
        )

    return df


def main():
    org = "GB_GCA_000814675.1"  # accession
    data_dir = "../data/genome_example"
    genome_path = os.path.join(data_dir, "genome", f"{org}_genomic.fna")
    rrna_path = os.path.join(data_dir, "rrna", f"{org}.fa")
    rrna_li_path = os.path.join(data_dir, "rrna", f"{org}_rrna_linearfold.fa")
    trna_path = os.path.join(data_dir, "trna", f"{org}.fa")
    trna_li_path = os.path.join(data_dir, "trna", f"{org}_trna_linearfold.fa")
    cds_faa_path = os.path.join(data_dir, "cds", f"{org}_protein.faa")
    cds_fna_path = os.path.join(data_dir, "cds", f"{org}_protein.fna")
    columns_file = os.path.join(data_dir, "calculated_features_example", "columns.txt")
    output_file = os.path.join(
        data_dir, "calculated_features_example", f"{org}_features.csv"
    )

    with open(columns_file, "r") as f:
        column_names = [line.strip() for line in f]

    # make dataframe
    df = pd.DataFrame(index=[org], columns=column_names)

    # calculate features
    df = calc_genome(df, genome_path)
    df = calc_rrna(df, rrna_path)
    df = calc_rrna_mfe(df, rrna_li_path)
    df = calc_trna(df, trna_path)
    df = calc_trna_mfe(df, trna_li_path)
    df = calc_protein(df, cds_faa_path)
    df = calc_codon(df, cds_fna_path)
    df = calc_fraction(df, cds_fna_path, cds_faa_path)

    # output
    df.iloc[0, :].to_csv(output_file)
    print(f"Features saved to {output_file}")


if __name__ == "__main__":
    main()
