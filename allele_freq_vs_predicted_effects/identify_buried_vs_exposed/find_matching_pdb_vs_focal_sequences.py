from Bio import SeqIO
import pandas as pd

def parse_fasta(path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}

def pct_identity(s1, s2):
    if len(s1) != len(s2):
        return 0.0
    matches = sum(a == b for a, b in zip(s1, s2))
    return matches / len(s1) * 100

def find_match(focal_seq, alphafold_seqs, threshold=99.0):
    matches = []
    for af_id, af_seq in alphafold_seqs.items():
        if len(focal_seq) != len(af_seq):
            continue
        pid = pct_identity(focal_seq, af_seq)
        if pid >= threshold:
            matches.append((af_id, pid))
    if len(matches) == 0:
        return None, None, None
    if len(matches) > 1:
        best = max(matches, key=lambda x: x[1])
        warning = f"MULTI-MATCH: {[m[0] for m in matches]}"
        return best[0], best[1], warning
    return matches[0][0], matches[0][1], None

# Load AlphaFold sequences
af_ecoli = parse_fasta("/Users/gavin/Desktop/alphafold_working/UP000000625_83333_ECOLI_v4.faa")
af_human = parse_fasta("/Users/gavin/Desktop/alphafold_working/UP000005640_9606_HUMAN_v4.faa")

# Load focal sequences
focal_ecoli = parse_fasta("/Users/gavin/Desktop/alphafold_working/focal_seqs/ecoli_translated_aligned_consensus.faa")
focal_human = parse_fasta("/Users/gavin/Desktop/alphafold_working/focal_seqs/human_gencode_w_snv_noambig.faa")

# Load RaSP-Gencode map (human only)
rasp_map = pd.read_csv("/Users/gavin/Desktop/alphafold_working/focal_seqs/rasp_gencode_map.tsv", sep="\t")
uniprot_to_ensembl = {}
for _, row in rasp_map.iterrows():
    uniprot = row["RaSP"].split(";")[0]
    ensembl = row["Gencode"].split("|")[0]
    uniprot_to_ensembl[uniprot] = ensembl

results = []

# E. coli matches
for focal_id, focal_seq in focal_ecoli.items():
    af_id, pid, multi_warn = find_match(focal_seq, af_ecoli)
    results.append({
        "species": "ecoli",
        "focal_id": focal_id,
        "af_id": af_id,
        "pct_identity": pid,
        "match_found": af_id is not None,
        "multi_match_warning": multi_warn,
        "rasp_check": None
    })

# Human matches
for focal_id, focal_seq in focal_human.items():
    af_id, pid, multi_warn = find_match(focal_seq, af_human)
    rasp_check = None
    if af_id is not None:
        uniprot = af_id.split("-")[1]
        expected_ensembl = uniprot_to_ensembl.get(uniprot)
        if expected_ensembl is None:
            rasp_check = f"WARNING: {uniprot} not in RaSP map"
        else:
            focal_id_noversion = focal_id.split(".")[0]
            expected_noversion = expected_ensembl.split(".")[0]
            if focal_id_noversion == expected_noversion:
                rasp_check = "OK"
            else:
                rasp_check = f"MISMATCH: expected {expected_ensembl}, got {focal_id}"
    results.append({
        "species": "human",
        "focal_id": focal_id,
        "af_id": af_id,
        "pct_identity": pid,
        "match_found": af_id is not None,
        "multi_match_warning": multi_warn,
        "rasp_check": rasp_check
    })

df = pd.DataFrame(results)
print(df[df["match_found"]].to_string())
print(f"\nTotal matches: {df['match_found'].sum()} / {len(df)}")
print(f"Multi-matches: {df['multi_match_warning'].notna().sum()}")
print(f"Human RaSP mismatches: {(df['rasp_check'].str.startswith('MISMATCH', na=False)).sum()}")
print(f"Human RaSP missing: {(df['rasp_check'].str.startswith('WARNING', na=False)).sum()}")

df.to_csv("/Users/gavin/Desktop/alphafold_working/matching_results.csv", index=False)