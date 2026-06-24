from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import glob
import os

# Parse out amino acid sequences from PDB files. This is done to cross-reference with the protein sequences
# I already have (in the next step), to ensure I have 1:1 mapping of PDB to protein sequences of interest.

def extract_sequences(pdb_glob, fasta_out):
    parser = PDBParser(QUIET=True)
    with open(fasta_out, "w") as f:
        for pdb_file in glob.glob(pdb_glob):
            struct_id = os.path.basename(pdb_file).replace(".pdb", "")
            struct = parser.get_structure(struct_id, pdb_file)
            chains = list(struct[0].get_chains())
            if len(chains) != 1 or chains[0].id != "A":
                raise ValueError(f"Unexpected chains in {pdb_file}: {[c.id for c in chains]}")
            residues = [res for res in chains[0].get_residues() if res.id[0] == " "]
            seq = seq1("".join(res.resname for res in residues))
            f.write(f">{struct_id}\n{seq}\n")

#extract_sequences("/Users/gavin/Desktop/alphafold_working/UP000000625_83333_ECOLI_v4/*.pdb",
#                  "/Users/gavin/Desktop/alphafold_working/UP000000625_83333_ECOLI_v4.faa")

extract_sequences("/Users/gavin/Desktop/alphafold_working/UP000005640_9606_HUMAN_v4/*.pdb",
                  "/Users/gavin/Desktop/alphafold_working/UP000005640_9606_HUMAN_v4.faa")
