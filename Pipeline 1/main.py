import argparse
import subprocess
import pandas as pd

def run_blast(user_fasta, evalue, identity, num_threads=8):
    """Run BLASTp against reference and filter results."""
    
    # Run BLAST
    blast_cmd = f"blastp -query {user_fasta} -db transporters_db -outfmt 6 -evalue {evalue} -num_threads {num_threads} -out results.txt"
    subprocess.run(blast_cmd, shell=True)

    # Parse results
    blast_results = pd.read_csv("results.txt", sep="\t", header=None)
    blast_results.columns = ["QueryID", "SubjectID", "Identity", "AlignLength", 
                             "Mismatches", "GapOpens", "QStart", "QEnd", 
                             "SStart", "SEnd", "Evalue", "BitScore"]

    # Filter results
    filtered_hits = blast_results[
        (blast_results["Evalue"] <= evalue) & 
        (blast_results["Identity"] >= identity)
    ]
    filtered_hits_cp = filtered_hits.copy()
    filtered_hits_cp[["UID", "TCID"]] = filtered_hits["SubjectID"].str.split("|", expand=True)
    filtered_hits_cp = filtered_hits_cp.drop(columns=["SubjectID"])
    cols = ["UID", "TCID"] + [col for col in filtered_hits_cp.columns if col not in ["UID", "TCID"]]
    filtered_hits_cp = filtered_hits_cp[cols]
    return filtered_hits_cp

def extract_reactions(filtered_hits):
    """Extract reactions from DataFrame based on matched sequence IDs."""
    df = pd.read_csv("transporters_df.tsv", sep="\t")
    results = pd.merge(df, filtered_hits, on=["UID", "TCID"], how="inner")
    results.to_csv("identified_proteins.tsv", sep="\t", index=False)
    print("Done BLASTing!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST pipeline for transporter proteins")
    parser.add_argument("--user_fasta", required=True, help="User-provided genome FASTA")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Identity threshold")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of CPU threads for BLAST")

    args = parser.parse_args()

    # Run BLAST and extract reactions
    filtered_hits = run_blast(args.user_fasta, args.evalue, args.identity, args.num_threads)
    extract_reactions(filtered_hits)