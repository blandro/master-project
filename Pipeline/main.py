import argparse
import subprocess
import pandas as pd

def run_blast(user_fasta, reference_db, evalue, identity):
    """Run BLASTp against reference and filter results."""
    
    # Run BLAST
    blast_cmd = f"blastp -query {user_fasta} -db {reference_db} -outfmt 6 -evalue {evalue} -out results.txt"
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
    filtered_hits[["UID", "TCID"]] = filtered_hits["QueryID"].str.split("|", expand=True)
    filtered_hits = filtered_hits.drop(columns=['QueryID'])
    cols = ['UID', 'TCID'] + [col for col in filtered_hits.columns if col not in ['UID', 'TCID']]
    filtered_hits = filtered_hits[cols]
    return filtered_hits

def extract_reactions(filtered_hits, df_path):
    """Extract reactions from DataFrame based on matched sequence IDs."""
    df = pd.read_csv(df_path, sep="\t")
    results = pd.merge(df, filtered_hits, on=["UID", "TCID"], how="inner")
    results.to_csv("Identified_reactions.tsv", sep="\t", index=False)
    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST pipeline for transporter proteins")
    parser.add_argument("--user_fasta", required=True, help="User-provided genome FASTA")
    parser.add_argument("--reference_db", required=True, help="Reference protein database")
    parser.add_argument("--df_path", required=True, help="Path to transport reaction DataFrame TSV")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Identity threshold")

    args = parser.parse_args()

    # Run BLAST and extract reactions
    filtered_hits = run_blast(args.user_fasta, args.reference_db, args.evalue, args.identity)
    extract_reactions(filtered_hits, args.df_path)
