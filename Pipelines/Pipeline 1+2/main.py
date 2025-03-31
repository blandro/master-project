import argparse
import subprocess
import pandas as pd
import os

def run_blast(user_fasta, evalue, identity, num_threads=8):
    # Run BLASTp against reference and filter results.
    fasta_base = os.path.basename(user_fasta).rsplit('.', 1)[0]

    raw_file = f"outp/{fasta_base}_raw.txt"
    identified_proteins_file = f"outp/{fasta_base}_identified_proteins.tsv"
    transporters_reactions_file = f"outp/{fasta_base}_transporters_reactions.tsv"
    
    # Run BLASTp
    subprocess.run(
        ["blastp", "-query", user_fasta, "-db", "transporters_db", "-outfmt", "6", 
         "-evalue", str(evalue), "-num_threads", str(num_threads), "-out", raw_file],
        check=True, text=True
    )

    # Parse results
    column_names = ["QueryID", "SubjectID", "Identity", "AlignLength", 
                    "Mismatches", "GapOpens", "QStart", "QEnd", 
                    "SStart", "SEnd", "Evalue", "BitScore"]
    blast_results = pd.read_csv(raw_file, sep="\t", names=column_names, low_memory=False)

    # Filter results, and extract AID and TCID from SubjectID
    filtered_hits = blast_results.query("Evalue <= @evalue and Identity >= @identity").copy()
    filtered_hits[["AID", "TCID"]] = filtered_hits["SubjectID"].str.split("|", expand=True)
    filtered_hits.drop(columns=["SubjectID"], inplace=True)
    filtered_hits = filtered_hits[["AID", "TCID"] + [col for col in filtered_hits.columns if col not in ["AID", "TCID"]]]

    # Extract reactions from DataFrame based on matched sequence IDs
    df = pd.read_csv("transporters_df.tsv", sep="\t", dtype=str)
    results = df.merge(filtered_hits, on=["AID", "TCID"], how="inner")
    results.to_csv(identified_proteins_file, sep="\t", index=False)

    # Return reactions, relevant data and tags (Modify cols_to_use in order to selct more columns (like ChEBIs++))
    cols_to_use = [1, 2, 9, 11, 14, 15, 23]
    proteins = results.iloc[:, cols_to_use].drop_duplicates()
    proteins.rename(columns={"Reaction":"TCDB:Reaction", "R:Equation":"Rhea:Reaction"}, inplace=True)
    proteins = proteins[["UID", "TCID", "QueryID", "TCDB:Reaction", "Rhea:Reaction", "Evalue", "Identity"]]
    proteins.to_csv(transporters_reactions_file, sep="\t", index=False)

    print("Done - The BLAST is in the past!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST pipeline for transporter proteins and their reactions")
    parser.add_argument("--user_fasta", required=True, help="User-provided genome FASTA")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Identity threshold")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of CPU threads for BLAST")

    args = parser.parse_args()

    # Run BLAST and extract reactions
    run_blast(args.user_fasta, args.evalue, args.identity, args.num_threads)