import argparse
import subprocess
import pandas as pd

def run_blast(user_fasta, evalue, identity, num_threads=8):
    # Run BLASTp against reference and filter results.
    
    # Run BLASTp
    blast_cmd = f"blastp -query {user_fasta} -db transporters_db -outfmt 6 -evalue {evalue} -num_threads {num_threads} -out outp/results.txt"
    subprocess.run(blast_cmd, shell=True, check=True)

    # Parse results
    column_names = ["QueryID", "SubjectID", "Identity", "AlignLength", 
                    "Mismatches", "GapOpens", "QStart", "QEnd", 
                    "SStart", "SEnd", "Evalue", "BitScore"]
    blast_results = pd.read_csv("outp/results.txt", sep="\t", names=column_names, low_memory=False)

    # Filter results
    filtered_hits = blast_results.loc[
        (blast_results["Evalue"] <= evalue) & 
        (blast_results["Identity"] >= identity)
    ].copy()

    # Extract AID and TCID from SubjectID
    filtered_hits[["AID", "TCID"]] = filtered_hits["SubjectID"].str.split("|", expand=True)
    filtered_hits.drop(columns=["SubjectID"], inplace=True)

    # Reorder columns
    reordered_columns = ["AID", "TCID"] + [col for col in filtered_hits.columns if col not in ["AID", "TCID"]]
    filtered_hits = filtered_hits[reordered_columns]

    # Extract reactions from DataFrame based on matched sequence IDs
    df = pd.read_csv("transporters_df.tsv", sep="\t", dtype=str)

    results = pd.merge(df, filtered_hits, on=["AID", "TCID"], how="inner")
    results.to_csv("outp/identified_proteins.tsv", sep="\t", index=False)

    # Return reactions, relevant data and tags in tsv format (proteins_reactions.tsv)
    # .ipynb > test > insert here > test
    # It now works, but some minor bug fixes are necessary to improve program.
    # Some warnings are stilll dispplayed. And they stem from the section below.
    cols_to_use = [1, 2, 9, 11, 14, 15, 23]
    proteins = results.iloc[:, cols_to_use]
    proteins.drop_duplicates(inplace=True)
    proteins.rename(columns={"Reaction":"TCDB:Reaction", "R:Equation":"Rhea:Reaction"}, inplace=True)
    new_order = ["UID", "TCID", "QueryID", "TCDB:Reaction", "Rhea:Reaction", "Evalue", "Identity"]
    proteins = proteins[new_order]
    proteins.to_csv("outp/transporters_reactions.tsv", sep="\t", index=False)


    print("Done - The BLAST is in the past!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST pipeline for transporter proteins")
    parser.add_argument("--user_fasta", required=True, help="User-provided genome FASTA")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Identity threshold")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of CPU threads for BLAST")

    args = parser.parse_args()

    # Run BLAST and extract reactions
    run_blast(args.user_fasta, args.evalue, args.identity, args.num_threads)