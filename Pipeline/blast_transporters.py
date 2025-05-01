import argparse
import subprocess
import pandas as pd
import os
from io import StringIO


def run_blast(user_fasta, evalue, identity, num_threads=8):

    # Various path and folder works
    script_dir = os.path.dirname(os.path.realpath(__file__))
    db_path = os.path.join(script_dir, "transporters_db")
    output_dir = os.path.join(os.getcwd(), "outp")
    os.makedirs(output_dir, exist_ok=True)

    # Run BLASTp and capture stdout
    try:
        blast_process = subprocess.run(["blastp", "-query", user_fasta, "-db", db_path, "-outfmt", "6", 
     "-evalue", str(evalue), "-num_threads", str(num_threads)], check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print("BLASTp failed with error:\n", e.stderr)
        raise


    column_names = ["QueryID", "SubjectID", "Identity", "AlignLength", 
                    "Mismatches", "GapOpens", "QStart", "QEnd", 
                    "SStart", "SEnd", "Evalue", "BitScore"]

    # Parse output from stdout
    blast_results = pd.read_csv(StringIO(blast_process.stdout), sep="\t", names=column_names)

    # Filter results and extract UID + TCID
    filtered_hits = blast_results.query("Evalue <= @evalue and Identity >= @identity").copy()
    filtered_hits["QUID"] = filtered_hits["QueryID"].str.split("|", expand=True)[1]
    filtered_hits[["UID", "TCID"]] = filtered_hits["SubjectID"].str.split("|", expand=True)
    filtered_hits.drop(columns=["QueryID", "SubjectID"], inplace=True)

    # Merge with annotated transporter DF
    df_path = os.path.join(script_dir, "transporters_df.tsv")
    df = pd.read_csv(df_path, sep="\t", dtype=str)
    results = df.merge(filtered_hits, on=["UID", "TCID"], how="inner")

    cols_to_use = ["UID", "TCID", "QUID", "Reaction", "Evalue", "Identity"]
    proteins = results.loc[:, cols_to_use].drop_duplicates()

    fasta_base = os.path.splitext(os.path.basename(user_fasta))[0]
    output_file = os.path.join(output_dir, f"{fasta_base}_transporters_reactions.tsv")
    proteins.to_csv(output_file, sep="\t", index=False)

    print(f"Done - The BLAST is in the past!")


def main():
    parser = argparse.ArgumentParser(description="BLAST pipeline for transporter proteins and their reactions")
    parser.add_argument("--user_fasta", required=True, help="User-provided genome FASTA")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Identity threshold")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of CPU threads for BLAST")

    args = parser.parse_args()

    run_blast(args.user_fasta, args.evalue, args.identity, args.num_threads)


if __name__ == "__main__":
    main()