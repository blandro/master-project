import os
import json
import pandas as pd
import ast
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import re
from collections import defaultdict
from matplotlib.lines import Line2D

script_dir = os.path.dirname(os.path.abspath(__file__))

def resolve_path(*parts, mkdir=False):
    path = os.path.join(script_dir, *parts)
    if mkdir:
        os.makedirs(os.path.dirname(path) if os.path.splitext(path)[1] else path, exist_ok=True)
    return path

tdb_path = resolve_path("transporters_df.tsv")
json_path = resolve_path("iML1515.json")
model_path = resolve_path("iML1515_transport_rxs.tsv")
c2b_path = resolve_path("chebi2bigg.tsv")

outp_dir = resolve_path("outp", mkdir=True)
figs_dir = resolve_path("figs", mkdir=True)
data_dir = resolve_path("data", mkdir=True)

tdb = pd.read_csv(tdb_path, sep="\t")

with open(json_path, "r") as f:
    iml1515 = json.load(f)
gene_uniprot_map = {g["id"]: g.get("annotation", {}).get("uniprot", ["Unknown"]) for g in iml1515["genes"]}

model = pd.read_csv(model_path, sep="\t")
model_transport_uids = set(model["UIDs"].str.split(", ").explode())
model_transport_uids.remove("Unknown")

c2b_df = pd.read_csv(c2b_path, sep="\t")
c2b = {row["CHEBI"].upper().replace("CHEBI:", "CHEBI:"): row["BIGG"] for _, row in c2b_df.iterrows()}

PR_points = []

PR_points_genes = []

def chunked_lines(items, chunk_size=10):
    return '\n'.join([', '.join(items[i:i+chunk_size]) for i in range(0, len(items), chunk_size)])


def save_venn3_plot(set1, set2, set3, labels, title, filepath):
    plt.figure(figsize=(10, 6))
    out = venn3([set1, set2, set3], set_labels=labels)
    for t in out.set_labels:
        t.set_fontsize(14)
    if out.subset_labels:
        for label in out.subset_labels:
            if label:
                label.set_fontsize(14)
    plt.title(title, fontsize=20)
    plt.tight_layout()
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()

def save_venn2_plot(set1, set2, labels, title, filepath):
    plt.figure(figsize=(6, 6))
    out = venn2([set1, set2], set_labels=labels)
    for t in out.set_labels:
        t.set_fontsize(14)
    for t in out.subset_labels:
        t.set_fontsize(14)
    plt.title(title, fontsize=18)
    plt.tight_layout()
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()


for fname in os.listdir(outp_dir):
    e_val, identity, *_ = fname.replace("_e_coli_K12_MG1655_transporters_reactions.tsv", "").split("_")
    e_val_latex = e_val.replace("e", r"$^{") + "}$"

    e_coli_path = os.path.join(outp_dir, fname)
    e_coli = pd.read_csv(e_coli_path, sep="\t")

    # Only keeping the QUID entry with the lowest E-val, as others are homologs and orthologs hits
    e_coli = e_coli.loc[e_coli.groupby("QUID")["Evalue"].idxmin()].reset_index(drop=True)
    e_coli_uids = set(e_coli["QUID"])

    all_uids_model = {uid
        for uid_list in gene_uniprot_map.values()
        for uid in uid_list
        if uid != "Unknown"
    }

    save_venn3_plot(e_coli_uids, all_uids_model, model_transport_uids,
        ("Transporter genes - BLASTp", "All genes - iML1515", "Transporter genes - iML1515"),
        f"E-value: {e_val_latex}\nIdentity: {identity}%",
        resolve_path("figs", f"transport_genes_tdb_iML1515_{e_val}_{identity}.pdf"))
    
    TP_genes = len(e_coli_uids & model_transport_uids)
    FN_genes = len(model_transport_uids - e_coli_uids)
    FP_genes = len(e_coli_uids - model_transport_uids)

    P_genes = TP_genes / (TP_genes + FP_genes) * 100
    R_genes = TP_genes / (TP_genes + FN_genes) * 100

    PR_points_genes.append({
        "Evalue": e_val,
        "Identity": identity,
        "Precision": P_genes,
        "Recall": R_genes,
    })

    set_e_coli_bigg = set()

    for reaction_list in e_coli["Reaction"]:
        reaction_list = ast.literal_eval(reaction_list)
        
        for _, chebi_str in reaction_list:
            if chebi_str != "no_reaction_identified":

                matches = re.findall(r'CHEBI:\d+', chebi_str)
                if not matches:
                    continue
                
                bigg_ids = {c2b.get(chebi.upper()) for chebi in matches}
                bigg_ids.discard(None)
                
                if bigg_ids:
                    set_e_coli_bigg.add(frozenset(bigg_ids))


    unique_model_sets = {
        frozenset(ast.literal_eval(s)) 
        for s in model["Substrates (BiGG)"] 
        if pd.notnull(s)
    }

    reactions_bigg_model = [sorted(list(s)) for s in unique_model_sets]
    set_model_bigg = set(frozenset(r) for r in reactions_bigg_model)

    save_venn2_plot(set_e_coli_bigg, set_model_bigg,
                    ("BLASTp reactions", "Model reactions"),
                    f"E-value: {e_val_latex}\nIdentity: {identity}%",
                    resolve_path("figs", f"reactions_tdb_iML1515_{e_val}_{identity}.pdf"))
    

    TP = len(set_e_coli_bigg & set_model_bigg)
    FP = len(set_e_coli_bigg) - TP
    FN = len(set_model_bigg) - TP

    precision = TP / (TP + FP) * 100
    recall = TP / (TP + FN) * 100

    PR_points.append({
        "Evalue": e_val,
        "Identity": identity,
        "Precision": precision,
        "Recall": recall,
    })

    model["Substrates (BiGG) frozen"] = model["Substrates (BiGG)"].apply(
        lambda x: frozenset(ast.literal_eval(x)) if isinstance(x, str) else frozenset(x))

    only_model_reactions = set_model_bigg - set_e_coli_bigg
    uids_by_set = {}

    for substrate_set in only_model_reactions:
        matching_rows = model[model["Substrates (BiGG) frozen"] == substrate_set]
        uid_lists = matching_rows["UIDs"].tolist()
        uids = set()
        for uid_entry in uid_lists:
            if isinstance(uid_entry, str):
                split_uids = [uid.strip() for uid in uid_entry.split(",")]
                uids.update(split_uids)
            elif isinstance(uid_entry, list):
                uids.update(uid_entry)

        uids_by_set[substrate_set] = list(uids)

######## Only used for misc. identification of reactions ########

    reactions_by_set = defaultdict(dict)

    for substrate_set, uids in uids_by_set.items():
        for uid in uids:
            match = tdb[tdb["UID"] == uid]
            if match.empty:
                reactions_by_set[substrate_set][uid] = "Gene not in TDB"
            else:
                rxns = match["Reaction"].dropna().tolist()
                parsed_rxns = []
                for rx_group in rxns:
                    try:
                        parsed_rxns.extend(ast.literal_eval(rx_group))
                    except Exception:
                        parsed_rxns.append((rx_group, ""))
                reactions_by_set[substrate_set][uid] = parsed_rxns

    file = os.path.join(data_dir, f"only_reactions_in_model_comp_tdb_{e_val}_{identity}.txt")

    with open(file, "w") as f:
        f.write("This is all reaction sets in the model,\n"
                "that do not exist in the BLAST results.\n"
                "A mapping of the corresponding genes was done, and then \n"
                "mapped over to the reactions these genes facilitate\n"
                "according to TDB.\n"
                f"E-value: {e_val}\n"
                f"Identity: {identity}\n\n")
        
        for i, (substrate_set, uid_data) in enumerate(reactions_by_set.items(), 1):
            f.write(f"--- Reaction Set {i} ---\n")
            f.write(f"Substrates (BiGG IDs): {', '.join(sorted(substrate_set))}\n\n")
            f.write("Genes facilitating this reaction according to iML1515:\n")
            f.write(f"{', '.join(uid_data.keys())}\n\n")

            model_matches = model[model["Substrates (BiGG) frozen"] == substrate_set]
            reaction_rows = model_matches[["Reaction (Names)", "Reaction (CHEBI)"]].dropna()

            f.write("Reactions from iML1515:\n")
            for _, row in reaction_rows.iterrows():
                name_rxn = row["Reaction (Names)"]
                chebi_rxn = row["Reaction (CHEBI)"]
                f.write(f"  - {name_rxn}  |  {chebi_rxn}\n")
                
            f.write("\nReactions from TDB:\n")
            for uid, rxns in uid_data.items():
                f.write(f"Gene: {uid}\n")
                if isinstance(rxns, str):
                    if uid != "Unknown":
                        f.write(f"  - {rxns}\n")
                else:
                    for rxn in rxns:
                        f.write(f"  - {rxn[0]}  |  {rxn[1]}\n")
                f.write("\n")
#################################################################

def plot_precision_recall(
    df,
    output_filename,
    marker_size=150,
    font_size=14,
    xlabel="Recall (%)",
    ylabel="Precision (%)"
):
    df = pd.DataFrame(df)
    df["Identity"] = df["Identity"].astype(int)
    
    unique_evalues = sorted(df["Evalue"].unique())
    unique_identities = sorted(df["Identity"].unique())

    marker_styles = ["o", "s", "D", "^", "v", ">", "<", "P", "*", "X"]
    evalue_markers = {
        e: marker_styles[i % len(marker_styles)]
        for i, e in enumerate(unique_evalues)
    }
    identity_colors = {
        id_: plt.cm.viridis(i / (len(unique_identities)-1))
        for i, id_ in enumerate(unique_identities)
    }

    plt.figure(figsize=(10, 7))
    for _, row in df.iterrows():
        plt.scatter(
            row["Recall"], row["Precision"],
            marker=evalue_markers[row["Evalue"]],
            color=identity_colors[row["Identity"]],
            s=marker_size,
            edgecolor="black"
        )

    legend_elements = []

    legend_elements.append(Line2D([0], [0], linestyle='none', label=r'$\bf{E\text{-}value}$', color='white'))
    legend_elements += [
        Line2D([0], [0], marker=marker, linestyle='None', color='w',
               label=f"{e}", markerfacecolor='gray', markersize=12, markeredgecolor='black')
        for e, marker in evalue_markers.items()
    ]
    legend_elements.append(Line2D([0], [0], linestyle='none', label=r'$\bf{Identity}$', color='white'))
    legend_elements += [
        Line2D([0], [0], marker='o', linestyle='None', color='w',
               label=f"{id_}%", markerfacecolor=color, markersize=12, markeredgecolor='black')
        for id_, color in identity_colors.items()
    ]

    plt.legend(
        handles=legend_elements,
        bbox_to_anchor=(1, 1.04),
        loc="upper left",
        fontsize=font_size,
        borderaxespad=1.2,
        handlelength=2,
        handletextpad=0.8,
        labelspacing=0.8,
        borderpad=1.0
    )

    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(resolve_path("figs", output_filename))
    plt.close()

plot_precision_recall(PR_points, "precision_recall_reactions.pdf")
plot_precision_recall(PR_points_genes, "precision_recall_genes.pdf")