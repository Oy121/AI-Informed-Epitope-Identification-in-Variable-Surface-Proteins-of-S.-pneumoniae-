import pandas as pd
import os
import glob

#Filters
N_TERMINAL_DOMAINS = {"PF04650", "PS50847"}  
C_TERMINAL_DOMAINS = {"PF01473", "PF19127"} 
LPXTG_DOMAIN = "PF00746"

#Load
def load_discotope_data(base_path, protein):
    folder = os.path.join(base_path, f"{protein}_Discotope_output")
    all_csvs = glob.glob(os.path.join(folder, "*.csv"))
    if not all_csvs:
        print(f"No CSV files found for {protein}")
        return None
    df = pd.concat([pd.read_csv(f).assign(source_file=os.path.basename(f)) for f in all_csvs], ignore_index=True)
    df["protein"] = protein
    df["pdb_trimmed"] = df["pdb"].str.split(".").str[0]
    return df

#Domain Parsing

def load_interproscan_domains(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", header=None, dtype=str)
    df.columns = ["seq_id", "md5", "length", "source", "acc", "desc", "start", "end", "score", "status",
                  "date", "ipr", "ipr_desc", "go", "pathway"]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    domain_info = {}

    for seq_id, group in df.groupby("seq_id"):
        info = {
            "n_term_end": None,
            "c_term_start": None,
            "lpxtg_pos": None,
        }

        for _, row in group.iterrows():
            if row["acc"] in N_TERMINAL_DOMAINS:
                info["n_term_end"] = max(info["n_term_end"] or 0, row["end"])
            elif row["acc"] in C_TERMINAL_DOMAINS:
                info["c_term_start"] = min(info["c_term_start"] or int(row["length"]), row["start"])
            elif row["acc"] == LPXTG_DOMAIN:
                info["lpxtg_pos"] = row["end"] 

        domain_info[seq_id] = info

    return domain_info

#Filtering function
def filter_epitopes(df, domain_info):
    def keep_row(row):
        seq_id = row["pdb_trimmed"].split("_Nterminus")[0]
        try:
            res_id = int(row["res_id"])
        except ValueError:
            return False

        prot = row["protein"]
        bounds = domain_info.get(seq_id, {})
        n_cut = bounds.get("n_term_end")
        c_cut = bounds.get("c_term_start")
        lpxtg = bounds.get("lpxtg_pos")

        if n_cut and res_id <= n_cut:
            return False
        if c_cut and res_id >= c_cut:
            return False
        if lpxtg:
            lpxtg = int(lpxtg)
            if prot.startswith("zmp") and res_id < lpxtg:
                return False 
            elif prot.startswith("psp") and res_id >= lpxtg:
                return False
        return True
    
    #Change according to desired filters
    #(df["rsa"] > 0.5)
    #(df["pLDDTs"] > 50) &
    high_conf = df[
        (df["epitope"] == True) &
        (df["pLDDTs"] > 0)
    ].copy()

    high_conf["keep"] = high_conf.apply(keep_row, axis=1)
    return high_conf[high_conf["keep"] == True].drop(columns=["keep"])

# Run
def run_pipeline():
    base_discotope = "/Users/omeryurttutmus/Desktop/taxa/Discotope"
    base_tsv = "/Users/omeryurttutmus/Desktop/taxa/Initial_Data"
    output_dir = "./Filtered_Discotope_LPXTG"

    os.makedirs(output_dir, exist_ok=True)

    for protein in ["pspA", "pspC", "zmpA", "zmpB"]:
        discotope_df = load_discotope_data(base_discotope, protein)
        if discotope_df is None:
            continue

        tsv_path = os.path.join(base_tsv, f"filtered.{protein}_all_seqs.aa.tsv")
        if not os.path.exists(tsv_path):
            continue

        domain_data = load_interproscan_domains(tsv_path)
        filtered_df = filter_epitopes(discotope_df, domain_data)

        out_csv = os.path.join(output_dir, f"{protein}_filtered_epitopes_Interpro_Epitopetrue.csv")
        filtered_df.to_csv(out_csv, index=False)

if __name__ == "__main__":
    run_pipeline()
