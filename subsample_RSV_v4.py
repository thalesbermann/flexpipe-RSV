import pandas as pd
import numpy as np

metadata_file = "metadata.tsv"

accession_output = "subsample_accessions.txt"
metadata_output = "metadata_subsample.tsv"

# -----------------------------
# parâmetros
# -----------------------------
focal_country = "Brazil"

# controle por país+ano
max_per_country_year = 20
max_per_lineage_country_year = 4

# balanceamento temporal global
target_per_year_global = 200

# NOVO: teto global por país (sem Brasil)
max_per_country_total = 100

# truncamento de lineage
lineage_levels = 3

# seed para reprodutibilidade
random_seed = 42

# -----------------------------
# funções auxiliares
# -----------------------------
def truncate_lineage(lineage, levels=3):
    parts = str(lineage).strip().split(".")
    parts = [p for p in parts if p != ""]
    return ".".join(parts[:levels]) if parts else ""

def build_metadata_score(df):
    score = pd.Series(0.0, index=df.index)

    date_str = df["sampleCollectionDate"].astype(str).str.strip()

    date_is_full = date_str.str.match(r"^\d{4}-\d{2}-\d{2}$")
    date_is_year_month = date_str.str.match(r"^\d{4}-\d{2}$")

    score += date_is_year_month.astype(int) * 2.0
    score += date_is_full.astype(int) * 4.0

    if "geoLocCountry" in df.columns:
        score += (df["geoLocCountry"].astype(str).str.strip() != "").astype(int)
    if "geoLocAdmin1" in df.columns:
        score += (df["geoLocAdmin1"].astype(str).str.strip() != "").astype(int)
    if "geoLocAdmin2" in df.columns:
        score += (df["geoLocAdmin2"].astype(str).str.strip() != "").astype(int)
    if "geoLocCity" in df.columns:
        score += (df["geoLocCity"].astype(str).str.strip() != "").astype(int)

    if "completeness" in df.columns:
        score += (df["completeness"].astype(str).str.strip() != "").astype(int) * 0.5

    if "length" in df.columns:
        length_num = pd.to_numeric(df["length"], errors="coerce")
        max_length = length_num.max()
        if pd.notna(max_length) and max_length > 0:
            score += length_num.fillna(0) / max_length

    return score

def select_country_year(group, max_total, max_per_lineage, random_seed=42):
    if group.empty:
        return group.copy()

    rng = np.random.default_rng(random_seed)

    group = group.copy()
    group["tie"] = rng.random(len(group))

    group = group.sort_values(
        by=["metadata_score", "tie"],
        ascending=[False, True]
    )

    selected = []
    selected_ids = set()

    # 1. melhor por lineage
    for lin, lin_df in group.groupby("lineage_trunc"):
        top = lin_df.head(1)
        selected.append(top)
        selected_ids.update(top["accessionVersion"])

    selected = pd.concat(selected)

    # 2. completar por lineage
    for lin, lin_df in group.groupby("lineage_trunc"):
        current = selected[selected["lineage_trunc"] == lin].shape[0]
        need = max_per_lineage - current

        if need <= 0:
            continue

        candidates = lin_df[~lin_df["accessionVersion"].isin(selected_ids)]
        add = candidates.head(need)

        selected = pd.concat([selected, add])
        selected_ids.update(add["accessionVersion"])

    # 3. completar geral
    if len(selected) < max_total:
        remaining = group[~group["accessionVersion"].isin(selected_ids)]
        add = remaining.head(max_total - len(selected))
        selected = pd.concat([selected, add])

    return selected.head(max_total).drop(columns=["tie"])

# -----------------------------
# carregar metadata
# -----------------------------
df = pd.read_csv(metadata_file, sep="\t", dtype=str).fillna("")

required_cols = ["geoLocCountry", "accessionVersion", "sampleCollectionDate", "lineage"]
df = df[df[required_cols].ne("").all(axis=1)].copy()

# -----------------------------
# FILTRO DE DATA
# -----------------------------
date_str = df["sampleCollectionDate"].str.strip()

is_full = date_str.str.match(r"^\d{4}-\d{2}-\d{2}$")
is_year_month = date_str.str.match(r"^\d{4}-\d{2}$")

df = df[is_full | is_year_month].copy()

df.loc[is_year_month, "sampleCollectionDate"] = df.loc[is_year_month, "sampleCollectionDate"] + "-01"

df["year"] = df["sampleCollectionDate"].str.extract(r"^(\d{4})")[0]
df = df[df["year"].astype(int) >= 2000].copy()

# -----------------------------
# score
# -----------------------------
df["metadata_score"] = build_metadata_score(df)

# -----------------------------
# lineage truncada
# -----------------------------
df["lineage_trunc"] = df["lineage"].apply(lambda x: truncate_lineage(x, lineage_levels))
df.loc[df["lineage_trunc"] == "", "lineage_trunc"] = df["lineage"]

# -----------------------------
# deduplicação
# -----------------------------
rng = np.random.default_rng(random_seed)
df["tie"] = rng.random(len(df))

df = df.sort_values(
    by=["accessionVersion", "metadata_score", "tie"],
    ascending=[True, False, True]
).drop_duplicates("accessionVersion")

df = df.drop(columns=["tie"])

# -----------------------------
# separar
# -----------------------------
focal_df = df[df["geoLocCountry"] == focal_country].copy()
context_df = df[df["geoLocCountry"] != focal_country].copy()

selected_focal = focal_df.copy()

# -----------------------------
# subsampling país+ano
# -----------------------------
selected_context_parts = []

for (country, year), group in context_df.groupby(["geoLocCountry", "year"]):
    chosen = select_country_year(
        group,
        max_total=max_per_country_year,
        max_per_lineage=max_per_lineage_country_year,
        random_seed=random_seed
    )
    selected_context_parts.append(chosen)

selected_context = pd.concat(selected_context_parts) if selected_context_parts else context_df.iloc[0:0]

# -----------------------------
# balanceamento temporal global
# -----------------------------
balanced_parts = []

for year, group in selected_context.groupby("year"):
    if len(group) > target_per_year_global:
        group = group.sort_values(
            by="metadata_score",
            ascending=False
        ).head(target_per_year_global)
    balanced_parts.append(group)

selected_context = pd.concat(balanced_parts)

# -----------------------------
# NOVO: teto global por país (sem Brasil)
# -----------------------------
country_limited_parts = []

for country, group in selected_context.groupby("geoLocCountry"):

    if len(group) > max_per_country_total:
        group = group.sort_values(
            by="metadata_score",
            ascending=False
        ).head(max_per_country_total)

    country_limited_parts.append(group)

selected_context = pd.concat(country_limited_parts)

# -----------------------------
# juntar
# -----------------------------
selected_df = pd.concat([selected_focal, selected_context])

selected_df = selected_df.sort_values(
    by="metadata_score",
    ascending=False
).drop_duplicates("accessionVersion")

# -----------------------------
# ordenar
# -----------------------------
selected_df = selected_df.sort_values(
    by=["geoLocCountry", "year", "lineage_trunc", "metadata_score"],
    ascending=[True, True, True, False]
)

# -----------------------------
# salvar
# -----------------------------
selected_df["accessionVersion"].to_csv(accession_output, index=False, header=False)

selected_df.to_csv(metadata_output, sep="\t", index=False)

# -----------------------------
# resumo
# -----------------------------
print("Sequências finais:", len(selected_df))
print("Países:", selected_df["geoLocCountry"].nunique())
print("Por ano:")
print(selected_df["year"].value_counts().sort_index())
