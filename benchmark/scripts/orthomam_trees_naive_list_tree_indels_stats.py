
from pathlib import Path

import pandas as pd


data_path = Path("benchmark/assets/data/algorithms_comparison_indels_orthomam.csv").resolve()


df_labels = ["tree", "Root length", "Method", "Time"]

times_df = pd.read_csv(data_path)

times_df["Time"] *= 1000

times_df_naive = times_df[times_df["Method"] == "Naive"]
print(times_df_naive.groupby("Root length").describe())

times_df_blocklist = times_df[times_df["Method"] == "Block list"]
print(times_df_blocklist.groupby("Root length").describe())

times_df_blocktree = times_df[times_df["Method"] == "Block tree"]
print(times_df_blocktree.groupby("Root length").describe())
