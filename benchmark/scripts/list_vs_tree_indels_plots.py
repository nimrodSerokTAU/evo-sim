from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt

branches_df = pd.read_csv("benchmark/assets/data/list_vs_tree_indel_only.csv", index_col=0)


# branches_df = branches_df.reset_index(drop=True)
OUTPUT_DIR = Path("benchmark/assets/plots/indel_only_comparisons").resolve()


styles=['C0*-', 'C1D-']
for group in branches_df.groupby("branch_scale"):
    group[1][["blocklist_time", "blocktree_time", "insertion_rate"]].plot(x ="insertion_rate",
                                                                          style=styles, markersize=9)
    plt.title(f"Multiplication factor: {group[0]}", size=18)
    plt.legend(labels=["Block list", "Block tree"])
    
    plt.gca().set_xlabel("Insertion rate", size=16)
    plt.gca().set_ylabel("Runtime (seconds)", size=16)

    plt.savefig(OUTPUT_DIR / f"blocks_runtime_factor_{group[0]}.svg" ,bbox_inches="tight", dpi=300)
    plt.savefig(OUTPUT_DIR / f"blocks_runtime_factor_{group[0]}.png" ,bbox_inches="tight", dpi=300)


