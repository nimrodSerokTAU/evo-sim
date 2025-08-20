
from pathlib import Path

from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

DATA_PATH = Path("benchmark/assets/data/algorithms_comparison_indels_orthomam.csv").resolve()
OUTPUT_PATH = Path("benchmark/assets/plots/indel_only_comparisons").resolve()

if not DATA_PATH.exists():
    print(f"Missing data at {DATA_PATH}")
    exit(1)

df_scores = pd.read_csv(DATA_PATH)

df_scores["Time [ms]"] = df_scores["Time"]*1000


ax = sns.barplot(data=df_scores, x="Root length", y="Time [ms]", hue="Method")
ax.figure.set_size_inches(7.1, 5)
sns.despine()
hatches = ['//', '//', '//','//', 
           '..', '..', '..', '..', 
           'xx', 'xx', 'xx', 'xx',
            '///','...', 'xxx']
# Loop over the bars

for i,thisbar in enumerate(ax.patches):
    # Set a different hatch for each bar
    thisbar.set_hatch(hatches[i])

handles, old_labels = ax.get_legend_handles_labels()
# new_labels = ['Block list', 'Block tree', 'Naive']  # Your custom labels here
plt.legend(handles, old_labels, title='Methods')

plt.tight_layout()
plt.savefig(OUTPUT_PATH / "orthomam_benchmark_root_length.svg",bbox_inches="tight", dpi=300)
plt.savefig(OUTPUT_PATH / "orthomam_benchmark_root_length.png",bbox_inches="tight", dpi=300)



df_scores.groupby(["Root length", "Method"]).describe()["Time [ms]"]

df_scores


df_scores.columns.tolist()


model = ols('Time ~ C(Method) * C(Q("Root length"))', data=df_scores).fit()
anova_table = anova_lm(model, typ=2)
print(anova_table)
anova_table


tukey_algorithm = pairwise_tukeyhsd(endog=df_scores['Time'],
                                    groups=df_scores['Method'],
                                    alpha=0.05)
print(tukey_algorithm)



# For scenario main effect
tukey_scenario = pairwise_tukeyhsd(endog=df_scores['Time'],
                                   groups=df_scores['Root length'],
                                   alpha=0.05)
print(tukey_scenario)


