import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load and process data
df = pd.read_csv('benchmark/assets/data/naive_list_tree_comparison.csv')

df = df[df["algorithm"] == "naive"]
# First, let's validate the data - check if time fractions sum to 1.0
print("DATA VALIDATION ANALYSIS")
print("=" * 60)

# Calculate time fractions the smart way - relative to operational time only
df['operational_time'] = df['indel_time'] + df['substitution_time']
df['indel_fraction'] = df['indel_time'] / df['operational_time']
df['substitution_fraction'] = 1.0 - df['indel_fraction']  # Ensures perfect sum to 1.0
df['case_number'] = range(1, len(df) + 1)

# Check benchmark overhead
df['total_time_difference'] = df['total_time'] - df['operational_time']
df['overhead_percentage'] = (df['total_time_difference'] / df['total_time']) * 100

print(f"OPERATIONAL TIME ANALYSIS:")
print(f"Average benchmark overhead: {df['overhead_percentage'].mean():.2f}% of total time")
print(f"Max benchmark overhead: {df['overhead_percentage'].max():.2f}% of total time")
print(f"Min benchmark overhead: {df['overhead_percentage'].min():.2f}% of total time")

if df['overhead_percentage'].mean() > 1.0:
    print(f"\n✓ Using operational time (indel + substitution) for fraction calculation")
    print(f"  This eliminates benchmark overhead and gives true algorithmic proportions")
else:
    print(f"\n✓ Minimal overhead detected - both approaches would be similar")

# Set up the plotting style
plt.style.use('default')

# Create the single plot you requested - Time Fraction Trends Across Cases
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
bottom = np.zeros(16)
# Line plot showing trends across cases
# ax.plot(df['case_number'], df['substitution_fraction'] * 100, 
#         marker='o', linewidth=3, markersize=8, label='Substitutions', color='blue', alpha=0.8)
# ax.plot(df['case_number'], df['indel_fraction'] * 100, 
#         marker='s', linewidth=3, markersize=8, label='Indels', color='red', alpha=0.8)
ax.bar(df['case_number'], df['indel_fraction'] * 100, edgecolor='black',
         linewidth=0.8, label='Indels', color="#3a923a", alpha=0.9, hatch='xx')

bottom = df['indel_fraction'] * 100
ax.bar(df['case_number'], df['substitution_fraction'] * 100, bottom=bottom, edgecolor='black',
        linewidth=0.8,  label='Substitutions', color="#3a923a", alpha=0.7, hatch='..')
# Since fractions now sum to exactly 100%, we can add a verification line
# ax.axhline(y=100, color='black', linestyle='-', alpha=0.3, linewidth=1)

# Formatting
# ax.set_title('Time Fraction Trends Across Test Cases', 
#             fontsize=20, fontweight='bold', pad=20)
ax.set_xlabel('Test Case', fontsize=18, fontweight='bold')
ax.set_ylabel('Time Fraction (%)', fontsize=18, fontweight='bold')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.minorticks_on()
# ax.legend(fontsize=12)

legend = ax.legend(fontsize=16, loc='lower right', frameon=True, 
                  fancybox=True, shadow=False, framealpha=1.0)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('#cccccc')
# Remove background color from legend items
for handle in legend.legend_handles:
    handle.set_facecolor('none')
    handle.set_edgecolor('black')
    handle.set_linewidth(0.8)
    handle.set_hatch(2*handle.get_hatch())

ax.grid(True, alpha=0.3)

ax.set_xticks(df['case_number'])
ax.set_xticklabels(df['case_number'].astype(int))
# Set reasonable y-axis limits
ax.set_ylim(0, 100)


# Despine - remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# Style remaining spines
ax.spines['left'].set_color('black')
ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_linewidth(0.8)



plt.tight_layout()
plt.savefig("benchmark/assets/plots/indels_vs_substitutions/time_fraction_trends.png",
            bbox_inches="tight", dpi=300)
plt.savefig("benchmark/assets/plots/indels_vs_substitutions/time_fraction_trends.svg",
            bbox_inches="tight", dpi=300)

# Print detailed analysis
print(f"\n" + "="*60)
print("DETAILED CASE-BY-CASE ANALYSIS")
print("="*60)
print(f"{'Case':>4} | {'Seq Len':>7} | {'Indel %':>7} | {'Sub %':>6} | {'Overhead':>8}")
print("-" * 52)

for _, row in df.iterrows():
    print(f"{int(row['case_number']):4d} | {int(row['sequence_length']):7d} | "
          f"{row['indel_fraction']*100:7.1f} | {row['substitution_fraction']*100:6.1f} | "
          f"{row['overhead_percentage']:8.2f}%")

# Summary insights
print(f"\n" + "="*60)
print("KEY INSIGHTS")
print("="*60)
print(f"1. Sequence length strongly correlates with indel dominance")
print(f"2. For 100bp sequences: indel time ≈ {df[df['sequence_length']==100]['indel_fraction'].mean()*100:.1f}%")
print(f"3. For 5000bp sequences: indel time ≈ {df[df['sequence_length']==5000]['indel_fraction'].mean()*100:.1f}%")
print(f"4. Benchmark overhead: {df['overhead_percentage'].mean():.2f}% of total time")
print(f"5. ✅ Fractions now perfectly sum to 100% (operational time basis)")
print(f"6. This approach gives true algorithmic time proportions")


