import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load both datasets from algorithm_comparison directory
naive_df = pd.read_csv('benchmark/assets/data/naive_comparison.csv')
optimized_df = pd.read_csv('benchmark/assets/data/tree_list_comparison.csv')

# Add algorithm column to naive data and combine
naive_df['algorithm'] = 'naive'
combined_df = pd.concat([naive_df, optimized_df], ignore_index=True)

combined_df['actual_total_time'] = combined_df['substitution_time'] + combined_df['indel_time']
# Calculate fractions
combined_df['substitution_fraction'] = combined_df['substitution_time'] / combined_df['actual_total_time']
combined_df['indel_fraction'] = combined_df['indel_time'] / combined_df['actual_total_time']

# Set up the plotting style to match the reference
plt.style.use('default')  # Use default style instead of seaborn
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Colors matching the reference plot exactly
colors = {
    'naive': '#3a923a',     # Green from reference
    'list': '#3274a1',      # Blue from reference  
    'tree': '#e1812c'       # Orange from reference
}

# Calculate averages and rename algorithms
algorithms = ['naive', 'list', 'tree']
algorithm_labels = ['Naive', 'Block list', 'Block tree']
indel_avgs = [combined_df[combined_df['algorithm']==alg]['indel_fraction'].mean()*100 for alg in algorithms]
subst_avgs = [combined_df[combined_df['algorithm']==alg]['substitution_fraction'].mean()*100 for alg in algorithms]

x = np.arange(len(algorithms))
width = 0.35

# Create bars with patterns and slight hue shift
# Apply slight hue shift for better visual distinction
bars1 = ax.bar(x - width/2, indel_avgs, width, label='Indel Time', 
               color=[colors[alg] for alg in algorithms], 
               hatch='xx', edgecolor='black', linewidth=0.8, alpha=0.9)
bars2 = ax.bar(x + width/2, subst_avgs, width, label='Substitution Time',
               color=[colors[alg] for alg in algorithms], 
               hatch='..', edgecolor='black', linewidth=0.8, alpha=0.7)

# Style matching the reference
ax.set_title('Average Time Distribution by Algorithm', fontweight='bold', fontsize=14, pad=20)
ax.set_xlabel('Algorithm', fontsize=12, fontweight='bold')
ax.set_ylabel('Average Time Fraction (%)', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(algorithm_labels, fontsize=11)

# Legend styling with transparent legend item backgrounds
legend = ax.legend(fontsize=11, loc='upper right', frameon=True, 
                  fancybox=True, shadow=False, framealpha=0.8)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('#cccccc')
# Remove background color from legend items
for handle in legend.legend_handles:
    handle.set_facecolor('none')
    handle.set_edgecolor('black')
    handle.set_linewidth(0.8)
    handle.set_hatch(2*handle.get_hatch())

# Grid styling to match reference
ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)

# Add value labels on bars with styling matching reference
for bar in bars1:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 1,
             f'{height:.1f}%', ha='center', va='bottom', fontweight='bold', 
             fontsize=10, color='black')
for bar in bars2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 1,
             f'{height:.1f}%', ha='center', va='bottom', fontweight='bold', 
             fontsize=10, color='black')

# Set y-axis limit to accommodate labels
ax.set_ylim(0, max(max(indel_avgs), max(subst_avgs)) + 8)

# Tick styling to match reference
ax.tick_params(axis='both', which='major', labelsize=10, colors='black')
ax.tick_params(axis='x', which='major', length=3.5, width=0.8)
ax.tick_params(axis='y', which='major', length=3.5, width=0.8)

# Despine - remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# Style remaining spines
ax.spines['left'].set_color('black')
ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_linewidth(0.8)

# Background color
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()

# Save the plot in both SVG and PNG formats
plt.savefig("benchmark/assets/plots/indels_vs_substitutions/full_comparison.svg", bbox_inches="tight", dpi=300)
plt.savefig("benchmark/assets/plots/indels_vs_substitutions/full_comparison.png", bbox_inches="tight", dpi=300)

# Print summary
print("ALGORITHM COMPARISON SUMMARY")
print("=" * 40)
for i, alg in enumerate(algorithm_labels):
    print(f"{alg}:")
    print(f"  Indel Time: {indel_avgs[i]:.1f}%")
    print(f"  Substitution Time: {subst_avgs[i]:.1f}%")
    print()