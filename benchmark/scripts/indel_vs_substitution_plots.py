import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load and process data
df = pd.read_csv('algorithm_comparison/naive_comparison.csv')

# Calculate time fractions
df['substitution_fraction'] = df['substitution_time'] / df['total_time']
df['indel_fraction'] = df['indel_time'] / df['total_time']
df['case_number'] = range(1, len(df) + 1)

# Set up the plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Create figure with subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
fig.suptitle('Naive Method: Substitutions vs Indels Time Analysis', fontsize=16, fontweight='bold')

# 1. Stacked bar chart - most effective for showing the dominance
bottom_vals = np.zeros(len(df))
ax1.bar(df['case_number'], df['substitution_fraction'] * 100, 
        label='Substitution', color='skyblue', alpha=0.8)
ax1.bar(df['case_number'], df['indel_fraction'] * 100, 
        bottom=df['substitution_fraction'] * 100,
        label='Indel', color='coral', alpha=0.8)
ax1.set_title('Time Composition by Test Case', fontweight='bold')
ax1.set_xlabel('Test Case')
ax1.set_ylabel('Time Fraction (%)')
ax1.legend()
ax1.set_ylim(0, 100)

# 2. Scatter plot showing inverse relationship
ax2.scatter(df['substitution_fraction'] * 100, df['indel_fraction'] * 100, 
           s=80, alpha=0.7, color='purple', edgecolors='black', linewidth=0.5)
ax2.set_title('Substitution vs Indel Time Fractions', fontweight='bold')
ax2.set_xlabel('Substitution Time Fraction (%)')
ax2.set_ylabel('Indel Time Fraction (%)')
ax2.grid(True, alpha=0.3)

# 3. Line plot showing trends
ax3.plot(df['case_number'], df['substitution_fraction'] * 100, 
         marker='o', linewidth=2, markersize=6, label='Substitution', color='blue')
ax3.plot(df['case_number'], df['indel_fraction'] * 100, 
         marker='s', linewidth=2, markersize=6, label='Indel', color='red')
ax3.set_title('Time Fraction Trends Across Cases', fontweight='bold')
ax3.set_xlabel('Test Case')
ax3.set_ylabel('Time Fraction (%)')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Box plot comparing distributions
data_for_box = [df['substitution_fraction'] * 100, df['indel_fraction'] * 100]
bp = ax4.boxplot(data_for_box, labels=['Substitution', 'Indel'], 
                 patch_artist=True, showmeans=True)
bp['boxes'][0].set_facecolor('lightblue')
bp['boxes'][1].set_facecolor('lightcoral')
ax4.set_title('Distribution Comparison', fontweight='bold')
ax4.set_ylabel('Time Fraction (%)')
ax4.grid(True, alpha=0.3)

# Add statistics text
stats_text = f"""Key Statistics:
Substitution: {df['substitution_fraction'].min()*100:.1f}% - {df['substitution_fraction'].max()*100:.1f}% (avg: {df['substitution_fraction'].mean()*100:.1f}%)
Indel: {df['indel_fraction'].min()*100:.1f}% - {df['indel_fraction'].max()*100:.1f}% (avg: {df['indel_fraction'].mean()*100:.1f}%)"""

fig.text(0.02, 0.02, stats_text, fontsize=10, 
         bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)
plt.savefig("algorithm_comparison/plots/temp.png", dpi=300)

# Print summary
print("Summary: Substitutions vs Indels in Naive Method")
print("=" * 50)
print(f"Substitution fraction: {df['substitution_fraction'].min()*100:.1f}% to {df['substitution_fraction'].max()*100:.1f}%")
print(f"Indel fraction: {df['indel_fraction'].min()*100:.1f}% to {df['indel_fraction'].max()*100:.1f}%")
print(f"Cases where indels > 90% of time: {(df['indel_fraction'] > 0.9).sum()}")
print(f"Cases where substitutions < 10% of time: {(df['substitution_fraction'] < 0.1).sum()}")
