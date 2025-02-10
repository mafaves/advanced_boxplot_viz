import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Dictionary for axis and title names
biomarker_names = {
    "TREM2_pgml": "CSF sTREM2 (pg/mL)",
    "Fractalkina_pgml": "CSF Fractalkine (pg/mL)",
    "YKL_ngml": "CSF YKL40 (ng/mL)",
    "S100b_pgml": "CSF S100β (pg/mL)",
    "LCR_GFAP_SIMOA": "CSF GFAP (pg/mL)",
    "GFAp pgmL_plasma baseline_SIMOA": "Plasma GFAP (pg/mL)"
}

# Default color palettes
palette_1 = {0: 'blue', 1: 'green', 2: 'red', 3: 'orange'}
palette_2 = {0: 'blue', 1: 'red'}

def collect_p_values(df, group_col, biomarker_list, correction_method='fdr_bh'):
    """
    Collects p-values for pairwise comparisons across biomarkers and applies p-value correction.

    Parameters:
    - df (pd.DataFrame): The dataset.
    - group_col (str): The column representing groups.
    - biomarker_list (list): List of biomarker column names.
    - correction_method (str): P-value correction method ('bonferroni', 'fdr_bh', etc.).

    Returns:
    - dict: A dictionary mapping biomarkers to significant comparisons.
    """
    p_values, comparisons, biomarker_map = [], [], []

    for biomarker in biomarker_list:
        groups = sorted(df[group_col].unique())
        group_combinations = list(combinations(groups, 2))

        for comb in group_combinations:
            group1 = df[df[group_col] == comb[0]][biomarker].dropna().values
            group2 = df[df[group_col] == comb[1]][biomarker].dropna().values
            _, p_value = ttest_ind(group1, group2, equal_var=True)

            p_values.append(p_value)
            comparisons.append(comb)
            biomarker_map.append(biomarker)

    # Apply multiple testing correction
    _, corrected_p_values, _, _ = multipletests(p_values, alpha=0.05, method=correction_method)

    # Store significant results
    significance_dict = {}
    for biomarker, comb, p_corr in zip(biomarker_map, comparisons, corrected_p_values):
        if p_corr < 0.05:
            significance_dict.setdefault(biomarker, []).append((comb, p_corr))

    return significance_dict


def generate_boxplots_with_significance(
    df, group_col, biomarker_list, palette,
    subplots_x=1, subplots_y=1, fig_size=(10, 6), 
    xtick_labels=None, image_name="boxplot.png",
    bar_height_factor=0.02, bar_tips_factor=0.005, 
    ytop_factor=0.05, yrange_factor=0.1, asterisk_factor=0.02, 
    title=True, biomarker_title_names, y_labels=True, biomarker_y_label_names. correction_method='fdr_bh',
    iqr_min=0.05, iqr_max=0.95, jitter_size=8, alpha=0.8, showfliers=False
):
    """
    Generates boxplots with significance bars across multiple biomarkers.

    Parameters:
    - df (pd.DataFrame): The dataset.
    - group_col (str): The column representing groups.
    - biomarker_list (list): List of biomarker column names.
    - palette (dict): Color mapping for groups.
    - subplots_x (int): Number of rows in subplot grid.
    - subplots_y (int): Number of columns in subplot grid.
    - fig_size (tuple): Figure size (width, height).
    - xtick_labels (list): Custom x-axis labels.
    - image_name (str): Filename for saving the figure.
    - bar_height_factor (float): Height of significance bars.
    - bar_tips_factor (float): Offset for bar tips.
    - ytop_factor (float): Space above the highest data point.
    - yrange_factor (float): Extra space above the highest data point.
    - asterisk_factor (float): Position of significance stars.
    - title (bool): Whether to show titles.
    - biomarker_title_names (dict): Dictionary title names
    - y_labels (bool): Whether to show y-axis labels.
    - biomarker_y_label_names (dict): Dictionary for axis names
    - correction_method (str): P-value correction method.
    - iqr_min (float): Lower percentile for filtering outliers.
    - iqr_max (float): Upper percentile for filtering outliers.
    - jitter_size (int): Size of jittered points in the scatterplot.
    - alpha (float): Transparency level of jittered points.
    - showfliers (bool): Whether to show outliers in the boxplot.

    Returns:
    - None
    """
    significance_dict = collect_p_values(df, group_col, biomarker_list, correction_method)

    fig, axes = plt.subplots(subplots_x, subplots_y, figsize=fig_size)
    axes = axes.flatten()

    for ax, biomarker in zip(axes, biomarker_list):
        filtered_data = df.groupby(group_col)[biomarker].apply(
            lambda x: x[(x >= x.quantile(iqr_min)) & (x <= x.quantile(iqr_max))]
        ).reset_index(level=0, drop=True)

        sns.stripplot(
            data=df[df.index.isin(filtered_data.index)], x=group_col, y=biomarker, ax=ax,
            jitter=True, palette=palette, hue=group_col, legend=False, dodge=False,
            alpha=alpha, zorder=2, size=jitter_size
        )
        sns.boxplot(
            data=df[df.index.isin(filtered_data.index)], x=group_col, y=biomarker, ax=ax,
            palette=palette, hue=group_col, legend=False, showcaps=True, showfliers=showfliers,
            boxprops={'zorder': 1, 'alpha': 0.35}
        )

        if title:
            ax.set_title(biomarker_title_names.get(biomarker, biomarker), fontsize=20)
        if y_labels:
            ax.set_ylabel(biomarker_y_label_names.get(biomarker, biomarker), fontsize=18)
        ax.set_xticks(range(len(xtick_labels)))
        ax.set_xticklabels(xtick_labels, fontsize=16)

        # Add significance bars
        y_range = filtered_data.max() - filtered_data.min()
        top = filtered_data.max() + (y_range * ytop_factor)

        for i, (comb, p_corr) in enumerate(significance_dict.get(biomarker, [])):
            sig_symbol = '***' if p_corr < 0.001 else '**' if p_corr < 0.01 else '*' if p_corr < 0.05 else 'ns'
            x1, x2 = sorted([comb[0], comb[1]])
            bar_height = (y_range * bar_height_factor * (i + 1)) + top
            ax.plot([x1, x1, x2, x2], [bar_height, bar_height, bar_height, bar_height], lw=2, c='k')
            ax.text((x1 + x2) / 2, bar_height + (y_range * asterisk_factor), sig_symbol, ha='center', va='bottom', c='k', fontsize=18)

    plt.tight_layout()
    plt.savefig(image_name, format='png', dpi=300)
    plt.show()
