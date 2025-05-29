import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests

def collect_p_values(df, group_col, biomarker_list, normality_df, alpha = 0.05, correction_method='fdr_bh'):
    """
    Collects significant p-values for biomarker comparisons across groups using t-test or Mann-Whitney 
    depending on normality, and applies multiple testing correction.

    Parameters:
    - df (pd.DataFrame): Input data
    - group_col (str): Column name for group labels
    - biomarker_list (list): List of biomarker column names
    - normality_df (pd.DataFrame): DataFrame with normality test p-values (rows=biomarkers, cols=groups)
    - alpha (float): Significance level for normality and corrected p-values
    - correction_method (str): Correction method for multiple testing ('fdr_bh', 'bonferroni', etc.)

    Returns:
    - significance_dict (dict): {biomarker: [(group_pair, corrected_p)]}
    - detailed_results (list of dicts): All tests and raw p-values/statistics
    """
    raw_p_values = []
    comparisons = []
    biomarker_map = []
    results = []

    for biomarker in biomarker_list:
        groups = sorted(df[group_col].unique())
        group_combinations = list(combinations(groups, 2))

        for g1, g2 in group_combinations:
            data1 = df[df[group_col] == g1][biomarker].dropna().values
            data2 = df[df[group_col] == g2][biomarker].dropna().values

            normal1 = normality_df.loc[biomarker, g1] > alpha
            normal2 = normality_df.loc[biomarker, g2] > alpha

            if normal1 and normal2:
                stat, p = ttest_ind(data1, data2, equal_var=False)
                test = "t-test"
                print(f"t-test for {biomarker} between {g1} and {g2}: statistic={stat}, p-value={p}")
            else:
                stat, p = mannwhitneyu(data1, data2, alternative="two-sided")
                test = "Mann-Whitney"
                print(f"Mann-Whitney for {biomarker} between {g1} and {g2}: statistic={stat}, p-value={p}")

            # Store for correction
            raw_p_values.append(p)
            comparisons.append((g1, g2))
            biomarker_map.append(biomarker)

            # Detailed logging
            results.append({
                "biomarker": biomarker,
                "group1": g1,
                "group2": g2,
                "test": test,
                "statistic": stat,
                "p_value": p
            })

    # Apply multiple testing correction
    _, corrected_p_values, _, _ = multipletests(raw_p_values, alpha=alpha, method=correction_method)

    # Organize significant results
    significance_dict = {}
    for biomarker, comb, p_corr in zip(biomarker_map, comparisons, corrected_p_values):
        if p_corr < alpha:
            significance_dict.setdefault(biomarker, []).append((comb, p_corr))

    return significance_dict, results
