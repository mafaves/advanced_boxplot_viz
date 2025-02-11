import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

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
