from scipy.stats import shapiro
from statsmodels.stats.multitest import multipletests
import pandas as pd

def check_normality(df, biomarker_list, group_col, correction_method):
	"""
	Checks the normality of the data for each biomarker across different groups and adjusts the p-values.
	
	Parameters:
	df (DataFrame): The DataFrame containing the data.
	biomarker_list (list): List of biomarker names.
	group_col (string): Name of the column containing the groups.
	
	Returns:
	DataFrame: A DataFrame with the adjusted p-values for the normality test for each biomarker and group.
	"""
	normality_results = {}
	all_p_values = []
	biomarker_group_pvalues = []

	# Get the different groups
	groups = df[group_col].unique()

	for biomarker in biomarker_list:
		biomarker_pvalues = []
		for group in groups:
			# Filter the DataFrame to get only the data for the current group
			group_data = df[df[group_col] == group][biomarker].dropna()

			# Perform the Shapiro-Wilk test
			stat, p_value = shapiro(group_data)
			biomarker_pvalues.append(p_value)
			all_p_values.append(p_value)
		
		biomarker_group_pvalues.append(biomarker_pvalues)
		normality_results[biomarker] = biomarker_pvalues

	# Adjust the p-values using Benjamini-Hochberg method
	_, p_values_adj, _, _ = multipletests(all_p_values, method=correction_method)

	# Create the results DataFrame
	results_df = pd.DataFrame(biomarker_group_pvalues, columns=groups, index=biomarker_list)

	# Add the adjusted p-values to the DataFrame
	p_value_index = 0
	for biomarker in biomarker_list:
		for group in groups:
			results_df.at[biomarker, group] = p_values_adj[p_value_index]
			p_value_index += 1

	return results_df

from scipy.stats import levene
from statsmodels.stats.multitest import multipletests

def check_variance_homogeneity(df, biomarker_list, group_col, correction_method):
	"""
	Checks the homogeneity of variance for each biomarker across different groups using Levene's test and adjusts the p-values.
	
	Parameters:
	df (DataFrame): The DataFrame containing the data.
	biomarker_list (list): List of biomarker names.
	group_col (string): Name of the column containing the groups.
	
	Returns:
	DataFrame: A DataFrame with the adjusted p-values for the variance homogeneity test for each biomarker.
	"""
	variance_results = {}
	all_p_values = []
	biomarker_group_pvalues = []

	# Get the different groups
	groups = df[group_col].unique()

	for biomarker in biomarker_list:
		group_data_list = []
		for group in groups:
			# Filter the DataFrame to get only the data for the current group
			group_data = df[df[group_col] == group][biomarker].dropna()
			group_data_list.append(group_data)
		
		# Perform Levene's test across all groups for the current biomarker
		stat, p_value = levene(*group_data_list)
		biomarker_group_pvalues.append(p_value)
		all_p_values.append(p_value)
		variance_results[biomarker] = p_value

	# Adjust the p-values using Benjamini-Hochberg method
	_, p_values_adj, _, _ = multipletests(all_p_values, method=correction_method)

	# Create the results DataFrame
	results_df = pd.DataFrame(index=biomarker_list, columns=['Levene p-value'])

	# Add the adjusted p-values to the DataFrame
	for i, biomarker in enumerate(biomarker_list):
		results_df.at[biomarker, 'Levene p-value'] = p_values_adj[i]

	return results_df
