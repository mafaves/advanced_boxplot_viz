import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from advanced_boxplot_viz.utils import collect_p_values


def generate_boxplots_with_significance(
	df, group_col, biomarker_list, palette, biomarker_title_names,
	biomarker_y_axis_names, subplots_x=1, subplots_y=1, fig_size=(10, 6), 
	xtick_labels=None, image_name="boxplot.png",
	bar_height_factor=0.02, bar_tips_factor=0.005, 
	y_top_factor=0.05, y_range_factor=0.1, asterisk_factor=0.02, 
	title=True, y_labels=True, correction_method='fdr_bh', p_value_format = 'text',
	iqr_min=0.05, iqr_max=0.95, jitter_size=8, alpha=0.8, showfliers=False
):
	"""
	Generates boxplots with significance bars across multiple biomarkers.

	Parameters:
	- df (pd.DataFrame): The dataset.
	- group_col (str): The column representing groups.
	- biomarker_list (list): List of biomarker column names.
	- biomarker_y_axis_names (dict): Dictionary for axis names
	- palette (dict): Color mapping for groups.
	- biomarker_title_names (dict): Dictionary title names
	- subplots_x (int): Number of rows in subplot grid.
	- subplots_y (int): Number of columns in subplot grid.
	- fig_size (tuple): Figure size (width, height).
	- xtick_labels (list): Custom x-axis labels.
	- image_name (str): Filename for saving the figure.
	- bar_height_factor (float): Height of significance bars.
	- bar_tips_factor (float): Offset for bar tips.
	- y_top_factor (float): Space above the highest data point.
	- y_range_factor (float): Extra space above the highest data point.
	- asterisk_factor (float): Position of significance stars.
	- title (bool): Whether to show titles.
	- y_labels (bool): Whether to show y-axis labels.
	- correction_method (str): P-value correction method.
 	- p_value_format (str): 'text' or 'asterisk'.
	- iqr_min (float): Lower percentile for filtering outliers.
	- iqr_max (float): Upper percentile for filtering outliers.
	- jitter_size (int): Size of jittered points in the scatterplot.
	- alpha (float): Transparency level of jittered points.
	- showfliers (bool): Whether to show outliers in the boxplot.

	Returns:
	- None
	"""
	print("Running...")
	print("Plotting data between interquartiles 0.05 and 0.95")
	significance_dict = collect_p_values(df, group_col, biomarker_list, correction_method)

	fig, axes = plt.subplots(subplots_x, subplots_y, figsize=fig_size)
	axes = axes.flatten()

	for ax, biomarker in zip(axes, biomarker_list):
		filtered_data = df.groupby(group_col)[biomarker].apply(
			lambda x: x[(x >= x.quantile(iqr_min)) & (x <= x.quantile(iqr_max))]
		).reset_index(level=0, drop=True)
		
		# Stripplot settings
		sns.stripplot(
			data=df[df.index.isin(filtered_data.index)], x=group_col, y=biomarker, ax=ax,
			jitter=True, palette=palette, hue=group_col, legend=False, dodge=False,
			alpha=alpha, zorder=2, size=jitter_size
		)
		
		# Boxplot settings
		sns.boxplot(
			data=df[df.index.isin(filtered_data.index)], x=group_col, y=biomarker, ax=ax,
			palette=palette, hue=group_col, legend=False, showcaps=True, showfliers=showfliers,
			boxprops={'zorder': 1, 'alpha': 0.35}
		)
		ax.set_xlabel(None)

		# Subplots settings
			# Title labels
		if title:
			ax.set_title(biomarker_title_names.get(biomarker, biomarker), fontsize=20)
			# y axis labels
		if y_labels:
			ax.set_ylabel(biomarker_y_axis_names.get(biomarker, biomarker), fontsize=18)
			# x axis labels
		ax.set_xlabel(None)
		ax.set_xticks(range(len(xtick_labels)))	
		ax.set_xticklabels(xtick_labels, fontsize=16)

		# Significance bars settings
		category_positions = {category: pos for pos, category in enumerate(sorted(df[group_col].unique()))}
		
		y_range = filtered_data.max() - filtered_data.min()
		y_min, y_max = ax.get_ylim()		
		top = filtered_data.max() + (y_range * y_top_factor)
		ax.set_ylim(y_min, top + (y_range * y_range_factor))
		
		significant_combinations = significance_dict.get(biomarker, [])

		for i, (comb, p_corr) in enumerate(significance_dict.get(biomarker, [])):
			# Determine the height of the significance bar
			level = len(significant_combinations) - i
			bar_height = (y_range * bar_height_factor * level) + top
			bar_tips = bar_height - (y_range * bar_tips_factor)

			if p_value_format == "asterisk":
				sig_symbol = '***' if p_corr < 0.001 else '**' if p_corr < 0.01 else '*' if p_corr < 0.05 else 'ns'
			
			elif p_value_format == "text":
				sig_symbol = "p ≤ 0.001" if p_corr < 0.001 else "p = {:.3f}".format(p_corr)
				
			else:
				raise Exception("Sorry, invalid p-value format. It should be 'text' or 'asterik'.")
			x1, x2 = category_positions[comb[0]], category_positions[comb[1]]
			bar_height = (y_range * bar_height_factor * (i + 1)) + top
			ax.plot([x1, x1, x2, x2], [bar_tips, bar_height, bar_height, bar_tips], lw=2, c='k')
			text_height = bar_height + (y_range * asterisk_factor)
			ax.text((x1 + x2) / 2, text_height, sig_symbol, ha='center', va='bottom', c='k', fontsize=18)
	
	# Hide extra subplots
	if len(biomarker_list) < len(axes):
		for j in range(len(biomarker_list), len(axes)):
			fig.delaxes(axes[j])
	
	plt.tight_layout()
	plt.savefig(image_name, format='png', dpi=300)
	plt.show()
