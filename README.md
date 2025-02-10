# Boxplot Significance 📊

## Overview

This repository provides a Python tool for creating boxplots with statistical significance bars. It supports multiple p-value correction methods (Bonferroni, FDR), customizable interquartile ranges, and color palettes for better visualization of statistical comparisons.

## Features

📌 Statistical significance bars with corrected p-values

🎨 Customizable boxplots (colors, labels, and layouts)

🔬 Supports multiple p-value corrections (Bonferroni, FDR, etc.)

📊 Ideal for biomarker analysis and research visualization

## Usage
```
import pandas as pd
from boxplot_significance import generate_boxplots_with_significance

# Load example data
df = pd.read_csv("example_data.csv")

# Define parameters
group_col = "Group"
biomarker_list = ["Biomarker1", "Biomarker2"]
palette = {0: 'blue', 1: 'red'}
xtick_labels=["Control", "Disease"]

biomarker_title_names = {
    "TREM2_pgml": "CSF sTREM2 (pg/mL)",
    "Fractalkina_pgml": "CSF Fractalkine (pg/mL)",
    "YKL_ngml": "CSF YKL40 (ng/mL)",
    "S100b_pgml": "CSF S100β (pg/mL)",
    "LCR_GFAP_SIMOA": "CSF GFAP (pg/mL)",
    "GFAp pgmL_plasma baseline_SIMOA": "Plasma GFAP (pg/mL)"
}

biomarker_y_axis_names = {
    "TREM2_pgml": "CSF sTREM2 (pg/mL)",
    "Fractalkina_pgml": "CSF Fractalkine (pg/mL)",
    "YKL_ngml": "CSF YKL40 (ng/mL)",
    "S100b_pgml": "CSF S100β (pg/mL)",
    "LCR_GFAP_SIMOA": "CSF GFAP (pg/mL)",
    "GFAp pgmL_plasma baseline_SIMOA": "Plasma GFAP (pg/mL)"
}

# Generate boxplots
generate_boxplots_with_significance(df, group_col, biomarker_list, palette,
                                      subplots_x=1, subplots_y=2, fig_size=(10, 6),
                                      xtick_labels=["Control", "Disease"], image_name="plot.png",
                                      bar_height_factor=0.05, bar_tips_factor=0.01,
                                      y_top_factor=0.1, y_range_factor=0.15, asterisk_factor=0.02,
                                      title= True, biomarker_title_names, y_labels= True, biomarker_y_label_names
                                      correction_method= "fdr_bh", iqr_min, iqr_max, jitter_size, alpha, showfliers)
```

## Parameters explained

**df (DataFrame):** The input dataset containing biomarker values.

**group_col (str):** The column in df that defines groups for comparison.

**biomarker_list (list of str):** List of biomarker column names to analyze.

**palette (dict):** A dictionary mapping group labels to colors.

**subplots_x (int, default=1):** Number of rows in the subplot grid.

**subplots_y (int, default=2):** Number of columns in the subplot grid.

**fig_size (tuple, default=(10,6)):** Figure size in inches.

**xtick_labels (list, default=["Control", "Disease"]):** Labels for the x-axis.

**image_name (str, default="plot.png"):** Filename for saving the generated plot.

**bar_height_factor (float, default=0.05):** Height factor for significance bars.

**bar_tips_factor (float, default=0.01):** Reduction factor for bar tips.

**y_top_factor (float, default=0.1):** Scaling factor to adjust y-axis top margin.

**y_range_factor (float, default=0.15):** Scaling factor to adjust y-axis range.

**asterisk_factor (float, default=0.02):** Offset factor for asterisk positioning.

**title (bool, default=True)** Whether to display titles for each biomarker plot.

**biomarker_title_names (dict, optional)** Custom titles for biomarkers.

**y_labels (bool, default=True)** Whether to display y-axis labels.

**biomarker_y_label_names (dict, optional)** Custom y-axis labels for biomarkers.

**correction_method (str, default="fdr_bh")** Method for p-value correction (e.g., "bonferroni", "fdr_bh").

**iqr_min (float, optional)** Lower bound for interquartile range filtering.

**iqr_max (float, optional)** Upper bound for interquartile range filtering.

**jitter_size (float, optional)** Size of jitter points in the strip plot.

**alpha (float, optional)** Transparency level for strip plot points.

**showfliers (bool, optional)** Whether to display outliers in the boxplot.
