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
from boxplot_significance import generate_boxplots_with_significance_2

# Load example data
df = pd.read_csv("example_data.csv")

# Define parameters
group_col = "Group"
biomarker_list = ["Biomarker1", "Biomarker2"]
palette = {0: 'blue', 1: 'red'}

# Generate boxplots
generate_boxplots_with_significance_2(df, group_col, biomarker_list, palette,
                                      subplots_x=1, subplots_y=2, fig_size=(10, 6),
                                      xtick_labels=["Control", "Disease"], image_name="plot.png",
                                      barheightfactor=0.05, bartipsfactor=0.01,
                                      ytopfactor=0.1, yrangefactor=0.15, asterisk_factor=0.02)
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

**barheightfactor (float, default=0.05):** Height factor for significance bars.

**bartipsfactor (float, default=0.01):** Reduction factor for bar tips.

**ytopfactor (float, default=0.1):** Scaling factor to adjust y-axis top margin.

**yrangefactor (float, default=0.15):** Scaling factor to adjust y-axis range.

**asterisk_factor (float, default=0.02):** Offset factor for asterisk positioning.
