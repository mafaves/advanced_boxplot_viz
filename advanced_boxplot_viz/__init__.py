# -*- coding: utf-8 -*-
# @Time    : 11/2/25 10:20 AM
# @Author  : Marcos Aguilella
# @Affiliation  : IDIVAL
# @Email   : marcos.aguilella@idival.org
# @File    : __init__.py

# __init__.py in the main folder
from .boxplot_with_significance_bars import generate_boxplots_with_significance
from .utils import collect_p_values
from .stats import check_variance_homogeneity, check_normality

__all__ = ['boxplot_with_significance_bars', 'collect_p_values', 'check_variance_homogeneity', 'check_normality']


