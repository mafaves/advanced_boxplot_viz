from setuptools import setup, find_packages

# if __name__ == '__main__':
#     setup()


setup(
    name='advanced_boxplot_viz',  # Name of your package
    version='1.0.0',                # Version of your package
    description='Python tool for creating customizable boxplots with statistical significance bars that are ready for paper presentation',
    long_description=open('README.md').read(), 
    long_description_content_type='text/markdown',  
    author='Marcos Aguilella Fabregat',
    author_email='marcos.aguilella@idival.org',
    url='https://github.com/mafaves/advanced_boxplot_viz', 
    packages=find_packages(),  # Automatically finds all submodules and subpackages
    install_requires=[  # List of dependencies
        'matplotlib', 
        'seaborn',
        'panda',
        'numpy',
        'itertools',
        'scipy',
        'statsmodels'
    ],
    classifiers=[  # Optional classifiers to help categorize your package
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Minimum Python version required
)


