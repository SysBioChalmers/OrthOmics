# CHASSY_multiOmics_Analysis
This repository contains a collection of scripts for mining, analysing and visualizing several omics data layers:
- Transcriptomics
- Relative and absolute proteomics
- Metabolomics

These datasets were produced for three different fungal organisms:
- *S. cerevisiae*
- *K. marxianus*
- *Y. lipolytica*

Growing on different experimental conditions:

| Conditions | *S. cervisiae* CEN.PK 113 /-D |	*K. marxianus* CBS6556 |	*Y. lipolytica* W29 
| ------------- |:-------------:|:-------------:|:-------------:|
| Reference |	30°C / pH 5.5 |	30°C / pH 5.5 |	25°C / pH 5.5 |
| High temperature |	36°C / pH 5.5 |	44°C / pH 5.5 |	32°C / pH 5.5 |
| Low pH |	30°C / pH 3.5 |	30°C / pH 3.5 |	25°C / pH 3.5 |
| Osmotic stress |	30°C / pH 5.5 / KCL [1M] |	36°C / pH 5.5 / KCL [1M] |	--- |


The proteomics datasets have been incorporated to enzymatically constrained GEMs for the three organisms, available at:

| Organism | Model ID |	URL |
| ------------- |:-------------:|:-------------:|
| *S. cervisiae* |	ecYeastGEM |	https://github.com/SysBioChalmers/GECKO |
| *K. marxianus* |	ecKmarxGEM |	https://github.com/SysBioChalmers/EnzymeConstrained-Kmarx-GEM |
| *Y. lipolytica* |	ecYaliGEM |	https://github.com/SysBioChalmers/EnzymeConstrained-Yali-GEM |


- KeyWords

**Repo Category:** Data Analysis; **Utilisation:** Multi-omics/multi-organisms datasets analysis; **Field:** Metabolic engineering, Omics, Evolutionary studies;**Omic Source:** Transcriptomics, Proteomics, Metabolomics; **Taxonomy:** *S. cervisiae* CEN.PK 113 /-D,	*K. marxianus* CBS6556, *Y. lipolytica* W29  
Last update: 2018-11-13


This repository is administered by [@IVANDOMENZAIN](https://github.com/IVANDOMENZAIN), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation
### Required Software
- MATLAB (version 2013b)
- R studio (version 1.0.136 or later)
### Installation Instructions
* Clone master branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers/CHASSY_multiOmics_Analysis).

## Development Guidelines

Anybody is welcome to contribute to the development of multiOmics analysis Toolbox, but please abide by the following guidelines.

Each function should start with a commented section describing the function and explaining the parameters. Existing functions can clarify what style should be used. When making *any* changes to an existing function (`*.R`-file), change the date and name of developer near the bottom of this commented section.

### Bugfixes, new features and functions
* For any development, whether bugfixes, introducing new functions or new/updated features for existing functions: make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`.
* Make commits to this branch while developing. Aim for backwards compatibility.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries (`R` workspaces), UniProt databases, omics data files, etc.
* `doc`: updating documentation (`README` files) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: add new proteins normalization method
chore: update UniProt database for CENPK113-7D
fix: variable name corrected in `load_ProtData` function
```
More detailed explanation or comments can be left in the commit description.

### Pull request
* No changes should be directly commited to the `master` or `devel` branches. Commits are made to side-branches, after which pull requests are made for merging with `master` or `devel`.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* A merge with the master branch invokes a new release.

