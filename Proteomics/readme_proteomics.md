# LC-MS/MS PROTEOMICS

## Results

This directory contains the proteomics data of the yeast strains _S. cerevisiae_ CEN.PK113-7D (**sce**), _K. marxianus_ CBS 6556 (**kma**) and _Y. lipolytica_ W29 (**yli**). 

The results obtained from the differental analysis pipeline can be found in the **.xlsx**  files:
- sce_Proteins_DE_BothMSRuns.xlsx 	
- kma_Proteins_DE_BothMSRuns.xlsx 	
## Data
**Data** directory contains the detected proteins in an Orbitrap Fusion^TM^ Lumos^TM^ Tribrid^TM^ Mass Spectrometer (Thermo Fischer Scientific) and identified by the search engine X!Tandempipeline [Langella _et al._ 2017](https://pubs.acs.org/doi/10.1021/acs.jproteome.6b00632), _J Proteome Res_. 

Label-free proteomics data is in **relative** and **absolute** quantification. 

## Relative quantification
Relative quantification can be found in the **txt** files:

The methods of quantification for each strain are:

* **emPAI**
The exponentially modified protein abundance index (emPAI) was proposed by [Ishihama _et al.,_ (2005)](https://www.mcponline.org/content/4/9/1265.long) as:

emPAI  = (10<sup>PAI</sup>) – 1

PAI = the number of observed peptides counts per protein over the number of observable peptides per protein

* **Spectral Counts**
Defined as the number of MS2 spectra assigned to a protein ([Liu _et al.,_, 2004](https://www.ncbi.nlm.nih.gov/pubmed/15253663)).

* **XIC** 
The average of the eXtracted Ion Current (XIC) or MS1 intensities of all the peptides associated to a protein.


## Absolute quantification
Absolute quantification was determined using the external standard UPS2 (Sigma) to estimate the proteins abundance in **fmol**.  These values were generated using two different quantification methods: Intensity Based Absolute Quantification (**iBAQ**) and Normalized Spectral Abundance Factor (**NSAF**).

iBAQ is defined as the sum of the MS1 intensities of all the peptides associated to a protein divided by the number of theoretically observable peptides ([Schwanhäusser _et al._, 2011](https://www.nature.com/articles/nature10098)).

NSAF was defined by [Zibailov _et al._ (2006)](https://pubs.acs.org/doi/abs/10.1021/pr060161n) as: 

NSAF<sub>k</sub> = (SC/L)<sub>k</sub> /  SUM (SC/L) <sub>N</sub>

where:

_k_ =  a given protein

_N_ = total quantified proteins

_SC_ = Spectral counts 

_L_ = protein length in aa
