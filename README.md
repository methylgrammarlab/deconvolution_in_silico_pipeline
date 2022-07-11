# In-silico experiment pipeline

The pipeline uses two types of configuration files:
1) config.json - general parameters such as run name, models to use, number of replicates
2) tabular config - each row represents an in-silico mixture. 
The tabular config files used to create paper figures can be found under "resources".

Two in-house packages are required to run this pipeline: 
1) https://github.com/methylgrammarlab/epiread-tools (install this first)
2) https://github.com/methylgrammarlab/deconvolution_models 

Deconvolution is based on two inputs: a reference atlas and a mixture sample. The reference atlas contains information
on known cell types. The mixture contains unknown proportions of the reference cell types. Deconvolution models will 
use the atlas to determine proportions in the mixture.

The in-silico experiments do not use the entire genome. Rather, they use a defined list of genomic 
intervals, hereby referred to as "regions". These regions should be differentially methylated across the reference cell
types to enable deconvolution. There are many approaches to selecting an appropriate set of regions. Here, 
if no user regions are supplied, [Tissue Informative Markers](https://github.com/christacaggiano/celfie) are called.

Regions file should be BED-formatted and have no header. Regions may overlap. If they do - the overlapping section
will be read twice. If this is not the desired behaviour, use [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to avoid overlap. 

For CelFiE and CelFiE+, the atlas format should be METH, COV, as in [CelFiE](https://github.com/christacaggiano/celfie).
A small pipeline to create such an atlas from bedgraph files is provided here. The atlas should include individual
CpGs (this is the input for CelFiE+). Summing of methylation and coverage per region for CelFiE happens in [deconvolution_models](https://github.com/methylgrammarlab/deconvolution_models). 

For epistate models, an epistate atlas is required. This is done separately, [here](https://github.com/methylgrammarlab/proj-epireads).

Whole-genome files are required to call TIM regions. If user-supplied regions are used, the atlas may only cover the regions.

| Input        | No atlas                                    | User atlas  |
|--------------|---------------------------------------------|-------------|
| No regions   | - call atlas from bedgraph<br/> - call TIMs | - call TIMs |
| User regions | - cell atlas from epiread                   |             |

## Usage

Clone the repository, adjust config.json as necessary. From the project root simply run
```
snakemake --cores 1
```

Output files will be created under "results". 