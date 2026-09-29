# :leaves: Ophrys aveyronensis - vitorica speciation :dna:

Scripts and files to generate results and output of the *Ophrys vitorica* - *O. aveyronensis* speciation study. 

## General scripts

:file_folder: Files:

:bar_chart: Scripts:
- `snakefile.py` : general script 

## Read filtering, trimming and mapping

:file_folder: Files:

:bar_chart: Scripts:
- `extract_fastp.py` : get content of reads after reads trimming and filtering
- `get_mapping_statistics.R` : get mapping statistics

## Estimation of recombination rate

:file_folder: Files:

:bar_chart: Scripts:
- `iSMC/create_tab.py` : create tab file to run iSMC
- `iSMC/make_input_iSMC.R` : create input to run iSMC
- `get_iSMC_windows.py` : get mean recombination rate inferred by iSMC per custom genomic window size
- `plot_iSMC.R` : plot recombination rate along the genome

## :wrench: Tools needed

* [fastp v.1.3.4](https://github.com/OpenGene/fastp)
* [bwa v.0.7.17](http://bio-bwa.sourceforge.net/bwa.shtml)
* [samtools v.1.19](https://github.com/samtools/samtools/releases/)
* [bcftools v.1.19](https://samtools.github.io/bcftools/bcftools.html)
* [vcftools v.0.1.16](https://vcftools.github.io/index.html)
* [snakemake v.9.4.0](https://github.com/snakemake/snakemake)
* [R v.4.4.1](https://cran.r-project.org/bin/windows/base/old/4.4.1/)
* [iSMC v.0.0.23](https://github.com/gvbarroso/iSMC)
* [mosdepth v.0.2.6](https://github.com/brentp/mosdepth)
* [ngsParalog v.1.3.4](https://github.com/tplinderoth/ngsParalog)
* [getorganelle v.1.7.5.0](https://github.com/kinggerm/getorganelle)
