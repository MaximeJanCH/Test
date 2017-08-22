
---
title: "BXD MetaData"
output: 
  html_document:
    code_folding: hide
    toc: true
    toc_float: true
    theme: spacelab
    df_print: paged
    rows.print: 6
---
![alt text](MyDiagram.png)

# Project

## Molecular data

### BXD Molecular SD

#### Description

BXD cohort for molecular quantification after sleep deprivation


#### Input

[BXD](Test.html#bxd-mouse); [Sleep Deprivation](Test.html#sleep-deprived)

#### Output

[fastq files Cortex SD](Test.html#fastq-cortex-sd); [fastq files Liver SD](Test.html#fastq-liver-sd); [Metabolite level](Test.html#metabolite-level)

### BXD Molecular NSD

#### Description

BXD cohort for molecular quantification after sleep deprivation


#### Input

[BXD](Test.html#bxd-mouse); [Control](Test.html#control)

#### Output

[fastq files Liver Control](Test.html#fastq-liver-nsd); [fastq files Cortex Control](Test.html#fastq-cortex-nsd); [Metabolite level](Test.html#metabolite-level)

## Systems Genetics of Sleep

#### Authors

Maxime Jan; Ioannis Xenarios; Paul Franken; Shanaz Diessler; Debra S.

#### Hyperlink

Link toward our publication: https://www.ncbi.nlm.nih.gov/pubmed/

Link toward Vital-IT site: https://www.vital-it.ch/

Link toward Franken lab: https://www.unil.ch/cig/en/home/menuinst/research/research-groups/prof-franken.html	

## BXD Behavioral/EEG

#### Description

Generation of data for EEG and behavioral phenotypes.
The [BXD](Test.html#bxd-mouse) are recorded for 4 days under light and dark condition. The 3rd day, mice underwent 6h of [sleep deprivation](Test.html#sleep-deprived). The recording follow-up for 2 day (recovery).
The 2 days of baseline are considered as [Control](Test.html#control) and the 2 days of recorvery are considered as the condition [sleep deprivted](Test.html#sleep-deprived)	


#### Input

[BXD](Test.html#bxd-mouse); [Control](Test.html#control); [Sleep deprivation](Test.html#sleep-deprived)

#### Output

[EEG Spectral Data](Test.html#spectral-data); [PIR motion Data](Test.html#pir-motion-sensor-data); [EEG manula annotation](Test.html#manual-annotation)

# Author

## Maxime Jan

#### Hyperlink

Adress: Maxime.Jan@sib.swiss

## Nicolas Guex

#### Hyperlink

Adress: Nicolas.Guex@sib.swiss

## Yann Emmenegger

#### Hyperlink

Adress: yann.emmenegger@unil.ch

## Paul Franken

#### Hyperlink

Adress: paul.franken@unil.ch

## Ioannis Xenarios

#### Hyperlink

Adress: Ioannis.Xenarios@sib.swiss

# File

## Normalized Counts

### CPM Liver NSD

#### Description

File name is: htseq-count_summary.Liver.txt.NSD.CPMnormalizedWITHNSDSD.txt
 
The file is generated using [Transcript normalization](Test.html#counts-normalization), each row is a [Gene](Test.html#refseq-dataset-2014) and each column is a [BXD line](Test.html#molecular-data).
 
The tissu come from the [Liver](Test.html#liver) during the [Control](Test.html#control). Counts are in Count Per Million [CPM]


#### Path

../GeneExpression/NormalizedExpression/

### CPM Liver SD

#### Description

File name is: htseq-count_summary.Liver.txt.SD.CPMnormalizedWITHNSDSD.txt
 
The file is generated using [Transcript normalization](Test.html#counts-normalization), each row is a [Gene](Test.html#refseq-dataset-2014) and each column is a [BXD line](Test.html#molecular-data).
 
The tissu come from the [Liver](Test.html#liver) during the [sleep deprivation](Test.html#sleep-deprived). Counts are in Count Per Million [CPM]


#### Path

../GeneExpression/NormalizedExpression/

## Transcript raw counts

### RawCounts.Liver.txt

#### Description

The file is generated using [Htseq-count](Test.html#htseq-count-union), each row is a [Gene](Test.html#refseq-dataset-2014) and each column is a [BXD line](Test.html#molecular-data).
 
The tissu come from the [Liver](Test.html#liver) during the [sleep deprivation or SD](Test.html#sleep-deprived) and [Control or NSD](Test.html#control).


#### Path

../GeneExpression/RawCounts/

### RawCounts.Cortex.txt

#### Description

The file is generated using [Htseq-count](Test.html#htseq-count-union), each row is a [Gene](Test.html#refseq-dataset-2014) and each column is a [BXD line](Test.html#molecular-data).
 
The tissu come from the [Cortex](Test.html#cortex) during the [sleep deprivation or SD](Test.html#sleep-deprived) and [Control or NSD](Test.html#control).


#### Path

../GeneExpression/RawCounts/

## Fastq Files

### Fastq Liver SD

#### Description

Fastq files from each [BXD lines](Test.html#bxd-mouse) generated using Illumina RNA-seq


#### Path

Archived

#### Input

[BXD molecular experiment](Test.html#molecular-data)

### Fastq Liver NSD

#### Description

Fastq files from each [BXD lines](Test.html#bxd-mouse) generated using Illumina RNA-seq


#### Path

Archived

#### Input

[BXD molecular experiment](Test.html#molecular-data)

### Fastq Cortex SD

#### Description

Fastq files from each [BXD lines](Test.html#bxd-mouse) generated using Illumina RNA-seq


#### Path

Archived

#### Input

[BXD molecular experiment](Test.html#molecular-data)

### Fastq Cortex NSD

#### Description

Fastq files from each [BXD lines](Test.html#bxd-mouse) generated using Illumina RNA-seq


#### Path

Archived

#### Input

[BXD molecular experiment](Test.html#molecular-data)

## MetaboFiles

### Metabolite level

#### Description

Metabolite level for each [BXD lines](Test.html#bxd-mouse)
 
* This file contain Metabolite data for New strains and parental and F1 during [Control](Test.html#control) and [Sleep deprivation](Test.html#sleep-deprived)
 
* A file containing also old strains is located: /home/mjan/PhD/Rdata/Metabolites_BXD_alldata.txt
 
* Two file containing Mean metabolite values can be found here: /home/mjan/PhD/Rdata/MetabolitesMean.NSD.txt (for control) and /home/mjan/PhD/Rdata/MetabolitesMean.SD.txt (for SD)


#### Path

/home/mjan/PhD/BXD_Paper_Public/Data/Metabolite_BXD.txt

## Bam Files

### Bam Files Cortex NSD

#### Description

Alignement File [Bam] for [Cortex](Test.html#cortex) during [Control](Test.html#control)


#### Path

frt.el.vital-it.ch:/scratch/cluster/monthly/mjan

## eQTLfiles

### eQTL Cortex NSD

#### Description

eQTL for [Cortex](Test.html#cortex) during [Control](Test.html#control)
Each row is a [Gene](Test.html#refseq-dataset-2014) , first column is the FDR adjusted pvalue using Rpackage qvalue
the last column is a [genetic marker](Test.html#bxd-genotypes)


#### Path

/home/mjan/PhD/Rdata/ciseQTL.Cortex.NSD.pvalcorrected.txt

## Partial Correlation Matrix

### Partial correlation matrix, liver SD

#### Path

/scratch/cluster/monthly/mjan/pcor_spearman_LiverSD_LiverSD.Rdata

## RefSeq DataSet 2014

#### Description

File that come from UCSC table browser. It was generated using RefSeq Reflat database on the 2014/01/29.


#### Path

/home/mjan/PhD/Rdata/RefSeq_20140129.gtf

#### Hyperlink

UCSC table browser: https://genome.ucsc.edu/cgi-bin/hgTables

## Hypothalamus Bam files

#### Description

RNA-seq aligned read from hypothalamus.
The Data were generated by Bernard Thorens group
Single-End Illumina reads 100bps


#### Path

Archived

#### Input

[BXD mice](Test.html#bxd-mouse)

#### Hyperlink

Bernard Thorens group: https://www.unil.ch/cig/en/home/menuinst/research/research-groups/prof-thorens.html

## Brain stem Bam files

#### Description

RNA-seq aligned read from brain stem.
The Data were generated by Bernard Thorens group
Single-End Illumina reads 100bps


#### Path

Archived

#### Input

[BXD mice](Test.html#bxd-mouse)

#### Hyperlink

Bernard Thorens group: https://www.unil.ch/cig/en/home/menuinst/research/research-groups/prof-thorens.html

## BXD variant [.vcf]

#### Description

Variant calling file generated with [GATK](Test.html#variant-calling-with-gatk)


#### Path

My Data ..

## GeneNetwork genotypes

#### Description

Genotype from GeneNetwork


#### Path

/home/mjan/PhD/Rdata/GNmm9.Fred.GNformat.geno

#### Version

2005

#### Hyperlink

Genenetwork: http://genenetwork.org/webqtl/main.py

## BXD genotypes

#### Description

Genotype of the BXD using the merge of GeneNetwork and Variant calling


#### Path

/home/mjan/PhD/Rdata/Genotype.FormatedName.geno

## Spectral Data

#### Description

Spectral data used for manual annotation and semi-automatic learning approach.
The files are transforme using Fast Fourier transform (FFT). Ask [Yann Emmenegger](Test.html#yann-emmenegger) for the processing details.


#### Path

Archived by [Paul Franken](Test.html#paul-franken)

## Manual Annotation

#### Description

manually annotated data on the [3rd day of recording (R1)](Test.html#bxd-behavioraleeg). Performed by [Yann Emmenegger](Test.html#yann-emmenegger). Manual annotation is performed on each mouse of each [BXD lines](Test.html#bxd-mouse).


#### Path

Archived by [Paul Franken](Test.html#paul-franken)

## PIR motion sensor data

#### Description

PIR motion sensor data recording activity for each mouse of each [BXD lines](Test.html#bxd-mouse)


#### Path

Archived by [Paul Franken](Test.html#paul-franken)

## Predicted EEG State

#### Description

Predicted state for each [BXD](Test.html#bxd-mouse). Prediction for 4 days every 4 seconds (86400 epochs). State are wake [w], NREM sleep [n], REM sleep [r] or artefact value of wake [1], NREM [2] or REM [3]. The 3rd day is replaced by the [manual annotation](Test.html#manual-annotation)


#### Path

Archived by [Paul Franken](Test.html#paul-franken)

## EEG and Behavioral Phenotypes

#### Description

EEG and Behavioral sleep phenotypes.
Phenotypes contain:
 	 	
* activity
* EEG
* State
 
The directory "Phenotypes/" contains the following information:
 
* [line (BXD lines)](Test.html#bxd-mouse)
* mouseID (in format "lines"_"mouseID")
* room (room used for recording)
* rec (recording time period)
* worm (worm detection yes:1 no:0)
* Phenotypes ...
 
A mean value per lines was computed and stored into the following files:
 	
* as text: AllPheno.txt
* as Rdata matrix: SleepPhenotype.Rdata


#### Path

/home/mjan/PhD/Rdata

## Phenotypes Categories

#### Description

Categories for each sleep phenotypes. The categories are related to EEG/State or Activity. The file contain also subcategories and condition related to each phenotypes (can be [bsl](Test.html#control), [rec: after sleep deprivation](Test.html#sleep-deprived) or $< SD: during sleep deprivation [Sd])	


#### Path

/home/mjan/PhD/Rdata/Phenotype.Categories.txt

# Workflow

## fastq to gene counts

#### Input

[Fastq Files](Test.html#fastq-files)

#### Output

[Bam Files](Test.html#bam-files)

#### CS

>fastq filtering <br>
STAR alignment <br>
samtools index <br>
samtools sort <br>

## Variant calling with GATK

#### Description

We use the standard GATK pipeline for variant calling on RNA-sequencing.
We use cosmic version X


#### Input

[Bam files](Test.html#bam-files); [Bam files hypothalamus](Test.html#hypothalamus-bam-files); [Bam files BrainStem](Test.html#brain-stem-bam-files)

#### Output

[VCF file](Test.html#bxd-variant-[vcf])

#### CS

>GATK RNA-seq <br>
GATK VQSR <br>
Filtering low quality variant <br>

#### Version

1.1.1

#### Hyperlink

GATK website: https://software.broadinstitute.org/gatk/

## Variant merging

#### Description

Workflow to merge the public genotype available in geneNetwork and variant called using
RNA-sequencing. The variant that are highly different from the RNA-seq data are tag as
unknown 'U'


#### Input

[VCF file](Test.html#bxd-variant-[vcf]); [GeneNetwork Genotype](Test.html#genenetwork-genotypes)

#### Output

[BXD Genotypes](Test.html#bxd-genotypes)

#### CS

>python2 Merge.py	 <br>
python2 Filtering.py <br>
python2 GenerateGeno.py <br>
python2 StrangeVariantTagging.py <br>

## eQTL dectection with FastQTL

#### Description

Generation of cis-eQTL files using FastQTL


#### Path

/home/mjan/PhD/Scripts/

#### Input

[Normalized counts](Test.html#normalized-counts); [Genotypes](Test.html#bxd-genotypes)

#### Output

[eQTL files](Test.html#eqtlfiles)

#### CS

>GenoToVcf.py <br>
bgzip & tabix <br>
Create_BED_UCSC.py <br>
bgzip & tabix <br>
ciseQTLAnalysis.sh <br>
ciseQTL_qvalueCompute.R <br>

#### Version

fastQTL.1.165.linux

#### Hyperlink

FastQTL: http://fastqtl.sourceforge.net/

## EEG & State & PIR analysis

#### Description

* EEG spectral data analysis for spectral power and spectral gain analysis.
* Activity is measure using PIR motion sensor.
* State concerns Wake [w], NREM[n] and REM[n]


#### Input

[spectral files](Test.html#spectral-data); [PIR data](Test.html#pir-motion-sensor-data); [Semi-Automatic Annotation](Test.html#predicted-eeg-state)

#### Output

[341 Phenotype Data](Test.html#eeg-and-behavioral-phenotypes)

#### Authors

Paul Franken

#### Hyperlink

Phenotypes full description: Phenotypes.docx

## Semi-Automatic EEG scoring

#### Description

 
* Supervised machine learning approach for EEG annotation.
 
* Multiple SVMs on 1s data resolution
 
* For reinstallation follow the README.txt
 
* Full description of the Workflow available in the README.txt
 
* Main bash script is: Mouse_SleepState_Prediction


#### Path

/data/ul/projects/bxd/PROD-V1.0.1a/

#### Input

[Manually annotated data](Test.html#manual-annotation); [Spectral data](Test.html#spectral-data)

#### Output

[Predicted Files](Test.html#predicted-eeg-state)

#### Version

1.0.1a

#### Authors

[Nicolas Guex](Test.html#nicolas-guex); [Maxime Jan](Test.html#maxime-jan)

# Script

## Counts normalization

#### Description

Script used for transcript normalization using edgeR


#### Path

../Scripts/Limma_RNA.R

#### Input

[Raw counts](Test.html#transcript-raw-counts)

#### Output

[Normalized counts](Test.html#normalized-counts)

## htseq-count union

#### Description

Script to generate count file from multiple bam file alignment


Script take as input a list of the processed bam file


#### Path

../Scripts/GenBSUB_htseq-count.py

#### Input

[Bam Files](Test.html#bam-files); [RefSeq file](Test.html#refseq-dataset-2014)

#### Output

[Raw counts](Test.html#transcript-raw-counts)

#### Arguments

-q -s reverse -t exon -m union

#### Version

0.5.4p3

## Genes partial correlations

#### Description

A Rscript to generate a partial correlation matrix.
 
* Example using the script to compute pcor for gene in cortex during Control (Cortex NSD vs Cortex NSD): Matrix<-BXDpcor(CNSD,cis_CNSD,CNSD,cis_CNSD,GenotypeMatrix,CORES=30,method="pearson")
 
* The input file (Expression matrix: CNSD and cis-eQTL matrix: cis_CNSD) need to be formated to contain the same gene, You can see an example here: /home/mjan/PhD/Analysis/Partial_Correlation/GetPartialCorrelationForBXD.R
 
* The method can be "spearman" or "pearson"
* "parallel","reshape2","matrixStats" and "RcppEigen" packages are required
* If you want to change between interaction or additive effect, edit model matrix in the CorResiduals() function


#### Path

/home/mjan/PhD/Scripts/BXDpcor_V1.0.1.R

#### Input

[Expression Data](Test.html#normalized-counts); [cis-eQTL File](Test.html#eqtlfiles)

#### Output

[Partial correlation Matrix](Test.html#partial-correlation-matrix)

#### Version

1.0.1

# Population

## BXD mouse

#### Description

Used to generate [Molecular data](Test.html#molecular-data) and [Behavioral data](Test.html#bxd-behavioraleeg)   


# Tissu

## Cortex

#### Description

Cortex sample extracted for RNA-seq sequencing


## Liver

#### Description

Liver sample extracted for RNA-seq sequencing


## Blood Metabolism

#### Description

Blood Metabolism quantification


# Condition

## Control

#### Description

Control sample from the [BXD](Test.html#bxd-mouse) Experiments


## Sleep deprived

#### Description

6h of sleep deprivation during the 3rd day of recording during the first our of light phase [ZT0-6]


# Rmarkdown

## Heritability calculation

#### Description

Rmarkdown for heritability calculation, details in the Rmarkdown


#### Path

/home/mjan/PhD/BXD_Paper_Figures/Heritability.html

#### Input

[Sleep phenotypes](Test.html#eeg-and-behavioral-phenotypes); [Sleep categories](Test.html#phenotypes-categories)

## Metabolite Ctl vs SD

#### Description

Rmarkdown for differential metabolite level, details in the Rmarkdown


#### Path

/home/maxime/Link_To_Cluster/BXD_Paper_Public/Analysis/Differential_Expression_Metabolites.html

#### Input

[Metabolite level](Test.html#metabolite-level)

