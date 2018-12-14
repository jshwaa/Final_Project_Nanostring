The NanoString nCounter technology is a medium-throughput method for the analysis of differentially expressed __pre-determined genes__ ("NanoString panels") that is unique in its direct, digital readout of mRNA counts from sample tissue. 
This juxtaposes the technique to those that require amplification steps (RNAseq), which may introduce greater false-positives, or those that depend on intensity measures for relative abundance (microarrays). As such, a number of methods have been designed to analyze this data.

Here, I use the R package NanoStringNorm, developed by Waggot et al(ref) and available on CRAN, to process, perform QC checks, and analyze NanoString gene expression counts in striatal brain tissue collected from wild-type or R6/2 mice, a model of Huntington's disease (HD).

This can be quickly done using R alone, and will be outlined here as command line input in the R console for readability.
```
> getwd()
[1] "C:/Users/Josh/Documents"
> install.packages("NanoStringNorm")
> library(NanoStringNorm)
```

## Importing Data

Read the Results sheet of the nCounter RCC aggregate sample excel spreadsheet containing the counts into R. 
This sheet also includes counts from housekeeping genes, negative probes and positive probes of varying concentration with each sample for later normalization:
```
> NanoString.mRNA <- read.xls.RCC("ee282_final/RCCCollector_R62.xls", sheet = 1)

You have chosen to import worksheet 1 named RCC Collection. Does that sound correct?
The other sheet names are: 
1:RCC Collection
2:Reporter Order

There were 6 samples imported. 
Note that spaces in sample names will be replaced by dots. 
The first and last 3 sample names found in the dataset are:
NTVeh82str NTVeh83str NTVeh84str HDVeh124str HDVeh114str HDVeh113str

There were 782 genes imported with the following Code Class breakdown:
  Endogenous Housekeeping     Negative     Positive 
         750           20            6            6 
```         

Extract the counts data frame from the list NanoString.mRNA, which also contains a "headers" data frame, with ```NanoString.mRNA <- NanoString.mRNA$x```.

## Gene expression normalization

Before normalizing the data and running QC, assign binary 'trait' values (i.e. 1 or 2) corresponding to experimental groups in a matrix:
```
> sample.names <- names(NanoString.mRNA)[-c(1:3)]; #gets all headers except annotation, which are just the sample names
> WT <- rep(1, times = (ncol(NanoString.mRNA)-3));
> WT[grepl('NT',sample.names)] <- 2;
> HD <- rep(1, times = (ncol(NanoString.mRNA)-3));
> HD[grepl('HD',sample.names)] <- 2;
> traits.geno <- data.frame(row.names = sample.names, WT = WT, HD = HD);

> traits.geno
           WT HD
NTVeh82str   2  1
NTVeh83str   2  1
NTVeh84str   2  1
HDVeh113str  1  2
HDVeh114str  1  2
HDVeh124str  1  2
```

Now use the normalization feature of the package, taking the geometric mean of positive and housekeeping controls as summary values for CodeCount and SampleContent, respectively, 
and performing a Background correction of mean +- 2 standard deviations to normalize the data. Any samples with housekeeping gene levels three standard deviations from the overall mean will be flagged for review.

```
> NanoString.mRNA.norm <- NanoStringNorm(x = NanoString.mRNA, CodeCount = 'geo.mean', Background = 'mean.2sd', SampleContent = 'housekeeping.geo.mean', round.values = TRUE, take.log = TRUE, traits = traits.geno);

##############################
### NanoStringNorm v1.2.1 ###
##############################

There are 6 samples and 750 Endogenous genes 

Background: After correction 6 samples and 554 
	Endogenous genes have less than 90% missing. 

log: Setting values less than 1 to 1 in order to calculate the log in positive space.
```

## Quality Control
NanoStringNorm's "Plot" feature generates a number of essential diagnostics for a quick examination of data integrity.
```
> png('NanoStringNorm_Example_Plots_%03d.png', units = 'in', height = 6, width = 6, res = 250, pointsize = 10);
> Plot.NanoStringNorm(x = NanoString.mRNA.norm, label.best.guess = TRUE, plot.type = c('cv', 'mean.sd', 'RNA.estimates', 'missing', 'norm.factors', 'positive.controls', 'batch.effects'));
> dev.off();
```

IMAGES

## Differential Expression
To examine differential expression, NanoStringNorm also generates volcano plots of the log p-value vs. log fold-change:
```
> png('NanoStringNorm_Example_Plots_Volcano.png', units = 'in', height = 6, width = 6, res = 250, pointsize = 10);
> Plot.NanoStringNorm(x = NanoString.mRNA.norm, label.best.guess = TRUE, plot.type = c('volcano'), title = FALSE);
> dev.off();
```

Here, the dashed lines indicate p=0.05 and after Bonferroni correction. Point size is proportional to expression level. 
