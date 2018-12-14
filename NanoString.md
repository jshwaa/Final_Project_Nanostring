The NanoString nCounter technology is a medium-throughput method for the analysis of differentially expressed __pre-determined__ genes ("NanoString panels") that is unique in its direct, digital readout of mRNA counts from sample tissue. 
This juxtaposes the technique to those that require amplification steps (RNAseq), which may introduce greater false-positives, or those that depend on intensity measures for relative abundance (microarrays). As such, a number of methods have been designed to analyze this data.

Here, I use the R package NanoStringNorm, developed by Waggot et al (1) and available on CRAN, to process, perform QC checks, and analyze NanoString gene expression counts in striatal brain tissue collected from adult wild-type or R6/2 mice, a model of Huntington's disease (HD). I was interested in looking at changes in inflammatory and protein clearance genes contained in the Cancer panel due to HD, in the brain region most susceptible to the disease.

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
and performing a Background correction of mean +- 2 standard deviations of the negative control values to normalize the data. Any samples with housekeeping gene levels three standard deviations from the overall mean will be flagged for review.

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
This indicates that, as we chose a more stringent background correction, 196 genes on our panel were indistinguishable from background (and thus zero-valued after normalization) in >90% of samples.

## Quality Control
NanoStringNorm's "Plot" feature generates a number of essential diagnostics for a quick examination of data integrity.
```
> png('NanoStringNorm_Example_Plots_%03d.png', units = 'in', height = 6, width = 6, res = 250, pointsize = 10);
> Plot.NanoStringNorm(x = NanoString.mRNA.norm, label.best.guess = TRUE, plot.type = c('cv', 'mean.sd', 'RNA.estimates', 'missing', 'norm.factors', 'positive.controls', 'batch.effects'));
> dev.off();
```

![NSM_mean.sd](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_mean.sd.png?raw=true)
Standard deviation of gene expression across samples plotted against mean expression shows that negative controls are undetectable while housekeeping genes have high means and little variation, as is expected. There is variation in a number of genes that could be due to genotype effects. 


![NSM_cv](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_cv.png?raw=true)
Normalization highly left-shifted the coefficient of variation, lowering the overall variation in the data across samples.

![NSM_missing](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_missing.png?raw=true)
This plot of missing gene counts following normalization indicates that the expression of a moderate number of genes within this panel was undetectable in the striatal tissue. However, if we run a more liberal correction, using only the mean the negative probes as background (as opposod to mean and two standard deviations):

```
>NanoString.mRNA.norm2 <- NanoStringNorm(x = NanoString.mRNA, CodeCount = 'geo.mean', Background = 'mean', SampleContent = 'housekeeping.geo.mean', round.values = TRUE, take.log = TRUE, traits = traits.geno);

##############################
### NanoStringNorm v1.2.1 ###
##############################

There are 6 samples and 750 Endogenous genes 

Background: After correction 6 samples and 712 
	Endogenous genes have less than 90% missing. 

log: Setting values less than 1 to 1 in order to calculate the log in positive space.

> png('NanoStringNorm_Example_Plots_%03d.png', units = 'in', height = 6, width = 6, res = 250, pointsize = 10);
> Plot.NanoStringNorm(x = NanoString.mRNA.norm2, label.best.guess = TRUE, plot.type = c('missing'));
> dev.off();
```

We get a ~20% increase in the number of nonzero gene counts after normalization:

![NSM_missinglib](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSM_missinglib.png?raw=true)
As NanoString is less prone to amplification bias and increased false positives as seen in RNAseq/microarrays, a more liberal background correction process may be argued for in subsequent analysis.  

![NSM_RNAcontent](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_RNAcontent.png?raw=true)
Returning to the origin normalization, RNAcontent plots of highly expressed house-keeping genes vs. other highly expressed endogenous genes show that no samples deviate far from the best-fit line, indicating good sample distribution with no outliers.

![NSM_BatchEffects](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_BatchEffects.png?raw=true)
This plot indicates that there are indeed batch effects from the binary trait groups supplied to the normalization function earlier. For instance, it appears that HD samples have lower over RNA content, mean gene expression, and more missing genes compared to WT. 

![NSM_Norm](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_Norm.png?raw=true)
Plots of normalizaiton factors supplied to the function (positive controls/negative controls) and adjusted RNA content allow for a per sample look at potential outliers that could be driving the previous graph's batch effects. Here, as none of the samples deviate too much from the mean (~<50%) it would appear that the batch effects are due to treatment and not individual sample outliers. The only potential confounder (sample 1) could be removed and the data reanalyzed to confirm this. 


![NSM_Controls](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSN_Controls.png?raw=true)
Finally, QC provided by NanoStringNorm include plots of observed vs expected concentrations of controls ran with each sample. Positive controls should hvae similarly linear slopes across samples, while negative controls should not trend in any direction and instead have values near zero. By these parameters, all of the samples run through this NanoString appear sufficient for downstream analysis.


## Differential Expression
To examine differential expression, NanoStringNorm also generates volcano plots of the log p-value vs. log fold-change of the groups defined earlier as binary traits. With ```label.best.guess = TRUE``` the package will label the most significant genes (# of which is equal to label.n, here defined as 20) in descending order, with fold-change being used for tie-breaking: 
```
> png('NanoStringNorm_Example_Plots_Volcano.png', units = 'in', height = 6, width = 6, res = 250, pointsize = 10);
> Plot.NanoStringNorm(x = NanoString.mRNA.norm, label.best.guess = TRUE, label.n = 20, plot.type = c('volcano'), title = FALSE);
> dev.off();
```

Here, the dashed lines indicate p=0.05 and after Bonferroni correction. Point size is proportional to expression level. Fold change is measured as HD genes relative to WT:
![NSM_volcano](https://github.com/jshwaa/Final_Project_Nanostring/blob/master/Images/NSM_volc.png?raw=true)

As we can see, the most significantly different genes between WT and HD striatal tissue don't display massive fold changes in expression. Additionally, of all the top 25 significantly different genes, few if any are inflammation-related. Indeed, markers for the primary immune cell in the brain, microglia, are not significantly elevated in this model as they are in other neurodegenerative contexts (i.e. Alzheimer's, Parkinson's - these markers include Cx3cr1, trem2, Csf1r, specific in the brain for microglia). This indicates a lack of neuroinflammatory response in HD mice, which I have confirmed separately with immunohistochemistry on the same tissue. While further analysis is called for, the downregulation of Adora2a (gene #15), a hallmark of HD striatal tissue, supports the valiity of these results.

1. Waggott D, Chu K, Yin S, Wouters BG, Liu FF, Boutros PC. NanoStringNorm: an extensible R package for the pre-processing of NanoString mRNA and miRNA data. Bioinformatics. 2012;28(11):1546-8.
