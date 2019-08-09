# Table of contents

[Session 1: Variant calling with GATK](https://github.com/prenderj/variant_tutorial#session-1-variant-calling-with-gatk)

[Session 2: Variant annotation and filtering](https://github.com/prenderj/variant_tutorial#session-2-variant-annotation-and-filtering)

[Session 3: Association study with plink](https://github.com/prenderj/variant_tutorial#session-3-association-study-with-plink)

# Session 1: Variant calling with GATK

## Dataset
After starting the interactive session the data you need for this session can be copied from */export/data/ilri/bioinformatics/workshop_dnaseq/Session1*. 

```
interactive -c 4
```
We first go to the home directory, then copy the required data and move to this new directory
```
cd
cp -R /export/data/ilri/bioinformatics/workshop_dnaseq/Session1/ ~/
cd ~/Session1
```

This session will take you through the process of calling variants with GATK using cattle whole genome sequencing data. There are 2 bam files in the Session1 folder you just copied (aligned with BWA), each corresponding to a different cow. These contain reads mapping to an approximately 4Mb region centred on the cow Leptin gene on chromosome 4 (you can view the region here: http://www.ensembl.org/Bos_taurus/Location/Overview?r=4%3A91249874-95266624).

Programs required:
*	R (https://cran.r-project.org/)
*	Samtools (http://www.htslib.org/doc/samtools.html)
* Picard (https://broadinstitute.github.io/picard/)
*	GATK (https://software.broadinstitute.org/gatk/)
*	VCFtools (https://vcftools.github.io/index.html)

These are already installed on the ILRI cluster. To make them available type:
```
module load samtools/1.3.1
module load picard/2.8.2
module load gatk/3.7.0
module load R/3.6
module load vcftools/0.1.15
```

## 1. Preparing the reference genome for use with GATK
Earlier we showed you how to map your reads against a reference genome, and these bams were aligned with BWA. To run GATK it is first necessary to create partner files to the fasta format reference genome used to do the alignment; a fasta index file (.fai) and a dictionary file (.dict). We have provided the reference genome used to align these data (Bos_taurus.UMD3.1.dna.toplevel.fa). The partner files can be generated from this.

First use samtools to generate the fasta index file:
```
samtools faidx Bos_taurus.UMD3.1.dna.toplevel.fa
```

Next use Picard to generate the sequence dictionary:
```
picard CreateSequenceDictionary R=Bos_taurus.UMD3.1.dna.toplevel.fa \
 O=Bos_taurus.UMD3.1.dna.toplevel.dict
```
## 2. Sorting and indexing the bam file
Now sort and index your bam files using Samtools. This allows downstream programs to access its contents more quickly.
```
samtools sort cow1.bam > cow1_sorted.bam
samtools index cow1_sorted.bam
```
## 3. Mark duplicates
As discussed PCR duplicates can cause problems when calling variants. Using Picard we can flag to GATK which reads are potential duplicates.
```
picard MarkDuplicates I=cow1_sorted.bam O=cow1_dupMarked.bam \
 M=cow1_dup_metrics.txt
```
Use samtools to index this new bam file.
```
samtools index cow1_dupMarked.bam
```
**_Question 1: What proportion of reads were duplicates?_**
## 4. Base Quality Score Recalibration (BQSR)
BQSR involves adjusting the base quality scores so that they more accurately represent the probability of errors in the reads.  BQSR requires a set of known variants so that it can filter these out. For this we are going to use the variants from this study https://www.ncbi.nlm.nih.gov/bioproject/262770. We can download this in VCF format using the wget program. As these bams only contain reads mapping to a region of chromosome 4 we can just download known variants on this chromosome.
```
wget ftp://ftp.ebi.ac.uk/pub/databases/nextgen/bos/variants/population_sites/UGBT.population_sites.UMD3_1.20140307.vcf.gz
```
As the bams we are working with only contain reads mapping to a region of chromosome 4 we can first extract just the variants on chromosome 4 from this larger file using vcftools.
```
vcftools --gzvcf UGBT.population_sites.UMD3_1.20140307.vcf.gz --chr 4 --recode \
 --out UGBT.population_sites.UMD3_1.20140307.chr4
```
Can see we specified four arguments. Have a look here https://vcftools.github.io/man_latest.html to see what each of these arguments mean.

We can then use this file of known sites in BQSR.
```
GenomeAnalysisTK -T BaseRecalibrator \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -I cow1_dupMarked.bam -knownSites UGBT.population_sites.UMD3_1.20140307.chr4.recode.vcf \
 -o cow1_recal_data.table
```
And apply recalibration
```
GenomeAnalysisTK -T PrintReads \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -I cow1_dupMarked.bam -BQSR cow1_recal_data.table \
 -o cow1_recal.bam
```
## 5. Single sample GVCF calling
As this bam file only contains reads mapped to a small region of the cow genome, to speed things up we are going to restrict variant calling to just this region. To do this we first need to make a target file, defining the region of interest. Using a text editor (e.g. nano) create a new file, targetRegion.bed, containing the coordinates of the region as three tab separated columns. Chromosome, start position, and end coordinate i.e.:

4 91249874  95266624 

We can then use this targetRegion.bed with the –L flag to restrict calling to this region. If you wanted to call variants across the whole genome you would just leave this out.
```
GenomeAnalysisTK -T HaplotypeCaller \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -L targetRegion.bed -I cow1_recal.bam \
 --emitRefConfidence GVCF --dbsnp UGBT.population_sites.UMD3_1.20140307.chr4.recode.vcf \
 -o cow1.g.vcf
```
**IMPORTANT. The steps 2-5 above would normally be run for each of your samples to generate one g.vcf file per sample. However to save time we have prepared the g.vcf file for the remaining sample for use with the next step.**
## 6. Raw variant calls
Run GenotypeGVCFs over the 2 g.vcf files to get a list of raw (unfiltered) variants. Each of the separate g.vcf files is specified with the --variant flag.
```
GenomeAnalysisTK -T GenotypeGVCFs \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -L targetRegion.bed \
 --variant cow1.g.vcf \
 --variant cow2.g.vcf \
 -o variants.vcf
```
##7.	Variant Quality Score Recalibration (VQSR)
The raw variant calls coming from the command above are expected to contain a large number of false positives. To try and identify the good from the bad variants we can use VQSR. This involves first building a SNP recalibration model using the following command. 
```
GenomeAnalysisTK \
    -T VariantRecalibrator \
    -R Bos_taurus.UMD3.1.dna.toplevel.fa \
    -input variants.vcf \
    -resource:bov1000G,known=false,training=true,truth=true,prior=12.0 bov1000G.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 UGBT.population_sites.UMD3_1.20140307.chr4.recode.vcf \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R
```
Each of the -an parameters specify which features we want to use to try and separate good from bad variants. The bov1000G.vcf is a list of variants generated by the 1000 bulls sequencing consortium. We are assuming here that they are good quality variants.  We could also, for example, use the set of variants on the Illumina bovine HD array as an even better truth set. The dbSNP VCF file we downloaded earlier is expected to contain many false positive variants and so is not used to train the model.

**_Question 2: How many novel variants were in the top quality tranche?_**

**_Question 3: How do the Transition/Transversion ratios differ between known and novel variants in each tranche?_**

**IMPORTANT. The command above just runs the recalibration model on SNPs. To run it on the indel calls you would normally need to change the -mode flag and provide appropriate training and truth datasets. However for the tutorial we will just continue with the SNPs only.**

Next we modify the vcf file containing our variants so that the top quality variants are indicated with a PASS label in the VCF file.

```
GenomeAnalysisTK -T ApplyRecalibration \
    -R Bos_taurus.UMD3.1.dna.toplevel.fa \
    -input variants.vcf \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o variants_VQSR.vcf
```
Congratulations. You have now called a set of variants from your bam file.

# Session 2: Variant annotation and filtering
## Dataset
In this session we will practice annotating and filtering variants using the variants_VQSR.vcf file from the previous session.

Programs required:
*	VEP (www.ensembl.org/vep)
*	VCFtools (https://vcftools.github.io/index.html)

These are already installed on the cluster. To make them available type:
```
module load vep/87
module load vcftools/0.1.14
```
##8.	Variant annotation with VEP
VEP (variant effect predictor) is a commonly used tool to annotate genetic variants. To speed up annotation we will run VEP specifying the location of a VEP local cache that has been installed on the ILRI cluster at /export/data/bio/vep/. This contains information on the location of genes etc. that Ensembl uses to annotate the variants, and means that VEP doesnt need to connect to the Ensembl servers to get this information. By specifying the --merged flag we can get annotations for both Ensembl and RefSeq genes.
```
variant_effect_predictor.pl --cache \
 --dir_cache /export/data/bio/vep/ --species bos_taurus \
 --vcf -i variants_VQSR.vcf -o variants_VEP.vcf --merged
```
Have a look at the html output file.

**_Question 4: How many missense variants were identified?_**

**_Question 5. Have a look at the manual page (http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) and work out how to add SIFT scores to the output vcf. You may need to delete the previous output file if using the same name. What proportion of the missense variants are predicted to be deleterious with higher confidence?_**

Ensembl also has a REST API. This allows variant annotation to be incorporated into your own programs/scripts. You can also use this to quickly get the annotation of individual variants on the command line. For example:
```
wget -q --header='Content-type:application/json' 'https://rest.ensembl.org/vep/cow/region/4:94947404-94947404:1/A?' -O -
```
This returns the consequence of a change to an A allele at position 4:94947404 in the cow genome. If you specify a filename after the -O parameter then the results, in json format, will be output to this file. Specifying "-" on the other hand leads the results to simply be printed to the screen. 


##9. Variant filtering with VEP
VEP comes with an accompanying script that can be used to filter for specific types of variants. For example, it is possible to use the filter_vep.pl script to just extract the missense variants from the VEP annotated file.
```
filter_vep.pl -I variants_VEP.vcf -format vcf -o missense_only.vcf \
 -filter "Consequence matches missense"
```
**_Question 6: Get a count of the number of coding variants in the RBM28 gene (hint take a look here http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html)_**

##10. Variant filtering and metrics with VCFtools
VCFtools is a program for manipulating VCF files. As well as getting various summary statistics it can also be used to convert VCF files to different formats (e.g. PLINK format). For example, to get the allele frequency of each variant:
```
vcftools --vcf variants_VEP.vcf --freq --out freqs
```
To get counts of each type of change e.g. (A<->C) and the total number of transitions and transversions
```
vcftools --vcf variants_VEP.vcf --TsTv-summary --out TsTv
```
**_Question 7: Which base changes (e.g. A<->C, A<->G etc) are observed most often in this data? Why might these be most common (hint: https://en.wikipedia.org/wiki/Deamination)_**

**_Question 8: How can you output a new vcf file with only the genotypes for the first cow (with ID 2637) (https://vcftools.github.io/man_latest.html)?_** 

# Session 3: Association study with plink
## Dataset
The data you need for this session can be copied from */export/data/ilri/bioinformatics/workshop_dnaseq/Session3*
```
cp -R /export/data/ilri/bioinformatics/workshop_dnaseq/Session3/ ~/
cd ~/Session3/
```
In this session we will use human data from the 1000 genomes dataset in a simulated GWAS analysis. The file we are using contains genotypes for approximately 2500 individuals from 26 different global populations. The file was originally downloaded from the 1000 genomes ftp site here (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

Programs required:
*	Plink (https://www.cog-genomics.org/plink2/)
* R

To load them type:
```
module load plink/1.9
module load R/3.3.2
```
##11. Converting VCF to plinks binary format
Although PLINK can read VCF files, downstream analyses will be quicker if we first convert the data to plink's binary file format. For this tutorial we will also restrict the analysis to just British and Kenyan (Luhya in Webuye) individuals in this cohort and variants with a minor allele frequency of at least 5%. We also exclude genotypes with a quality score less than 40 to make sure low quality genotypes do not affect our analyses.
```
plink --biallelic-only strict --vcf-min-gq 40 --pheno 1kg.ped --mpheno 4 --vcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep GBR_LWK_ids.txt --maf 0.05 --make-bed --out plinkFormatted --update-sex 1kg.ped 3
```
PLINK requires a ped file containing each samples phenotype individual. Have a look at the 1kg.ped file. The first two columns are the family and individual ids. As these individuals are all unrelated we just use their individual id as their family id. 
To specify the individuals to keep we just provide a list of ids in the GBR_LWK_ids.txt file.

##12. Summary statistics with PLINK 
Like VCFtools PLINK can generate various summary statistics. For example to get the allele frequency of each variant you can type:
```
plink --bfile plinkFormatted --freq
```
This reads in the binary files we generated with the previous command and outputs the frequencies of each variant. If you do not specify an output file with PLINK commands, the results file will be given the default name, in this case plink.frq.

**_Question 9: A common filter in GWAS is excluding variants that show substantial deviation from Hardy-Weinberg equilibrium (https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle). Have a look at the plink manual (https://www.cog-genomics.org/plink2/) and work out how to get a HWE p value for every variant._**
##13. Studying population stratification
A good place to start in most association studies is to check for population stratification in your data. As discussed population stratification can skew your analyses if not properly accounted for.
With large datasets you often want to prune your variants down to a subset that are largely in linkage equilibrium first. This reduces redundancy as well as making the data more manageable. To do this in PLINK:
```
plink --bfile plinkFormatted --indep-pairwise 500 50 0.2
```
The first parameter specifies the window size (in KB) within which pruning occurs and the second the step size by which the windows are slid by. The final parameter specifies the level of LD (r squared value) at which pruning is carried out. The lower this value the more stringent the pruning.
This command outputs a list of pruned variants. To generate new binary files just containing these variants type:
```
plink --bfile plinkFormatted --extract plink.prune.in --make-bed --out plinkFormatted_pruned
```
We can then use this pruned dataset with the PCA function in plink to get principal components:
```
plink --bfile plinkFormatted_pruned --pca header
```
We can examine the PCA results in R. For example we can plot the first two principal components using the following R code (the British samples are in rows 1 to 91 and the Kenyan samples in rows 92 to 190)

First start R:
```
R
```
Then run the following R code
```{r}
#read in the plink PCA results
dat<-read.table("plink.eigenvec", header=T)

#create an empty plot
plot(dat[,3], dat[,4], xlab="Principal component 1", ylab="Principal component 2", type="n")
#plot the two populations separately and add a legend
points(dat[1:91,3], dat[1:91,4],col="blue")
points(dat[92:190,3], dat[92:190,4],col="red")
legend("topright", c("British", "Kenyan"), col=c("blue", "red"), pch=1)

```
This code reads in the eigenvec file produced by plink and plots the first two principal components against each other (found in columns 3 and 4), colouring samples by their country of origin. If you look at the plot you can see two main clusters, with individuals separating by country of origin. If the incidence of the phenotype differs between samples this will cause problems, and variants will be associated to the phenotype due to population differences in the variant frequencies rather than real associations between the genetic locus and the disease. This is an extreme example as GWAS would rarely be performed on such genetically divergent populations together. However, in this tutorial we will carry on with both populations.

To quite R you can just type:
```
q()
```

##14. Basic association analysis
To perform a basic case control study in plink:
```
plink --bfile plinkFormatted --assoc fisher
```
We can then plot a Manhattan and Q-Q plot in **R**:
```{r}
#read in the association results
dat<-read.table("plink.assoc.fisher", header=T)
#open pdf to write graphs to

plot(dat$BP, -log10(dat$P), pch=20, col="grey" , xlab="Position on chr22", ylab=expression(-log[10](italic(p))))
#extract just those rows where p value was less than standard threshold
signif<-dat[which(dat$P < 5e-08),]
#plot these significant points in a different colour and indicate the threshold with a dashed line
points(signif$BP, -log10(signif$P),pch=20, col="red")
abline(h=7.3, col="red", lty=2)

#plot Q-Q plot
#log transform p values and generate expected distribution of p values
obs <- -log10(sort(dat$P,decreasing=F))
exp <- -log10(1:length(obs)/length(obs))
#plot observed versus expected p values
plot(exp,obs,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),
	xlim=c(0,max(c(obs,exp))),ylim=c(0,max(c(obs,exp))))
#add y=x line
lines(exp,exp,col="red")

```
In the Manhattan plot we have log transformed the p values so that small p values now have larger values on the plot i.e. the higher they are the more significant. Those above the p=5x10-8 significance threshold commonly used in GWAS are indicated in red. As can be seen a lot of sites are apparently significant. This though is how not to do a GWAS. As we have shown this data shows substantial population stratification that will likely be confounding these results. When performing GWAS you control for confounders by fitting what can be termed covariates. These are factors that you think may be important to control for when doing your analysis e.g. diet when looking at the genetics of body weight. Common confounders include sex, age, diet etc. as well as population stratification.

##15. Association analysis with covariates
Plink allows you to fit covariates when doing an association test. For example, to fit the PCA results to correct for population stratification we can do:
```
plink --bfile plinkFormatted --covar plink.eigenvec --covar-number 1-20 --logistic sex hide-covar
```
We are fitting as covariates the top twenty principal components as well as sex when fitting this additive model. --logistic indicates we are doing a case control study. If the phenotype was a continuous trait, such as height or weight, we could specify --linear instead of --logistic.

**_Question 10: Adapt the R code above to read this new results file and plot the p values. You can see that now nothing exceeds our significance threshold suggesting that the previous significant hits were likely false positives._**

**_Question 11: How does fitting a dominant or recessive model change the results (hint: see https://www.cog-genomics.org/plink2/assoc#linear)?_**

**_Question 12: Have a look at the same page and rerun the standard additive analysis (the code above) but with adjusted p values also outputted. What was the smallest UNADJusted p value and what was its corresponding Benjamini-Hochberg False Discovery Rate (FDR_BH)? What does this mean?_**

##16. Permutation derived p values
Permutation derived p values involve randomly shuffling the link between genotypes and phenotypes across samples to see how unusual the observed association is. P values derived in this way can be more robust than the nominal p values obtained from the standard tests however they come at the cost of being computationally more expensive to calculate. As an example though we can calculate permutation derived p values across just a small region of the chromosome:
```
plink --bfile plinkFormatted --covar plink.eigenvec --covar-number 1-20 --logistic sex perm --chr 22 --from-mb 43.05 --to-mb 43.15
```

