# Session 1: Variant calling with GATK

## Dataset
This tutorial will take you through the process of calling variants with GATK using cattle whole genome sequencing data. There are 2 bam files (aligned with BWA), each corresponding to a different cow. These contain all the reads mapping to an approximately 4Mb region centred around the cow Leptin gene on chromosome 4 (you can view the region here: http://www.ensembl.org/Bos_taurus/Location/Overview?r=4%3A91249874-95266624).

Programs required:
*	R (https://cran.r-project.org/)
*	samtools (http://www.htslib.org/doc/samtools.html)
* Picard (https://broadinstitute.github.io/picard/)
*	GATK (https://software.broadinstitute.org/gatk/)

These are already installed on the cluster. To make them available type:
```
module load samtools/1.3.1
module load picard/2.8.2
module load gatk/3.7.0
module load R/3.3.2
```

## 1. Preparing the reference genome for use with GATK
Earlier we showed you how to map your reads against a reference genome using BWA. To run GATK it is first necessary to create partner files to the fasta format reference genome you used to do the alignment. We have provided the reference genome used to do these alignments (Bos_taurus.UMD3.1.dna.toplevel.fa). You need to generate two partner files from this; a fasta index file and a sequence dictionary file.
Use samtools to generate the fasta index file:
```
samtools faidx Bos_taurus.UMD3.1.dna.toplevel.fa
```
Use Picard to generate the sequence dictionary:
```
picard CreateSequenceDictionary R=Bos_taurus.UMD3.1.dna.toplevel.fa \
 O=Bos_taurus.UMD3.1.dna.toplevel.dict
```
## 2. Sorting and indexing the bam file
First sort and index your bam files using samtools. This allows downstream programs to access its contents more quickly.
```
samtools sort cow1.bam > cow1_sorted.bam
samtools index cow1_sorted.bam
```
## 3. Mark duplicates
Use Picard to flag potential duplicate reads.
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
BQSR requires a set of known variants to filter against. We can download this in VCF format from dbSNP using the wget program. As these bams only contain reads mapping to a region of chromosome 4 we can just download known variants on this chromosome.
```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/cow_9913/VCF/vcf_chr_4.vcf.gz*
```
As this file is compressed we need to uncompress it first
```
gunzip vcf_chr_4.vcf.gz
```
We can then use this file of known sites in BQSR.
```
GenomeAnalysisTK -T BaseRecalibrator \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -I cow1_dupMarked.bam -knownSites vcf_chr_4.vcf \
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
As this bam file only contains reads mapped to a small region of the cow genome, to speed things up we are going to restrict variant calling to just this region. To do this we first need to make a target file, defining the region of interest. Using a text editor (e.g. nano) create a new file, targetRegions.bed, containing the coordinates of the region as three tab separated columns. Chromosome, start position, then end coordinate i.e.:

4 91249874  95266624 

We can then use this targetRegion.bed with the –L flag to restrict calling to this region. If you wanted to call variants across the whole genome you would just leave this out.
```
GenomeAnalysisTK -T HaplotypeCaller \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -L targetRegion.bed -I cow1_recal.bam \
 --emitRefConfidence GVCF --dbsnp vcf_chr_4.vcf \
 -o cow1.g.vcf
```
**IMPORTANT. The steps 2-5 above would normally be run for each of your samples to generate one g.vcf file per sample. However to save time we have prepared the g.vcf file for the remaining sample for use with the next step.**
## 6. Raw variant calls
Run GenotypeGVCFs over the 2 g.vcf files to get a list of raw (unfiltered) variants. Each of the separate g.vcf files is specified with the --variant flag.
```
GenomeAnalysisTK -T GenotypeGVCFs \
 -R Bos_taurus.UMD3.1.dna.toplevel.fa \
 -L targetRegion.bed -I cow1_recal.bam \
 --variant cow1.g.vcf \
 --variant cow2.g.vcf \
 -o variants.vcf
```
##7.	Variant Quality Score Recalibration (VQSR)
To build the SNP recaibration model run the following command. 
```
java -jar GenomeAnalysisTK.jar \ 
    -T VariantRecalibrator \ 
    -R Bos_taurus.UMD3.1.dna.toplevel.fa \ 
    -input variants.vcf \ 
    -resource:bov1000G,known=false,training=true,truth=true,prior=12.0 bov1000G.vcf \ 
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 vcf_chr_4.vcf \ 
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
The bov1000G.vcf is a list of variants generated by the 1000 bulls sequencing consortium. We are assuming here that they are good quality variants.  We could also, for example, use the set of variants on the Illumina bovine HD array as an even better truth set. The dbSNP vcf we downloaded earlier is expected to contain many false positive variants and so is not used to train the model.

**_Question 2: How many novel variants were in the top quality tranche?_**

**_Question 3: How do the Transition/Transversion ratios differ between known and novel variants in each tranche?_**

**IMPORTANT. The command above just runs the recalibration model on SNPs. To run it on the indel calls you would normally need to change the -mode flag and provide appropriate training and truth datasets. However for the tutorial we will just continue with the SNPs only.**

Next we modify the vcf file containing our variants so that the top quality variants are indicated with a PASS label.

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
We will use the variants_VQSR.vcf from the previous session.

Programs required:
*	VEP (www.ensembl.org/vep)
*	VCFtools (https://vcftools.github.io/index.html)

These are already installed on the cluster. To make them available type:
```
module load vep/87
module load vcftools/0.1.14
```
##8.	Variant annotation with VEP
To run VEP we specify the location of the VEP local cache that has been installed at /export/data/bio/vep/. By specifying the --merged flag we can get annotations for Ensembl and RefSeq genes.
```
variant_effect_predictor.pl --cache \
 --dir_cache /export/data/bio/vep/ --species bos_taurus \
 --vcf -i variants_VQSR.vcf -o variants_VEP.vcf --merged
```
Have a look at the html output file.

**_Question 4: How many missense variants were identified?_**

**_Question 5. Have a look at the manual page (http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html) and work out how to add SIFT scores to the output vcf. You may need to delete the previous output file if using the same name. What proportion of the missense variants are predicted to be deleterious with higher confidence?_**

##9. Variant filtering with VEP
Use the filter_vep.pl script to just extract the missense variants from the VEP annotated file.
```
filter_vep.pl -I variants_VEP.vcf -format vcf -o missense_only.vcf \
 -filter “Consequence matches missense”
```
**_Question 6: Get a count of the number of coding variants in the RBM28 gene (hint take a look here http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html)_**

Ensembl also has a REST API. This allows variant annotation to be incorporated into your own programs/scripts. You can also use this to quickly get the annotation of individual variants on the command line:
```
wget -q --header='Content-type:application/json' 'https://rest.ensembl.org/vep/cow/region/4:94947404-94947404:1/A?' -O -
```
This returns the consequence of a change to an A allele at position 4:94947404 in the cow genome. If you specify a filename after the -O parameter the results, in json format, will be output to this file. Specifying "-" on the other hand leads the results to simply be printed to screen. 

##10. Variant filtering and metrics with VCFtools
Use VCFtools to get the allele frequency of each variant:
```
vcftools --vcf variants_VEP.vcf --freq --out freqs
```
To get counts of each type of change e.g. (A<->C) and the total number of transitions and transversions
```
vcftools --vcf variants_VEP.vcf --TsTv-summary --out TsTv
```
**_Question 7: Which base changes (e.g. A<->C, A<->G etc) are observed most often in this data? Why might these be most common (hint: https://en.wikipedia.org/wiki/Deamination)_**

**_Question 8: How can you output a new vcf file with only the genotypes for the first cow with ID 2637 (https://vcftools.github.io/man_latest.html)?_** 
