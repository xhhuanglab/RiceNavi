### RiceNavi ###

### General Introduction

Advances in functional studies and extensive allelic variation in agronomically important genes serve as the basis of molecular breeding in rice. However, the power of molecular breeding is currently limited due to lack of effective integration of genetic findings and bioinformatics methods. We have constructed a comprehensive map of rice quantitative trait nucleotides (QTNs), which genotyped 404 rice accessions at collected rice causative variant sites for 225 QTGs. RiceNavi, was developed for QTN pyramiding and breeding route optimization with the integration of the constructed rice QTN map. The RiceNavi includes three packages, namely RiceNavi-QTNpick, RiceNavi-Sim and RiceNavi-SampleSelect. The detailed descriptions and usage of the three packages are described below.



* [RiceNavi-QTNpick](#RiceNavi-QTNpick)

* [RiceNavi-Sim](#RiceNavi-Sim)

* [RiceNavi-SampleSelect](#RiceNavi-SampleSelect)   

  

For ease of use, we also constructed a web-based version of RiceNavi:   

http://www.xhhuanglab.cn/tool/RiceNavi.html (supporting most browsers including Chrome, Firefox, and Safari, but not IE).   

---





### RiceNavi-QTNpick

Taking adavantage of the integrated QTN map, the RiceNavi-QTNpick package can take mapping BAM file of user's rice sample (receptor rice sample) as input, and call the genotype of the sample at the causative sites. The genotype of  the user's sample is compared to those of the samples in the QTN library, after which the QTN library samples harboring the alterative allele for each QTN could be provided. Users can pick the beneficial QTN(s). After picking the QTN(s), the donor sample list will be provided.

Before running RiceNavi-QTNpick package, the BAM file of the rice sample is need to be genenerated by User.
We here provide a script `pre_mapping_to_bam_pipeline.pl` (in the RiceNavi-QTNpick directory) for users to generate realigned.bam with Paired-end fastq files as  inputs.  Please use the rice genome [MSU v7](http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/) as the reference. 

Generate `realigned.bam` with User sample's paired-end fastq files as  inputs. 

CMD: `perl /Path_to_RiceNavi/RiceNavi-QTNpick/pre_mapping_to_bam_pipeline.pl RiceNavi-QTNpick.cfg fq1 fq2 SamplePrefix `   

Before running, please install the required softwares  ([bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [samtools](http://www.htslib.org/download/), [sambamba](https://lomereiter.github.io/sambamba/), [GATK3](https://anaconda.org/bioconda/gatk/files), [GATK4](https://github.com/broadinstitute/gatk/releases), [Manta](https://github.com/Illumina/manta/releases), [bam2fastq](https://github.com/jts/bam2fastq)).  Edit the PATHs of the softwares and rice genome (MSU v7) in the config file `RiceNavi-QTNpick.cfg`.

[RiceNavi-QTNpick.cfg]

\## e.g. ##

bowtie2 = /tool/bowtie2-2.3.2/bowtie2  

samtools = /tool/samtools-1.9/bin/samtools  
sambamba_soft = /tool/sambamba_v0.6.7  
GATK4 = /tool/GATK4/share/gatk4-4.1.3.0-0/gatk-package-4.1.3.0-local.jar  
GATK3 = /tool/GATKv3.7/GenomeAnalysisTK.jar  
configManta =/tool/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py  
bam2fastq = /tool/bam2fastq-1.0.0/bam2fastq  

Ref_genome = /rice_genome/MSUv7/Rice_MSUv7.fa  

#######

After BAM file is generated, further steps could be performed.  



**Step1: Genotyping causal variant sites for User's sample.**  

CMD: `perl /Path_to_RiceNavi/RiceNavi-QTNpick/step1_Indiv_CausalVar_Calling.pl RiceNavi-QTNpick.cfg <Sample_BAM> <SamplePrefix>`   

Then the genotype file (`SamplePrefix.RiceNavi_Causal_Var.site.geno`) of the User's sample will be generated in the folder `SamplePrefix`



**Users can directly upload the sample's genotype file `SamplePrefix.RiceNavi_Causal_Var.site.geno` to our online [RiceNavi](http://www.xhhuanglab.cn/tool/RiceNavi.html) (QTNpick (User Sample) ) to pick the beneficial QTN(s) and find the donor samples. The following steps for RiceNavi-QTNpick are not needed.**





**Step2: Find rice accessions in our QTNlib with different allele of each causative site.**  

CMD:`perl /Path_to_RiceNavi/RiceNavi-QTNpick/step2_samples_with_diff_CausalVar.pl <SamplePrefix.RiceNavi_Causal_Var.site.geno>`   
Output: 'SamplePrefix.RiceNavi_Causal_Var.site.geno.samples'  
The output format is like:  

| GeneName | Chr  | Pos_7.0 | Method_Genotyping | Ref_geno | Alt_geno | Alt_Allele_Func                | HuangHuaZhan | No. samples with different alleles | QTNlib samples with different alleles |
| -------- | ---- | ------- | ----------------- | -------- | -------- | ------------------------------ | ------------ | ---------------------------------- | ------------------------------------- |
| tms5     | Chr2 | 6397412 | GATK4             | C        | A        | thermosensitive male sterility | 0\|0         | 1                                  | A474                                  |
| GW2      | Chr2 | 8117283 | GATK4             | CA       | C        | larger grain width and weight  | 0\|0         | 3                                  | A104\|A214\|A450                      |





**Step3:  Users pick QTNs to find candidate rice donor lines in QTNlib based on specific rice breeding purpose.**  
`perl /Path_to_RiceNavi/RiceNavi-QTNpick/step3_pick_QTN.pl <SamplePrefix.RiceNavi_Causal_Var.site.geno.samples> <Selected_QTNlist>`  

#for example:  

[Selected_QTNlist]  

| #GeneID | Chrom | Site     |
| ------- | ----- | -------- |
| Badh2   | Chr8  | 20382858 |
| TAC1    | Chr9  | 20731844 |
| OsSOC1  | Chr3  | 1270327  |

The candidate accession IDs will be provided in the file `Candidate.Improvement_samples`  
The samples which harbor the selected alternative alleles are listed below.  
No. candidate accessions: 10   

| SampleID | Type |
| -------- | ---- |
| A138     | TRJ  |
| A202     | TRJ  |
| A210     | BAS  |
| A336     | BAS  |
| A389     | TEJ  |
| A396     | TRJ  |
| A397     | TRJ  |
| A437     | TRJ  |
| A482     | BAS  |
| A98      | BAS  |



---



### RiceNavi-Sim 

With the constructed rice genetic map, the genotype matrix of different generations (F1 to BCnF1) for breeding population can be simulated by RiceNavi-Sim. RiceNavi-Sim is taking advantage of the  [PedigreeSim](https://github.com/PBR/pedigreeSim) software, which can simulate the genotype of the offspring if the genotypes of the parents and the genetic map are given.  During each generation, RiceNavi-Sim adopted `Rice-SampleSelect` package to select the best candidate as parental lines for next generation. The simulation time can be set by users. After all simulations are performed, the likelihood can be estimated. In each generation, the likelihood was calculated based on the percentage of simulations that have the ‘ideal’ individuals with only heterozygous genotypes in the regions covering selected gene(s).

Before running the script, edit the parameters in the config file `RiceNavi-Sim.cfg`
and file for target genes `Selected_Genes.loci` based on needs.

[RiceNavi-Sim.cfg]  
#for example: setting the population size for each generation

BC1F1 = 100  
BC2F1 = 200  
BC3F1 = 300  
BCnF1 = 300  

#set the number of Backcrossing times  

BC_times = 5  

#set the number of Simulation times  

Sim_times = 100  

[Selected_Genes.loci]  
#for example:  
LOC_Os08g07740	DTH8	Chr8	4333717	4335434  



Script Usage:  
`perl /Path_to_RiceNavi/RiceNavi-Sim/RiceNavi-Sim_run_scripts.pl <RiceNavi-Sim.cfg> <Selected_Genes.loci>  `  



The simulation outputs will be stored in the `simulation_dir`  
After simulation process is finished, `cd` into directory `simulation_dir`, and run script `stat_likelihood.pl`  .

CMD: `perl Stat_Sim_likelihood.pl <Target_size> <backcrossing_times>`   
The <Target_size> is the size (Mb) of genomic region that covers the selected gene  
if <Target_size> is set to 2, the output files are: `stat_simulation.2M` & `stat_simulation.2M.Likelihood`  



<u>Note that simulations with large population size and simulation times could be time-consuming, running on  linux server is recommended. The simulations (population size less than 500) could also be performed [online](http://www.xhhuanglab.cn/tool/RiceNavi.html).</u>



############################################################

############################################################

Since the RiceNavi-Sim package is implemented taking advantage of the [PedigreeSim](https://github.com/PBR/pedigreeSim) software. Please see the copyright notice of PedigreeSim:

This distribution includes a copy of the jsci-core library which is subject to the following conditions:

JScience - Java(TM) Tools and Libraries for the Advancement of Sciences. Copyright (C) 2006 - JScience (http://jscience.org/) All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

```
* Redistributions of source code must retain the above copyright notice
  and include this license agreemeent. 
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
```

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

############################################################

############################################################





---

### RiceNavi-SampleSelect

The Rice-SampleSelect package can select suitable genotypes to facilitate the breeding process. The input file for this package is a genotyping matrix, where each column represents samples, while each row is the binned genotype (e.g. 0.3 Mb per bin). The genotyping matrix is generated by the genotyping pipeline [SEG-map](http://db.ncgr.ac.cn/SEG/) for skim genome sequencing.  The Rice-SampleSelect package can output the  summarized genotype characteristics for each individual of the population, such as No. recombination breakpoints, heterozygosity across the whole genome, the No. heterozygous genomic blocks, the size of the heterozygous regions covering the targeted genes, and etc.. The samples with heterozygous genotypes on target genes are ranked according to the whole genome heterozygosity level.

Usage: 
`perl /Path_to_RiceNavi/RiceNavi-SampleSelect/step1_RiceNavi-SampleSelect.pl <Genotyping Matrix> <Genelist> <OutPut(Indiv_info))>`

`perl /Path_to_RiceNavi/RiceNavi-SampleSelect/step2_pick_top_candidates.pl <OutPut(Indiv_info)>`

Input format:  
(1) \<Genotyping Matrix> 

The format is as follows:

| marker  | BC1F1_1 | BC1F1_2 | BC1F1_3 | BC1F1_4 | BC1F1_5 |
| ------- | ------- | ------- | ------- | ------- | ------- |
| 1_0     | 0       | 1       | 0       | 0       | 1       |
| 1_0.6   | 0       | 1       | 0       | 0       | 1       |
| 1_0.9   | 0       | 1       | 0       | 0       | 1       |
| …       | …       | …       | …       | …       | …       |
| 12_27   | 0       | 0       | 0       | 1       | 1       |
| 12_27.3 | 0       | 0       | 0       | 1       | 1       |
| 12_27.6 | 0       | 0       | 0       | 1       | 1       |

Note that the 1st column is the combination of chromosome ID and the start location of that window bin (Mb).

(2) \<Genelist>  
The locations (Rice MSUv7) of genes selected by Users.  
The format is like:  

| #GeneID        | GeneName | Chrom | Start   | End     |
| -------------- | -------- | ----- | ------- | ------- |
| LOC_Os08g07740 | DTH8     | chr8  | 4333717 | 4335434 |
| LOC_Os05g01710 | xa5      | Chr5  | 437010  | 443270  |



Output: Summarized genotype characteristics for each individual. The individuals are ordered based on Heterozygosity Percentage.  
        including  (1) No. genes targeted by heterozygous regions  
                   (2) No. recombination breakpoints  
                   (3) No. heterozygous genomic blocks  
                   (4) Heterozygosity Percentage across the whole genome  
                   (5) Homozygous (Donor genotype) percentage across the whole genome  
                   (6) Size of the heterozygous regions covering the targeted genes  
                   (7) Whether (Y or N) all the selected genes are covered by heterozygous regions   

The output for the best 5 candidates (`OutPut.top_candidates`) will  be like:

\## Information for Top 5 best candidates ##

|   Indiv   | Selected_genes | Breakpoint_Count | Het_regions_count | Heterozygosity_Percent | Homo_donor_Pct |   DTH8(8\|4.33)   |   xa5(5\|0.44)    | Y_or_N |
| :-------: | :------------: | :--------------: | :---------------: | :--------------------: | :------------: | :---------------: | :---------------: | :----: |
| BC1F1_17  |     2 \| 2     |        11        |        10         |        0.346119        |       0        |  (5.4)8_0\|8_5.4  |  (1.5)5_0\|5_1.5  |   Y    |
| BC1F1_43  |     2 \| 2     |        21        |        16         |        0.378995        |       0        | (19.2)8_0\|8_19.2 | (29.7)5_0\|5_29.7 |   Y    |
| BC1F1_194 |     2 \| 2     |        13        |        11         |        0.379909        |       0        |  (5.4)8_0\|8_5.4  | (29.7)5_0\|5_29.7 |   Y    |
| BC1F1_91  |     2 \| 2     |        16        |        14         |        0.383562        |       0        | (28.2)8_0\|8_28.2 |  (2.1)5_0\|5_2.1  |   Y    |
| BC1F1_188 |     2 \| 2     |        11        |        10         |        0.39726         |       0        | (10.8)8_0\|8_10.8 | (29.7)5_0\|5_29.7 |   Y    |