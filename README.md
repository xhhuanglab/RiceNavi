## RiceNavi



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

CMD: `perl pre_mapping_to_bam_pipeline.pl RiceNavi-QTNpick.cfg fq1 fq2 SamplePrefix `   

In the 'RiceNavi-QTNpick.cfg' file, please edit the PATH to the required softwares.

After BAM file is available, further steps could be performed.



**Step1: Genotyping causal variant sites for User's sample.**  

CMD: `perl step1_Indiv_CausalVar_genotyping.pl RiceNavi-QTNpick.cfg <Sample_BAM> <SamplePrefix>`   

Then the genotype file (`SamplePrefix.RiceNavi_Causal_Var.site.geno`) of the User's sample will be generated at the folder `SamplePrefix`



**Step2: Find rice accessions in our QTNlib with different allele of each causative site.**  

CMD:`perl step2_samples_with_diff_CausalVar.pl <SamplePrefix.RiceNavi_Causal_Var.site.geno>`   
Output: 'SamplePrefix.RiceNavi_Causal_Var.site.geno.samples'  
The output format is like:  

| GeneName | Chr  | Pos_7.0 | Method_Genotyping | Ref_geno | Alt_geno | Alt_Allele_Func                | HuangHuaZhan | No. samples with different alleles | QTNlib samples with different alleles |
| -------- | ---- | ------- | ----------------- | -------- | -------- | ------------------------------ | ------------ | ---------------------------------- | ------------------------------------- |
| tms5     | Chr2 | 6397412 | GATK4             | C        | A        | thermosensitive male sterility | 0\|0         | 1                                  | A474                                  |
| GW2      | Chr2 | 8117283 | GATK4             | CA       | C        | larger grain width and weight  | 0\|0         | 3                                  | A104\|A214\|A450                      |





**Step3:  Users pick QTNs to find candidate rice donor lines in QTNlib based on specific rice breeding purpose.**  
`perl step3_pick_QTN.pl <SamplePrefix.RiceNavi_Causal_Var.site.geno.samples> <Selected_QTNlist>`  

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

The RiceNavi-Sim package is implemented taking advantage of the [PedigreeSim](https://www.wur.nl/en/show/Software-PedigreeSim.htm) software. The PedigreeSim software can simulate the genotype of the offspring if the genotypes of the parents and the genetic map are given. With the constructed rice genetic map, the genotype matrix of different generations (F1 to BCnF1) for breeding population can be simulated by RiceNavi-Sim. During each generation, RiceNavi-Sim adopted `Rice-SampleSelect` package to select the best candidate as parental lines for next generation. The simulation time can be set by users. After all simulations are performed, the likelihood can be estimated. In each generation, the likelihood was calculated based on the percentage of simulations that have the ‘ideal’ individuals with only heterozygous genotypes in the regions covering selected gene(s).

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
`perl RiceNavi-Sim_run_scripts.pl <RiceNavi-Sim.cfg> <Selected_Genes.loci>  `  

The simulation outputs will be stored in 'simulation_dir'  
After simulation process is finished, `cd` into directory `simulation_dir`, and run script `stat_likelihood.pl`  .

CMD: `perl Stat_Sim_likelihood.pl <Target_size> <backcrossing_times>`   
The <Target_size> is the size (Mb) of genomic region that covers the selected gene  
if <Target_size> is set to 2, the output files are: `stat_simulation.2M` & `stat_simulation.2M.Likelihood`  





############################################################

############################################################

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

The Rice-SampleSelect package can select suitable genotypes to facilitate the breeding process. The input file for this package is a genotyping matrix, where each column represents samples, while each row is the binned genotype (e.g. 0.3 Mb per bin). The genotyping matrix is generated by the genotyping pipeline [SEG-map](http://db.ncgr.ac.cn/SEG/) for skim genome sequencing.  The Rice-SampleSelect package can output the  summarized genotype characteristics for each individual of that population, such as No. recombination breakpoints, heterozygosity across the whole genome, the No. heterozygous genomic blocks, the size of the heterozygous regions covering the targeted genes, and etc. The samples with heterozygous genotypes on target genes are ranked according to the whole genome heterozygosity level.

Usage: 
`perl RiceNavi-SampleSelect.pl <Genotyping Matrix> <Genelist> <OutPrefix> > OutFile`

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

Note that the 1st column is the combination of chromosome ID and 

(2) \<Genelist>  
The locations (Rice MSUv7) of genes selected by Users.  
The format is like:  

| GeneID         | GeneName | Chrom | Start   | End     |
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

The output will  be like:

|   Indiv   | Selected_genes | Breakpoint_Count | Het_regions_count | Heterozygosity_Percent | Homo_donor_Pct |    DTH8(8\|4.33)    |  Ghd7.1(7\|29.62)  | Y_or_N |
| :-------: | :------------: | :--------------: | :---------------: | :--------------------: | :------------: | :-----------------: | :----------------: | :----: |
| BC1F1_85  |     2 \| 2     |        14        |        12         |        0.355251        |       0        |  (25.2)8_0\|8_25.2  | (1.2)7_28.8\|7_30  |   Y    |
| BC1F1_91  |     2 \| 2     |        16        |        14         |        0.383562        |       0        |  (28.2)8_0\|8_28.2  | (0.6)7_29.4\|7_30  |   Y    |
| BC1F1_151 |     2 \| 2     |        12        |        10         |        0.391781        |       0        | (18.9)8_1.8\|8_20.7 | (1.8)7_28.2\|7_30  |   Y    |
| BC1F1_29  |     2 \| 2     |        15        |        12         |        0.414612        |       0        | (21.6)8_3.6\|8_25.2 | (27.9)7_2.1\|7_30  |   Y    |
|  BC1F1_2  |     2 \| 2     |        19        |        16         |        0.421918        |       0        | (25.2)8_2.1\|8_27.3 | (24.6)7_5.4\|7_30  |   Y    |
| BC1F1_140 |     2 \| 2     |        16        |        13         |        0.437443        |       0        |  (26.1)8_0\|8_26.1  | (8.4)7_21.6\|7_30  |   Y    |
| BC1F1_191 |     2 \| 2     |        21        |        17         |        0.441096        |       0        |   (4.5)8_0\|8_4.5   | (1.5)7_28.5\|7_30  |   Y    |
| BC1F1_77  |     2 \| 2     |        20        |        16         |        0.483105        |       0        |  (28.2)8_0\|8_28.2  | (15.9)7_14.1\|7_30 |   Y    |
| BC1F1_69  |     2 \| 2     |        17        |        13         |        0.485845        |       0        | (22.8)8_3.6\|8_26.4 | (27.6)7_2.4\|7_30  |   Y    |
| BC1F1_149 |     2 \| 2     |        19        |        15         |        0.485845        |       0        |  (3.9)8_0.9\|8_4.8  | (2.4)7_27.6\|7_30  |   Y    |

