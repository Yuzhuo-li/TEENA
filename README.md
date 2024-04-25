# TEENA  
**This is the code development repository of the TEENA web sever.**

![Static Badge](https://img.shields.io/badge/build-passing-brightgreen)  ![Static Badge](https://img.shields.io/badge/Python-3.8%2B-cornflowerblue)  <a href="https://sun-lab.yzu.edu.cn/TEENA/"> <img src="https://img.shields.io/badge/TEENA-@websever-sandybrown.svg" alt="TEENA"> </a>  

**For detailed instruction about TEENA, please visit: `https://sun-lab.yzu.edu.cn/TEENA/help/`**  

**For any question about TEENA, please contact `DX120230210@stu.yzu.edu.cn`**

----------------------------------------------------------------
# Building TEENA  

#### TEENA can run under various Linux Operating System. Here we showed all the required dependencies as tested with the CentOS Linux.  

#### The TEENA pipeline requires the following dependencies to be in your path:  

* bedtools(>2.30.0)  

* samtools(>1.10)  

* gffread(>0.12.8)  

* convert2bed(>2.4.39) is in bedops(>2.4.41)  

* Python (>3.8): pandas(>2.0.3), numpy(>1.24.4), matplotlib(>3.7.5), Bio(>1.6.2), scipy(1.10.1), openpyxl(>3.1.2), seaborn(>0.13.2) and scipy(>=1.4.1) libraries.

----------------------------------------------------------------
#### If you have any issues for installing these dependencies, please refer to their official websites:  

bedtools: `https://bedtools.readthedocs.io/en/latest/content/installation.html`  

samtools: `https://github.com/samtools/samtools`  

gffread: `https://github.com/gpertea/gffread`  

bedops：`https://bedops.readthedocs.io/en/latest/content/installation.html`  

Python3.8: `https://www.python.org/`  

----------------------------------------------------------------
# Usage  

TEENA has 12 optional arguments: 
```
python teena.py [options]  

  -q	The query bed file you want to analyze. Your query bed file only needs to include the chromosome name and its starting and ending positions in three columns(required=True).  

  -d 	The repbase annotation bed file we have already downloaded from UCSC, the first three columns are the same as the query file, but the 4th column contains the TE family names(required=True).  
  
  -uk 	Upstream from transcription initiation site(default=500, required=False).  
  
  -dk	Downstream from transcription initiation site(default=500, required=False).  
  
  -o 	Result file after running the process, and just enter the prefix name which you want, and the suffix xlsx will be automatically added(required=True).  
  
  -a	The GTF annotation file we have already downloaded from Ensembl(required=True).  
  
  -fa	The genome sequence file we have already downloaded from UCSC(required=False).  
  
  -ch	The chromosome length file we have already downloaded from UCSC(required=True).  
  
  -m	We believe that the midpoint of the query file is on the interval of the corresponding annotation file, that is, two intervals overlap; If your query file is certain broad-peak file (histone modified like H3K9me3, H3K27me3 etc.), then we consider the intersection of the comment file interval and the corresponding interval of the query file as overlap (required=False, default='True').  
  
  -n	The bed file of gap from genome sequence which we have preprocessed(required=False, default=None).  
  
  -rn	Choose whether to remove the gaps from genome sequence(required=False, default='True').  
  
  -rp	Choose whether to remove the promoter regions(required=False, default='True').
```

### Example  
------------------------------------------------------------------------------------------------
```
python ./teena.py -q GATA3_hg38.bed -d hg38.repbase.bed -ch hg38.chrom.sizes -fa hg38.fa -a Homo_sapiens.GRCh38.110.gtf -o test
```

### Run tests for animals:  
------------------------------------------------------------------------------------------------
```
## The first set of tests require bedops and bedtools to be in your path.
 
cd TEENA-master

cd ./test/animals/

bash step1_datadownload.sh

bash step2_fa.out2repbase.sh

bash step3_runteena.sh
```

### Run tests for plants:  
------------------------------------------------------------------------------------------------  
```
## The first set of tests require gffread, bedtools and samtools to be in your path.

cd TEENA-master

cd ./test/plants/

bash step1_datadownload.sh

bash step2-1_make_TE_file.sh

bash step2-2_makeTEfile_attention.sh  irgsp1.TE.bed

bash step3_create_chrom_sizes.sh  Oryza_sativa.IRGSP-1.0.dna.toplevel.fa

bash step4_gff2gtf.sh

bash step5_runteena.sh
```

# Citation  
Please cite this paper when using TEENA for your publications:  

*Yuzhuo Li, Renzhe Lyu, Shuai Chen, Yejun Wang, Ming-an Sun. TEENA: an integrated web server for transposable element enrichment analysis in various model and non-model organisms. 2024, Submitted.*  

See also `https://sun-lab.yzu.edu.cn/TEENA/cite/`  
