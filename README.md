# TEENA
![Static Badge](https://img.shields.io/badge/build-passing-brightgreen)  ![Static Badge](https://img.shields.io/badge/Conda-build-brightgreen)  
![Static Badge](https://img.shields.io/badge/Python-3.8%2B-cornflowerblue)

This is the code development repository of TEENA web sever：  
`https://sun-lab.yzu.edu.cn/TEENA/`

# Building TEENA  

#### TEENA can run in various Linux, here we showed how to install teena from CentOS, the following steps should provide all the required dependencies.  
------------------------------------------------------------------------------------------------  
```
yum install wget -y

wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh

chmod +x Anaconda3-2023.09-0-Linux-x86_64.sh

bash Anaconda3-2023.09-0-Linux-x86_64.sh

vim ~/.bashrc  

export PATH=$PATH:`pwd`   

source ~/.bashrc
```

#### The first set of tests require bedops, bedtools, samtools and gffread to be in your path.  
------------------------------------------------------------------------------------------------
```
conda install -c bioconda bedtools -y

conda install -c bioconda samtools openssl=1.0
```

```
wget https://github.com/gpertea/gffread/archive/refs/heads/master.zip

unzip gffread-master.zip

cd gffread-master

make release

./gffread
```

```
mkdir bedops

wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2

tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2

cd /bedops/bin

./rmsk2bed

# To prevent the existence of bugs, it is recommended to use anabsolute path of gffread and rmsk2bed (for example: /home/download/bedops/bin/rmsk2bed) 

vim ~/.bashrc  

export PATH=$PATH:`pwd`   

source ~/.bashrc
``` 

If there are any issues with your Bedops and gffread installation, please refer to the installation guide on the official website:  
`https://bedops.readthedocs.io/en/latest/content/installation.html`  
`https://github.com/gpertea/gffread`

#### It is recommended to have a version of Python above 3.8 to facilitate drawing.  
------------------------------------------------------------------------------------------------
```
conda search "^python$"

conda create --name py38 python=3.8 -c conda-forge -y

conda activate py38
```
#### Python3.8 can also be customized for installation:
------------------------------------------------------------------------------------------------
```
wget https://www.python.org/ftp/python/3.8.6/Python-3.8.6.tgz

tar -zxvf Python-3.8.6.tgz

yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel gcc make  

cd Python-3.8.6/  

./configure --prefix=`your path of Python3.8`  (for example: /home/download/python/python3.8/)

make && make install
```

#### To test whether python3.8 installation was successful.
------------------------------------------------------------------------------------------------
```
cd /home/download/python/python3.8/bin

ls

./python3.8

quit()
```

#### To establish a soft link for Python 3.8, first check the situation in the/usr/bin/directory before the establishment:
------------------------------------------------------------------------------------------------
```
ln -s /home/download/python/python3.8/bin/python3.8 /usr/bin/python3

ln -s /home/download/python/python3.8/bin/pip3.8 /usr/bin/pip3

ll /usr/local/bin
```

#### Check for availability:
------------------------------------------------------------------------------------------------
```
python3

quit()

vi ~/.bashrc

export PATH=$PATH:`pwd`  

source ~/.bashrc
```

#### Some dependency packages required for TEENA operation.  
------------------------------------------------------------------------------------------------
```
python -m pip install pandas  

python -m pip install matplotlib  

python -m pip install Bio  

python -m pip install scipy  

python -m pip install openpyxl  

python -m pip install seaborn  
```

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
  
  -m	We believe that the midpoint of the query file is on the interval of the corresponding annotation file, that is, two intervals overlap; If your query file is certain broad-peak file (histone modified like H3K9me3, H3K27me3 etc.), then we consider the intersection of the comment file interval and the corresponding interval of the query file as overlap(required=False, default='True').  
  
  -n	The bed file of gap from genome sequence which we have preprocessed(required=False, default=None).  
  
  -rn	Choose whether to remove the gaps from genome sequence(required=False, default='True').  
  
  -rp	Choose whether to remove the promoter regions(required=False, default='True').
```

#### Example  
------------------------------------------------------------------------------------------------
```
python ./teena.py -q GATA3_hg38.bed -d hg38.repbase.bed -ch hg38.chrom.sizes -fa hg38.fa -a Homo_sapiens.GRCh38.110.gtf -o test
```
For questions and discussion about TEENA please visit/join the mailing list: 
`https://sun-lab.yzu.edu.cn/TEENA/help/`


#### Run tests for animals:  
#### The first set of tests require bedops and bedtools to be in your path.
------------------------------------------------------------------------------------------------
```
cd TEENA-master

chmod +x *

cd ./test/animals/

chmod +x *

bash step1_datadownload.sh

bash step2_fa.out2repbase.sh

bash step3_runteena.sh
```

#### Run tests for plants:  
#### The first set of tests require gffread, bedtools and samtools to be in your path.
------------------------------------------------------------------------------------------------  
```
cd TEENA-master

chmod +x *

cd ./test/plants/

chmod +x *

bash step1_datadownload.sh

bash step2-1_make_TE_file.sh

bash step2-2_makeTEfile_attention.sh  irgsp1.TE.bed

bash step3_create_chrom_sizes.sh  Oryza_sativa.IRGSP-1.0.dna.toplevel.fa

bash step4_gff2gtf.sh

bash step5_runteena.sh
```

# Citing  

Please cite this paper when using TEENA for your publications:  

*Yuzhuo Li, Renzhe Lyu, Shuai Chen, Yejun Wang, Ming-an Sun. TEENA: an integrated web server for transposable element enrichment analysis in various model and non-model organisms. 2024, Submitted.*  

See also `https://sun-lab.yzu.edu.cn/TEENA/cite/`  

# Contact us
If you have any further questions, please feel free to contact `DX120230210@stu.yzu.edu.cn`
