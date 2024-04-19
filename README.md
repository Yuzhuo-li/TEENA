# TEENA
![Static Badge](https://img.shields.io/badge/build-passing-brightgreen)  

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

#### The first set of tests require bedops, bedtools and homer to be in your path.  
------------------------------------------------------------------------------------------------
```
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2

tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2

conda install -c bioconda bedtools -y

conda install -c bioconda homer -y  

vim ~/.bashrc  

export PATH=$PATH:`pwd`   

source ~/.bashrc  
```

#### Taking the human reference genome as an example:  
------------------------------------------------------------------------------------------------
```
cd homer

perl configureHomer.pl -install hg38  
```

#### Loading custom genomes, you need to refer to the FASTA and GTF files of the genome, and use the loadGenome.pl command to define homer.    

```
loadGenome.pl -gtf test.gtf -name test -fasta test.fa -org null
```
  
#### For more details about homer, please refer to the help page on the website:  

`http://homer.ucsd.edu/homer/introduction/update.html`


#### It is recommended to have a version of Python above 3.8 to facilitate drawing.  
------------------------------------------------------------------------------------------------
```
conda search "^python$"

conda create --name py38 python=3.8 -c conda-forge -y

conda activate py38

```
#### Python3.8 can also be customized for installation.
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
python -m pip install --upgrade pip  

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


#### Test  
------------------------------------------------------------------------------------------------
```
cd TEENA-master

cd test

chmod +x *

bash step1_datadownload.sh

bash step2_fa.out2repbase.sh

bash step3_runteena.sh
```
  
# Citing  

Please cite this paper when using TEENA for your publications:  

*Yuzhuo Li, Renzhe Lyu, Shuai Chen, Yejun Wang, Ming-an Sun. TEENA: an integrated web server for transposable element enrichment analysis in various model and non-model organisms. 2024, Submitted.*  

See also `https://sun-lab.yzu.edu.cn/TEENA/cite/`
