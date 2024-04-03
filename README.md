# TEENA
This is the code development repository in the background of TEENA web sever (https://sun-lab.yzu.edu.cn/TEENA/analysis/).  

# Building TEENA  

From a fresh install of CentOS, the following steps should provide all the required dependencies.  

yum install wget -y  

wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh  

chmod +x Anaconda3-2023.09-0-Linux-x86_64.sh  

bash Anaconda3-2023.09-0-Linux-x86_64.sh  

vim ~/.bashrc  

export PATH=$PATH:`pwd`   

source ~/.bashrc  

The first set of tests require bedtools and homer to be in your path.  

conda install -c bioconda bedtools -y  

conda install -c bioconda homer -y  

vim ~/.bashrc  

export PATH=$PATH:`pwd`   

source ~/.bashrc  

cd homer  



Taking the human reference genome as an example:  

perl configureHomer.pl -install hg38  

If it is custom HOMER data, you need to refer to the FASTA and GTF files of the genome, and use the loadGenome.pl command to define HOMER  

loadGenome.pl -gtf test.gtf -name test -fasta test.fa -org null  

For more details about HOMER, please refer to the help page on the website: http://homer.ucsd.edu/homer/introduction/update.html  


It is recommended to have a version of Python above 3.8 to facilitate drawing.  

wget https://www.python.org/ftp/python/3.8.6/Python-3.8.6.tgz  

tar -zxvf Python-3.8.6.tgz  

yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel gcc make  

cd Python-3.8.6/  

./configure --prefix=/usr/local/python38  

make && make install  

ln -s /usr/local/python38 /usr/local/bin/python3  

ls -l /usr/local/bin/  

vim /etc/profile  

PATH=/usr/local/python27/bin:/usr/local/python38/bin:$PATH  

export PATH  

source /etc/profile  
Some dependency packages required for TEENA operation.  

python -m pip install --upgrade pip  

python -m pip install pandas  

python -m pip install matplotlib  

python -m pip install Bio  

python -m pip install scipy  

python -m pip install openpyxl  

python -m pip install seaborn  



  
# Usage  
TEENA has 12 optional arguments:  

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
  


# Example  

python teena.py -q GATA3_hg38.bed -d hg38.repbase.bed -ch hg38.chrom.sizes -fa hg38.fa -a Homo_sapiens.GRCh38.110.gtf -o test1


  
# Citing  

Please cite this paper when using TEENA for your publications.  

Yuzhuo Li, Renzhe Lyu, Shuai Chen, Yejun Wang, Ming-an Sun. TEENA: an integrated web server for transposable element enrichment analysis in various model and non-model organisms. 2024, Submitted.  

See also https://sun-lab.yzu.edu.cn/TEENA/cite/
