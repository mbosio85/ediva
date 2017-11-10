FROM ubuntu:14.04
MAINTAINER Mattia Bosio
LABEL version="1.0" description="This is a sample desc"

# Expose port for mysql
EXPOSE 3306
ENV TERM linux
# Useful packages to install [tabix, bwa, mysql client, java, git]
RUN  apt-get update && apt-get install -y bwa 

RUN apt-get update && \
    apt-get install -y locales && \
    apt-get install -y wget && \
    apt-get install -y default-jre && \
    apt-get install -y default-jdk && \
    apt-get install -y git && \
    apt-get install -y tabix &&\
    apt-get install -y mysql-client &&\
    apt-get install  -y fastqc && \
    apt-get -f install
       
RUN  apt-get install -y  software-properties-common python-software-properties    
# install python 2.7
RUN add-apt-repository ppa:fkrull/deadsnakes && \
    apt-get update && \
    apt-get install -y python2.7  && \
    apt-get install -y python-tk 
    
    # &&\
    #cp -r /usr/bin/python2.7 /usr/bin/python
  
# install pip
RUN apt-get install -y python-pip python-dev build-essential && \
    pip install --upgrade pip && \
    pip install --upgrade virtualenv 

# install R
RUN sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list' && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | sudo apt-key add - && \
    apt-get update && \
    apt-get -y install r-base
   
# Install samtools
RUN cd /opt && \
    wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
    tar -xjf  samtools-1.5.tar.bz2 && \
    cd samtools-1.5 && \
    ./configure --prefix=/usr/bin && \
    make && \
    make prefix=/usr/ install && \
    rm /opt/samtools-1.5.tar.bz2 
    
# Install Bedtools
RUN cd /opt && \
     wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz && \
     tar -zxvf bedtools-2.25.0.tar.gz && \
     cd bedtools2 && \
     make && \
     cp ./bin/* /usr/bin/ && \
     rm ../bedtools-2.25.0.tar.gz
     
# Install Picard 1.119
RUN cd /opt && \
    wget -O picard-tools-1.119.zip  https://downloads.sourceforge.net/project/picard/picard-tools/1.119/picard-tools-1.119.zip && \
    unzip picard-tools-1.119.zip && \
    rm picard-tools-1.119.zip 

# Install eDiVA
RUN cd /opt && \
    git clone https://github.com/mbosio85/ediva.git

# Install Annovar
RUN cd /opt && \
    wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz &&\
    tar xzvf annovar.latest.tar.gz && \
    cp /opt/ediva/Resource/ediva_summarize_annovar.pl annovar/ && \
    cp /opt/ediva/Resource/maf2annovar.pl annovar/ && \
    cp /opt/ediva/Resource/convert2annovar.pl annovar/
#&& \
    
RUN cd /opt/ && cd annovar && \
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene  hg19  && \
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene  hg19  && \
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene  hg19  && \
    cd /opt/ediva/Annotate && \
    sed -i 's+/users/GD/tools/eDiVaCommandLine/lib/Annovar+/opt/annovar/+' annotate.py 
    



## SETUP environment ##

ENV EDIVA='/opt/ediva/'
ENV PYTHON='/usr/bin/python'
ENV SAMTOOLS='/usr/local/bin/samtools'
ENV PICARD='/opt/picard-tools-1.119/'
ENV GATK='/usr/bin/GenomeAnalysisTK.jar'
ENV BEDTOOLS='/usr/bin/'
ENV BWA='/usr/bin/bwa'
ENV FASTQC='/usr/bin/fastqc/'
####################################################################################################

#                   INSTALL python packages required  and R too

####################################################################################################
RUN apt-get install -y python-pip python-dev build-essential && \
    pip install numpy && \
    pip install scipy && \
    pip install pysam && \
    pip install drmaa && \
    pip install xlrd && \
    pip install xlsxwriter

RUN apt-get install -y libfreetype6-dev && \
    pip install matplotlib  
  # from https://stackoverflow.com/questions/25865270/how-to-install-python-mysqldb-module-using-pip
  
RUN apt-get install -y python-pip python-dev libmysqlclient-dev && \
    pip install MySQL-python && \
    pip install mysql-connector==2.1.4 && \
    pip install biopython
    
 RUN   R -e "install.packages('caret', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"


# Install GATK
##
## This file needs to be in Docker context (same folder) or do it from docker  run -ti image_name
 COPY GenomeAnalysisTK.jar  /usr/bin/
## 
##

    #cd /opt &&  \
    #git clone https://github.com/mysql/mysql-connector-python.git && \
    #cd mysql-connector-python && \
    #python ./setup.py build && \
ENV DISPLAY :0
