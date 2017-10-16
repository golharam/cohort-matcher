FROM centos:7

# Metadata
LABEL container.base.image="centos:7"
LABEL software.name="cohort-matcher"
LABEL software.version="0.1"
LABEL software.description="NGS Sample matching pipeline"
LABEL software.website="https://github.com/golharam/cohort-matcher"
LABEL software.documentation="README.md"
LABEL software.license="CCL 4.0"
LABEL tags="Genomics"

# System and library dependencies
#yum install -y gcc gcc-c++ git wget bzip2-devel xz-devel zlib-devel ncurses-devel
yum install -y git make gcc zlib-devel bzip2-devel xz-devel gcc-c++ wget bzip2 ncurses-devel python2-pip python-devel

# Get hg19
#RUN mkdir -p /ngs/reference/hg19/chromosomes
#RUN wget -O /ngs/reference/hg19/hg19.tar.gz \
#      http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz
#RUN cd /ngs/reference/hg19/chromosomes && tar zxvf ../hg19.tar.gz && \
#    cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa \
#        chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa \
#        chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa > ../hg19.fa
#RUN rm -rf /ngs/reference/hg19/hg19.tar.gz /ngs/reference/hg19/chromosomes
# Get GRCh37
#RUN mkdir -p /ngs/reference/GRCh37
#RUN wget -O /ngs/reference/GRCh37/GRCh37.fa.gz \
#      ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz && \
#    gunzip /ngs/reference/GRCh37/GRCh37.fa.gz

# Install Freebayes
RUN git config --global url.https://github.com/.insteadOf git://github.com/ && \
    git clone --recursive git://github.com/ekg/freebayes.git
RUN cd freebayes && make && make install
RUN rm -rf /freebayes

# Install R
yum install epel-release
yum install -y R

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 && \
    tar xvjf /samtools-1.6.tar.bz2 && \
    rm -f /samtools-1.6.tar.bz2
RUN cd /samtools-1.6 && make && make install
RUN rm -rf /samtools-1.6

# cohort-matcher
git clone https://github.com/golharam/cohort-matcher.git
RUN cd cohort-matcher && \
    pip install numpy && \
    pip install -r requirements.txt
VOLUME /scratch

# Application entry point
ENTRYPOINT ["python", "/cohort-matcher/run_cohort_matcher.py"]
