FROM python:2.7

RUN apt-get -y update && \
    apt-get clean

RUN git clone --recursive https://github.com/ekg/freebayes.git && cd /freebayes && git checkout 40155b407a4bead708aaafcadf7904854c411275
WORKDIR /freebayes
RUN make 
RUN make install
WORKDIR /
RUN rm -rf /freebayes

RUN git clone https://github.com/samtools/htslib.git
WORKDIR /htslib
RUN autoheader && autoconf && ./configure && make && make install
WORKDIR /
RUN rm -rf /htslib

RUN git clone --recursive https://github.com/samtools/samtools.git
WORKDIR /samtools
RUN autoheader && autoconf -Wno-syntax && ./configure && make && make install
WORKDIR /
RUN rm -rf /samtools

ENV LD_LIBRARY_PATH /usr/lib:/usr/local/lib:${LD_LIBRARY_PATH}
