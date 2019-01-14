FROM r-base:3.5.2
MAINTAINER Ryan Golhar <ryan.golhar@bms.com>

# Metadata
LABEL container.base.image="r-base:latest"
LABEL software.name="cohort-matcher"
LABEL software.version="0.2"
LABEL software.description="NGS Sample matching pipeline"
LABEL software.website="https://github.com/golharam/cohort-matcher"
LABEL software.documentation="README.md"
LABEL software.license="CCL 4.0"
LABEL tags="Genomics"

# System and library dependencies
RUN apt-get -y update && apt-get clean
RUN apt-get install libssl-dev libcurl4-openssl-dev 

# Pre-reqs
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('argparser', 'circlize', 'canvasXpress'))"

# Application entry point
