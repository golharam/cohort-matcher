FROM python:2.7
MAINTAINER Ryan Golhar <ryan.golhar@bms.com>

# Metadata
LABEL container.base.image="python:2:7"
LABEL software.name="cohort-matcher"
LABEL software.version="0.2"
LABEL software.description="NGS Sample matching pipeline"
LABEL software.website="https://github.com/golharam/cohort-matcher"
LABEL software.documentation="README.md"
LABEL software.license="CCL 4.0"
LABEL tags="Genomics"

# System and library dependencies
RUN apt-get -y update && apt-get clean

# Pre-reqs
RUN pip install numpy && \
    pip install argparse boto3 fisher pyvcf

# cohort-matcher
COPY compareGenotypes.py /
COPY common.py /
RUN chmod +x /compareGenotypes.py

# Application entry point
#ENTRYPOINT ["python", "/run_cohort_matcher.py"]
