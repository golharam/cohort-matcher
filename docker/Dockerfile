FROM python:2.7

# Metadata
LABEL container.base.image="python:2.7"
LABEL software.name="cohort-matcher"
LABEL software.version="0.1"
LABEL software.description="NGS Sample matching pipeline"
LABEL software.website="https://github.com/golharam/cohort-matcher"
LABEL software.documentation="README.md"
LABEL software.license="CCL 4.0"
LABEL tags="Genomics"

# System and library dependencies
RUN apt-get -y update && \
    apt-get -y install unzip && \
    apt-get clean

# Application installation
RUN wget -O /cohort-matcher.zip \
   https://github.com/golharam/cohort-matcher/archive/cohort-matcher.zip && \
   unzip /cohort-matcher.zip && mv cohort-matcher-cohort-matcher cohort-matcher && \
   rm /cohort-matcher.zip

# Other software dependencies
RUN pip install numpy
RUN pip install -r /cohort-matcher/requirements.txt

# Application entry point
RUN ln -s /cohort-matcher/cohort_matcher.py /cohort_matcher.py

ENTRYPOINT ["python", "cohort_matcher.py"]
