#!/bin/bash

if [ ! -d python_env ]; then
  virtualenv python_env
fi
source python_env/bin/activate
pip install -U setuptools nose mock pysam pyvcf fisher numpy coverage
#pip install -U pip setuptools  pylint  

python_env/bin/nosetests --with-coverage --cover-erase --cover-package cohort_matcher --cover-html
#pylint -f parseable cohort_matcher.py > pylint.out


