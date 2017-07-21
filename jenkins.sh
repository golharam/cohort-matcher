#!/bin/bash

if [ ! -d python_env ]; then
  virtualenv python_env
fi
source python_env/bin/activate
#pip install -U pip setuptools coverage mock pylint nose pysam pyvcf fisher numpy

#python_env/bin/nosetests --with-coverage --cover-erase --cover-package cohort_matcher --cover-html
#pylint -f parseable cohort_matcher.py > pylint.out


