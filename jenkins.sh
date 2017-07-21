#!/bin/bash

if [ ! -d python_env ]; then
  virtualenv python_env
fi
source python_env/bin/activate
pip install -U setuptools
pip install -U nose mock pysam pyvcf numpy fisher
pip install -U coverage pylint backports.functools-lru-cache singledispatch enum34 configparser

python_env/bin/nosetests --with-coverage --cover-erase --cover-package cohort_matcher --cover-html
pylint -f parseable cohort_matcher.py > pylint.out
