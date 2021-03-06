#!/bin/bash

if [ ! -d python_env ]; then
  virtualenv python_env
  source python_env/bin/activate
  pip install -U setuptools
  pip install -U nose mock pysam pyvcf numpy
  pip install -U fisher
  pip install -U coverage pylint backports.functools-lru-cache singledispatch enum34 configparser
else
  source python_env/bin/activate
fi

python_env/bin/nosetests --with-coverage --cover-erase --cover-package cohort_matcher --cover-xml
pylint -f parseable cohort_matcher.py | tee pylint.out
