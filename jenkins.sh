#!/bin/bash

#virtualenv python_env
#source python_env/bin/activate
#pip install -U pip setuptools coverage mock pylint
python_env/bin/nosetests --with-coverage --cover-erase --cover-package cohort_matcher --cover-html
pylint -f parseable cohort_matcher.py > pylint.out


