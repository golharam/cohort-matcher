#!/bin/bash
VERSION="latest"
opt="--build-arg http_proxy=http://proxy-server.bms.com:8080 --build-arg https_proxy=http://proxy-server.bms.com:8080"
docker build $opt -t 483421617021.dkr.ecr.us-east-1.amazonaws.com/cohort-matcher:$VERSION .

# Push
if [ $? -eq 0 ]; then
  aws ecr get-login --region us-east-1 | sh
  docker push 483421617021.dkr.ecr.us-east-1.amazonaws.com/cohort-matcher:$VERSION
fi
