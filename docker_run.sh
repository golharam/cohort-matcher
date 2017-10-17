#!/bin/bash

aws_access_key=`grep aws_access_key ~/.aws/credentials | cut -f 2 -d '='`
aws_secret_access_key=`grep aws_secret_access_key ~/.aws/credentials | cut -f 2 -d '='`

docker run \
  --rm -ti \
  -e "http_proxy=http://proxy-server:8080" \
  -e "https_proxy=http://proxy-server:8080" \
  -e "AWS_ACCESS_KEY_ID=$aws_access_key" \
  -e "AWS_SECRET_ACCESS_KEY=$aws_secret_access_key" \
  --privileged \
  483421617021.dkr.ecr.us-east-1.amazonaws.com/cohort-matcher:latest \
  /bin/bash

