NAME=cohort-matcher
REGISTRY=483421617021.dkr.ecr.us-east-1.amazonaws.com
VERSION=1.0
PROXY=--build-arg http_proxy=http://proxy-server.bms.com:8080 \
      --build-arg https_proxy=http://proxy-server.bms.com:8080 \
      --build-arg ftp_proxy=http://proxy-server.bms.com:8080
#PROXY=

all: build 

build:
	docker build ${PROXY} -t ${REGISTRY}/$(NAME):$(VERSION) \
                          -t ${REGISTRY}/$(NAME):latest \
                          -f Dockerfile .

push:
	aws ecr get-login --no-include-email --region us-east-1 | sh && \
	docker push ${REGISTRY}/${NAME}:${VERSION} && \
	docker push ${REGISTRY}/${NAME}:latest

test:
	docker run --rm -ti --entrypoint "/bin/bash" ${REGISTRY}/${NAME}:latest

run:
	docker run --rm -ti ${REGISTRY}/${NAME}:latest 

hg19:
	ln -s /ng18/galaxy/reference_genomes/hg19/hg19.fa .
	ln -s /ng18/galaxy/reference_genomes/hg19/hg19.fa.fai .
	tar cvjhf hg19-cohort-matcher.tar.bz2 hg19.*
	aws s3 cp hg19-cohort-matcher.tar.bz2 s3://bmsrd-ngs-repo/reference/ --sse
	rm hg19-cohort-matcher.tar.bz2 hg19.fa hg19.fa.fai

GRCh37:
	ln -s /ng18/galaxy/reference_genomes/GRCh37-lite/GRCh37-lite.fa ./GRCh37.fa
	ln -s /ng18/galaxy/reference_genomes/GRCh37-lite/GRCh37-lite.fa.fai ./GRCh37.fa.fai
	tar cvjhf GRCh37-cohort-matcher.tar.bz2 GRCh37.* 1kg.*
	aws s3 cp GRCh37-cohort-matcher.tar.bz2 s3://bmsrd-ngs-repo/reference/ --sse
	rm GRCh37-cohort-matcher.tar.bz2 GRCh37.fa GRCh37.fa.fai
