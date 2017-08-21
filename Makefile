NAME=cohort-matcher
REGISTRY=483421617021.dkr.ecr.us-east-1.amazonaws.com
VERSION=0.1
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
