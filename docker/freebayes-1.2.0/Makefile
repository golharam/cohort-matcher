NAME=freebayes
VERSION=1.2.0
REGISTRY=483421617021.dkr.ecr.us-east-1.amazonaws.com
PROXY=

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
	docker run --rm -ti --entrypoint "/bin/bash" ${REGISTRY}/${NAME}:${VERSION}

run:
	docker run --rm -ti ${REGISTRY}/${NAME}:${VERSION} --help