NAME=cohort-matcher
#REGISTRY=483421617021.dkr.ecr.us-east-1.amazonaws.com/${NAME}
TAG=0.1
#PROXY=--build-arg http_proxy=http://proxy-server.bms.com:8080 --build-arg https_proxy=http://proxy-server.bms.com:8080 --build-arg ftp_proxy=http://proxy-server.bms.com:8080
PROXY=

all: build 

build:
	docker build --no-cache ${PROXY} -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile ../..
	#docker tag $(NAME):latest $(REGISTRY):latest
	#docker tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker push $(REGISTRY):$(TAG)
	docker push $(REGISTRY):latest

test:
	docker run --rm -ti --entrypoint "/bin/bash" cohort-matcher

run:
	docker run --rm -ti cohort-matcher
