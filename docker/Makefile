NAME=cohort-matcher
TAG=0.1
PROXY=--build-arg http_proxy=http://proxy-server.bms.com:8080 --build-arg https_proxy=http://proxy-server.bms.com:8080

all: build push

build:
	docker build ${PROXY} -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile ../..
	docker tag $(NAME):latest $(REGISTRY):latest
	docker tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker push $(REGISTRY):$(TAG)
	docker push $(REGISTRY):latest
