# Build using `docker build -t pblat -f Dockerfile .`
# Run using `docker run -u $(id -u):$(id -g) -v $(pwd):/data pblat` 

FROM ubuntu:18.04
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends build-essential make libssl-dev zlib1g-dev
COPY . /app
RUN cd /app && make && cp ./pblat /usr/bin/pblat
CMD pblat
