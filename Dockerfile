FROM ubuntu:20.04
WORKDIR /data
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libssl-dev zlib1g-dev
COPY . /data
RUN make && cp ./pblat /usr/bin/pblat && rm -rf * .git .travis.yml
ENTRYPOINT ["pblat"]

