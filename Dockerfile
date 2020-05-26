FROM ubuntu:20.04
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libssl-dev zlib1g-dev
COPY . /app
RUN make && cp ./pblat /usr/bin/pblat && rm -rf * .git .travis.yml
CMD pblat

