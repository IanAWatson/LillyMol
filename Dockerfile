# builder stage
FROM python:3.11.8 AS build

RUN apt-get update && \ 
    apt-get upgrade -y

RUN apt-get install npm -y && \
    npm install -g @bazel/bazelisk && \
    apt-get install libblas-dev liblapack-dev libzmq3-dev xz-utils -y

RUN pip install pandas scipy absl-py pybind11 protobuf

COPY . ./LillyMol

WORKDIR /LillyMol/src

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux \
    BUILD_BDB=1 \
    BUILD_PYTHON=1 \
    BUILD_GO=1

# Protobuf-compiler needs to be installed before build_linux.sh is run.
RUN apt-get install -y golang protobuf-compiler 

RUN ./build_linux.sh

# Remove executables currently not being used.
RUN ./uninstall.sh

# final stage
FROM python:3.11.8-slim AS final

RUN apt-get update && \ 
    apt-get upgrade -y && \
    apt-get install build-essential libgomp1 ruby-dev protobuf-compiler -y && \
    rm -rf /var/lib/apt/lists/*

RUN gem install google-protobuf -v 3.21.12

RUN pip install pandas scipy absl-py pybind11 protobuf

COPY --from=build /LillyMol /LillyMol

ENV LILLYMOL_HOME=/LillyMol \
    BUILD_DIR=Linux

WORKDIR /LillyMol
