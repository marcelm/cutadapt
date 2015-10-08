FROM debian:jessie

RUN apt-get update -y && DEBIAN_FRONTEND=noninteractive apt-get install -y \
  python2.7-dev \
  cython

ADD . /cutadapt/

RUN cd /cutadapt/ && python setup.py install && python setup.py build_ext -i

ENTRYPOINT ["/cutadapt/bin/cutadapt"]
CMD ["--help"]

# git clone https://github.com/marcelm/cutadapt.git
# cd cutadapt
# docker build -t marcelm/cutadapt:latest .

