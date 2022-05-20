# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

COPY conda-linux-64.lock /tmp
RUN mamba create -n tardis --file /tmp/conda-linux-64.lock
RUN echo "conda activate tardis" >> /root/.bashrc \
    && rm /tmp/conda-linux-64.lock

COPY . /tmp/tardis 
RUN /opt/conda/envs/tardis/bin/pip install /tmp/tardis \
    && rm -rf /tmp/tardis

RUN mkdir /data
WORKDIR /data

CMD ["/bin/bash"]
