# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

COPY conda-linux-64.lock docs/tardis_example.yml /tmp/

RUN mamba create -n tardis --file /tmp/conda-linux-64.lock
RUN /opt/conda/envs/tardis/bin/pip install git+https://github.com/tardis-sn/tardis
RUN echo "conda activate tardis" >> /root/.bashrc

WORKDIR /tmp

CMD ["/bin/bash"]
