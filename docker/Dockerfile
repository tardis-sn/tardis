# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

ARG REPO_DIR=/tmp/tardis
COPY . $REPO_DIR 

WORKDIR $REPO_DIR
RUN mamba create -n tardis --file conda-linux-64.lock
RUN conda run -n tardis pip install . \
        && echo "conda activate tardis" >> ~/.bashrc

RUN conda clean --all \
        && apt-get autoremove --purge -y \
        && rm -rf /var/lib/apt/lists/* \
        && rm -rf /tmp/*

COPY ./docker/scripts/* /usr/bin
RUN chmod +x /usr/bin/entrypoint
RUN chmod +x /usr/bin/jupyter-server
WORKDIR /workspace

ENTRYPOINT ["/usr/bin/entrypoint"]
CMD ["/bin/bash"]
