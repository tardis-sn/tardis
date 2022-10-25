# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

ENV REPO_DIR=/tmp/tardis
COPY . $REPO_DIR

WORKDIR $REPO_DIR
RUN mamba create -n tardis --file conda-linux-64.lock
RUN conda run -n tardis pip install . \
    && echo "conda activate tardis" >> ~/.bashrc

RUN conda clean --all \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

COPY ./docker/scripts/jupyter-server.sh /usr/local/bin
RUN chmod +x /usr/local/bin/jupyter-server.sh

COPY ./docker/scripts/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

WORKDIR /workspace

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash"]
