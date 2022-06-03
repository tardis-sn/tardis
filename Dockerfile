# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

COPY conda-linux-64.lock /tmp
RUN mamba create -n tardis --file /tmp/conda-linux-64.lock

COPY . /tmp/repo
RUN conda run -n tardis pip install /tmp/repo \
    && echo "conda activate tardis" >> ~/.bashrc

RUN mkdir /workdir
COPY docs/quickstart.ipynb /workdir
WORKDIR /workdir

RUN conda clean --all \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/* \
    && rm /tmp/conda-linux-64.lock \
    && rm -rf /tmp/repo

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "tardis", "/bin/bash", "-c"]
CMD ["jupyter notebook --notebook-dir='/workdir' --ip='*' --no-browser --allow-root"]
