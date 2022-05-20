# syntax=docker/dockerfile:1
FROM condaforge/mambaforge
LABEL MAINTAINER="tardis.supernova.code@gmail.com"

COPY conda-linux-64.lock /tmp
RUN mamba create -n tardis --file /tmp/conda-linux-64.lock \
    && rm /tmp/conda-linux-64.lock

COPY . /tmp/tardis 
RUN conda run -n tardis pip install /tmp/tardis \
    && echo "conda activate tardis" >> ~/.bashrc \
    && rm -rf /tmp/tardis \
    && mkdir /workdir

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "tardis", "/bin/bash", "-c"]
CMD ["jupyter notebook --notebook-dir=/workdir --ip='*' --no-browser --allow-root"]
