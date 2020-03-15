# use ubuntu-16.04 image
FROM continuumio/miniconda3

# clone tardis-sn from master
RUN git clone https://github.com/tardis-sn/tardis.git

# install dependencies
RUN conda env create -f tardis/tardis_env3.yml

# set up conda path and activate tardis environment
RUN echo "conda activate tardis" >> ~/.bashrc
RUN echo "cd tardis" >> ~/.bashrc
ENV PATH /opt/conda/envs/tardis/bin:$PATH

# set conda default env as tardis
ENV CONDA_DEFAULT_ENV tardis

# success-message
RUN echo "tardis docker-image created"