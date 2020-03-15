# use ubuntu-16.04 image
FROM continuumio/miniconda3

# clone tardis-sn from master
RUN git clone https://github.com/tardis-sn/tardis.git

# install dependencies
RUN conda env create -f tardis/tardis_env3.yml

# build tardis-sn
RUN echo "conda activate tardis" >> ~/.bashrc
RUN echo "cd tardis" >> ~/.bashrc
ENV PATH /opt/conda/envs/tardis/bin:$PATH

ENV CONDA_DEFAULT_ENV tardis
RUN echo "tardis docker-image created"
