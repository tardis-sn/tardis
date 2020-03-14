# use ubuntu-16.04 image
FROM ubuntu:16.04

# set up miniconda env
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# update the image to the latest packages
RUN apt-get update && apt-get upgrade -y \
	curl \
	wget \
	git

# install miniconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
RUN conda --version

# clone tardis-sn from master
RUN git clone https://github.com/tardis-sn/tardis.git

# install dependencies
RUN conda env create -f tardis/tardis_env3.yml

# build tardis-sn
RUN source activate tardis
RUN python tardis/setup.py build_ext --inplace

RUN echo "tardis docker-image created"

CMD ['cd','tardis']











