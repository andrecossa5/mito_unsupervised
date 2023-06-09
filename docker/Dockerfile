# MI_TO image

################################################################

## Start Debian basic image: NBBB ricordarsi --platform linux/x86_64 in docker build
FROM python:3.9-bullseye

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libz-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

## Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && bash Mambaforge-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

## Clone repo and install code
RUN echo mito_utils
RUN git clone https://github.com/andrecossa5/mito_utils.git \
    && cd mito_utils \
    && mamba env update -n base -f envs/environment_mito_utils.yml \
    && mamba develop .. \
    && mamba develop ../mito_utils/ 

# RUN pipelines/genomics utils
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes pigz

RUN mamba install --yes -c bioconda samtools
RUN pip install dnaio pysam
RUN mamba install --yes -c conda-forge openjdk
RUN mamba install --yes -c bioconda snakemake fgbio 
RUN git clone https://github.com/caleblareau/maegatk.git 
RUN mamba install --yes -c bioconda bedtools bwa

# Clone tracer packages
RUN pip install shap circlify
RUN pip install pyro-ppl==1.7.0 torch==1.9.0
RUN pip install dill
COPY ./helper_functions.py /
RUN cd / && mamba develop /

## Default command to initialize the container
CMD [ "bash" ]

#################################################################

