FROM continuumio/miniconda
COPY container-requirements.txt /
RUN conda config --add channels defaults \
 && conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda install -y --file /container-requirements.txt
