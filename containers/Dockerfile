FROM continuumio/miniconda
COPY environment.yml /
RUN conda config --add channels defaults \
 && conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda env create -f /environment.yml