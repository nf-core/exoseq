FROM continuumio/miniconda
COPY environment.yml /
RUN conda env create -f /environment.yml
ENV PATH /opt/conda/envs/nfcore-exoseq/bin:$PATH