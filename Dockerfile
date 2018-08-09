FROM nfcore/base
MAINTAINER Alex Peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
LABEL authors="alexander.peltzer@qbic.uni-tuebingen.de, senthilkumar.paneerselvam@scilifelab.se, m.hoeppner@ikmb.uni-kiel.de" \
description="Docker image containing all requirements for the nfcore/ExoSeq pipeline."

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
