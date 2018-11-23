From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Lukas Lueftinger  <lukas.lueftinger@imp.ac.at>
    DESCRIPTION Container image containing all requirements for the nf-core/exoseq pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-exoseq-1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
