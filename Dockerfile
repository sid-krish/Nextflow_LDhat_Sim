FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get dist-upgrade --yes
RUN apt-get install build-essential procps --yes

RUN conda update conda
RUN conda install -c bioconda -c conda-forge pandas numpy seaborn python scipy future biopython pexpect seq-gen 'msprime>=1.1.1' gsl --yes
RUN conda clean --all --yes

WORKDIR /home
RUN git clone https://github.com/auton1/LDhat.git
WORKDIR /home/LDhat
RUN make
RUN mv complete ldhat_complete
RUN mv convert pairwise interval rhomap stat ldhat_complete lkgen fin /usr/local/bin/
WORKDIR /
RUN rm -r /home/LDhat
ENV PATH="/usr/local/bin/convert:${PATH}"
ENV PATH="/usr/local/bin/pairwise:${PATH}"
ENV PATH="/usr/local/bin/interval:${PATH}"
ENV PATH="/usr/local/bin/rhomap:${PATH}"
ENV PATH="/usr/local/bin/stat:${PATH}"
ENV PATH="/usr/local/bin/ldhat_complete:${PATH}"
ENV PATH="/usr/local/bin/lkgen:${PATH}"
ENV PATH="/usr/local/bin/fin:${PATH}"