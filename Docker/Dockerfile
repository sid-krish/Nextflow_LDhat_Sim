FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get upgrade --yes
RUN apt-get install samtools procps --yes

RUN conda install make pandas --yes
RUN conda install -c bioconda fastsimbac seq-gen --yes
RUN conda install -c conda-forge biopython --yes
RUN conda clean --all --yes

COPY art_illumina /usr/local/bin/
RUN chmod +x /usr/local/bin/art_illumina
ENV PATH="/usr/local/bin/art_illumina:${PATH}"

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