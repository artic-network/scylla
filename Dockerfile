FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .
COPY resources/k2_db k2_db

COPY tools/RNA_virus_detector /opt/RNA_virus_detector

# if spades does not work
# COPY tools/SPAdes-3.15.5-Linux /opt/SPAdes-3.15.5-Linux

# if 'genomad download_db' does not work:
# COPY tools/genomad_db /usr/genomad_db

RUN /opt/conda/bin/mamba env create -f /environment.yml
RUN /opt/conda/bin/mamba remove -n scylla genomad
RUN /opt/conda/bin/conda install -n scylla -c conda-forge -c bioconda genomad

ENV PATH=/opt/conda/envs/scylla/bin:${PATH}:/opt/RNA_virus_detector:/opt/SPAdes-3.15.5-Linux/bin
# ENV PATH=${PATH}:/opt/SPAdes-3.15.5-Linux/bin

RUN genomad download_db /usr

CMD ["/bin/bash"]
