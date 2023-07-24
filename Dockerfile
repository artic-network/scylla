FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .
COPY resources/k2_db k2_db

RUN /opt/conda/bin/mamba env create -f /environment.yml

RUN /opt/conda/bin/mamba remove -n scylla genomad
RUN /opt/conda/bin/conda install -n scylla -c conda-forge -c bioconda genomad
RUN genomad download-database .

CMD ["/bin/bash"]
