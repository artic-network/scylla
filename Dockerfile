FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .
COPY resources/k2_db k2_db

RUN /opt/conda/bin/mamba env create -f /environment.yml

ENV PATH=/opt/conda/envs/scylla/bin:$PATH

CMD ["/bin/bash"]