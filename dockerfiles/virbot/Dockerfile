FROM condaforge/miniforge3:latest

COPY environment.yml .

RUN sed -i 's/name: virbot/name: base/' environment.yml

RUN mamba env update -f environment.yml

RUN curl -O https://github.com/rmcolq/RNA_virus_detector/raw/main/ref.zip \
    && unzip ref.zip -d /opt/conda/lib/python3.11/site-packages/virbot/data/ref/ \
    && rm ref.zip

RUN conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

CMD ["/bin/bash"]