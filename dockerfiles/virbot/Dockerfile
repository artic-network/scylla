FROM condaforge/miniforge3:latest

COPY environment.yml .

RUN sed -i 's/name: virbot/name: base/' environment.yml

RUN mamba env update -f environment.yml

RUN git clone https://github.com/rmcolq/RNA_virus_detector.git

ENV PATH=/RNA_virus_detector/:$PATH

RUN wget https://github.com/rmcolq/RNA_virus_detector/raw/main/ref.zip \
    && unzip ref.zip -d RNA_virus_detector \
    && rm ref.zip

RUN conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

CMD ["/bin/bash"]