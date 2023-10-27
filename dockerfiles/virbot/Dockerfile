FROM condaforge/miniforge3:latest

COPY environment.yml .

RUN sed -i 's/name: virbot/name: base/' environment.yml

RUN mamba env update -f environment.yml

RUN conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

CMD ["/bin/bash"]