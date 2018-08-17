FROM vanallenlab/miniconda:3.6

WORKDIR /

COPY LICENSE /

COPY requirements.txt /

RUN conda install --yes -c conda-forge --file requirements.txt

COPY datasources/ /

COPY README.md /

COPY common_variant_filter.py /