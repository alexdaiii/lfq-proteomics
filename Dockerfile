FROM mambaorg/micromamba:debian-slim

ENV PYTHONPATH=/home/$MAMBA_USER/proteomics:${PYTHONPATH:-}
ENV USERNAME=$MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# fix R kernel
RUN sed -i 's|"R"|"/opt/conda/bin/R"|g' /opt/conda/share/jupyter/kernels/ir/kernel.json

# install GOPlot
RUN /opt/conda/bin/R -e "install.packages('GOplot', repos='http://cran.rstudio.com/', lib='/opt/conda/lib/R/library')"

# copy over the code
WORKDIR /home/$MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER proteomics /home/$MAMBA_USER/proteomics

ENTRYPOINT ["/opt/conda/bin/python", "proteomics/run_analysis.py"]