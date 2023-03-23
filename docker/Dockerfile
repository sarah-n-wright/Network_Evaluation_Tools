# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/neteval
COPY ./ /tmp/neteval/
RUN pip install /tmp/neteval

RUN rm -rf /tmp/neteval

ENTRYPOINT ["/opt/conda/bin/netevalcmd.py"]

CMD ["--help"]
