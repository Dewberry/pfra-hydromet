# Based on https://github.com/SCiO-systems/cgspatial-notebook/blob/master/Dockerfile
FROM jupyter/datascience-notebook:notebook-6.0.3

USER root

RUN apt-get update && apt-get install software-properties-common -y
RUN add-apt-repository ppa:ubuntugis/ubuntugis-unstable && apt-get update
RUN apt-get install gdal-bin -y && apt-get install libgdal-dev -y
RUN export CPLUS_INCLUDE_PATH=/usr/include/gdal && export C_INCLUDE_PATH=/usr/include/gdal

RUN pip3 install GDAL==$(gdal-config --version | awk -F'[.]' '{print $1"."$2}') 

RUN apt-get install unrar -y && \
    apt-get install lftp -y && \
    apt-get install libproj-dev -y && \
    apt-get install libgdal-dev -y && \
    apt-get install gdal-bin -y && \
    apt-get install proj-bin -y 

RUN apt-get remove pkg-config -y

ADD . /home/jovyan/work
WORKDIR /home/jovyan/work

RUN pip3 install branca && \
    pip3 install -r requirements.txt && \ 
    jupyter contrib nbextension install --user && \
    jupyter nbextension enable toc2/main

ENV PROJ_LIB="/opt/conda/share/proj"

RUN chown jovyan:users /home/jovyan/.local/share/jupyter