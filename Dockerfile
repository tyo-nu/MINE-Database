FROM continuumio/miniconda3
MAINTAINER James Jeffryes <jamesgjeffryes@gmail.com>

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

#patch to downgrade from 3.6 until supported by RDKit
RUN conda install python=3.5
# install the RDKit:
RUN conda config --add channels  https://conda.anaconda.org/rdkit
# note including jupyter in this brings in rather a lot of extra stuff
RUN conda install -y cairo \
                     cairocffi \
                     nomkl \
                     pandas \
                     pymongo \
                     rdkit

COPY minedatabase/ /mine/minedatabase/
ENV PYTHONPATH $PYTHONPATH:/mine
WORKDIR /mine/minedatabase
ENTRYPOINT ["python", "pickaxe.py"]