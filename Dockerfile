FROM continuumio/miniconda3
MAINTAINER James Jeffryes <jamesgjeffryes@gmail.com>

ENV PATH /opt/conda/bin:$PATH
ENV LANG C

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