FROM python:3.9-slim
RUN pip install --no-cache --upgrade pip
RUN pip install --no-cache notebook jupyterlab
RUN pip install --no-cache matplotlib tqdm numpy

ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}


COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}
