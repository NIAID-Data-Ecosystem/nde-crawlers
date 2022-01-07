# syntax=docker/dockerfile:1.2
ARG PYTHON_VERSION=3.8
# This leverages a multistage build so that
# dependencies are not always downloaded from the internet, speeding up the process
# see https://docs.docker.com/develop/develop-images/multistage-build/
FROM python:${PYTHON_VERSION} AS builder
WORKDIR /build/wheels
RUN apt-get -y -qq update && \
    apt-get install -y --no-install-recommends \
      rustc \
      cargo && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
RUN /usr/local/bin/pip install --upgrade pip wheel
ADD files/nde-hub/requirements.txt .
RUN pip wheel \
    pip \
    wheel
# cryptography is used for 0) SDK dependency 1) generate SSH keys 2) install w/o Rust
RUN pip wheel cryptography
RUN pip wheel -r requirements.txt

FROM python:${PYTHON_VERSION}-slim
# BioThings Hub depends on git (the executable)
RUN apt-get -y -qq update && \
    apt-get install -y --no-install-recommends git && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
# in config.py it will read environment variables
# so that while we use containerized ES and Mongo now, we can switch to
# any other host in the future.
# Modify config.py as you like,
ENV ES_HOST=elasticsearch:9200
ENV MONGO_HOST=mongodb:27017
# create user and group
# pick your favorite number just make sure it's the same across projects
ARG BIOTHINGS_UID=454
ARG BIOTHINGS_GID=454
RUN groupadd \
      --gid ${BIOTHINGS_GID}  \
      --system \
      biothings && \
    useradd \
      -g biothings \
      --system \
      --create-home \
      --home  /home/biothings \
      --shell /bin/bash \
      --uid ${BIOTHINGS_UID} \
      biothings
USER biothings
WORKDIR /home/biothings/nde-hub
# Create virtual env and install prebuilt packages
COPY --from=builder --chown=biothings:biothings /build/wheels /home/biothings/wheels
RUN python -m venv /home/biothings/venv && \
    /home/biothings/venv/bin/pip install --upgrade \
      /home/biothings/wheels/pip-*.whl \
      /home/biothings/wheels/wheel-*.whl && \
    /home/biothings/venv/bin/pip install /home/biothings/wheels/*.whl && \
    rm -rf /home/biothings/wheels
# Copy project files (hub stuff) to container
COPY --chown=biothings:biothings ./files/nde-hub .
# create symlinks to ssh key
# because the location is hardcoded and we don't care to fix it
# commented out because symlinks are files that we can check into source control
# so it was automatically added by the COPY instruction
# RUN ln -sv /data/nde-hub/ssh_host_key bin/ssh_host_key
# create dirs and change ownership
USER root
RUN mkdir -p /data/nde-hub/{datasources,plugins,dataupload,diff,logs,release,cache,run,esbackup} && \
    chown -R biothings:biothings /data
ADD ./files/docker-entrypoint.sh /

# export persistent data
VOLUME [ "/data" ]

# Hub only listens on the API port and SSH port
EXPOSE 7080/tcp
EXPOSE 7022/tcp

# Hub only responds to SIGINT (Ctrl-C)
# but does not respond to SIGSTOP.
# Maybe someone can fix it one day
STOPSIGNAL SIGINT
# run code as biothings user
USER biothings
ENTRYPOINT [ "/docker-entrypoint.sh" ]
