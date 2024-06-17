# when running in a development environment instead of a production environment specify crontab buildarg
# docker build --build-arg CRONTAB=dev-crontab .
FROM python:3.10 AS builder
WORKDIR /build/wheels
RUN apt-get -y -qq update && \
    apt-get install -y --no-install-recommends \
    rustc \
    cargo && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*
RUN /usr/local/bin/pip install --upgrade pip wheel
ADD files/requirements.txt .
RUN pip wheel pip wheel
RUN pip wheel -r requirements.txt


FROM python:3.10-slim
# fetch & install cron
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    cron && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# setup environment
# create user and group
# pick your favorite number just make sure it's the same across projects
ARG BIOTHINGS_UID=454
ARG BIOTHINGS_GID=454
RUN groupadd \
    --gid ${BIOTHINGS_GID} \
    --system \
    biothings && \
    useradd \
    -g biothings \
    --system \
    --create-home \
    --home /home/biothings \
    --shell /bin/bash \
    --uid ${BIOTHINGS_UID} \
    biothings
USER biothings
WORKDIR /home/biothings
# Create virtual env and install prebuilt packages
COPY --from=builder --chown=biothings:biothings /build/wheels /home/biothings/wheels
RUN python -m venv /home/biothings/venv && \
    /home/biothings/venv/bin/pip install --upgrade \
    /home/biothings/wheels/pip-*.whl \
    /home/biothings/wheels/wheel-*.whl && \
    /home/biothings/venv/bin/pip install /home/biothings/wheels/*.whl && \
    rm -rf /home/biothings/wheels
# add files
ADD files .
# crontab has to run as root
# I really don't want to bother with reinventing the wheel
USER root
# make arg to decide whether to use dev-crontab or prod-crontab, defaults as prod (--build-arg CRONTAB=dev-crontab)
ARG CRONTAB=prod-crontab
# move files
RUN mv docker-entrypoint.sh / && \
    mv ${CRONTAB} /etc/crontab
RUN mkdir -p /data && chown -R biothings:biothings /data && \
    mkdir -p /cache && chown -R biothings:biothings /cache
ENV PYTHONPATH="."
VOLUME [ "/data", "/cache" ]
STOPSIGNAL SIGINT
ENTRYPOINT ["/docker-entrypoint.sh"]
