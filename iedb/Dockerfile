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
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    cron && \
    apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

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
COPY --from=builder --chown=biothings:biothings /build/wheels /home/biothings/wheels
RUN python -m venv /home/biothings/venv && \
    /home/biothings/venv/bin/pip install --upgrade \
    /home/biothings/wheels/pip-*.whl \
    /home/biothings/wheels/wheel-*.whl && \
    /home/biothings/venv/bin/pip install /home/biothings/wheels/*.whl && \
    rm -rf /home/biothings/wheels
ADD files .
USER root
ARG CRONTAB=prod-crontab
RUN mv docker-entrypoint.sh / && \
    mv ${CRONTAB} /etc/crontab
RUN mkdir -p /data && chown -R biothings:biothings /data
ENV PYTHONPATH="."
VOLUME [ "/data" ]
STOPSIGNAL SIGINT
ENTRYPOINT ["/docker-entrypoint.sh"]
