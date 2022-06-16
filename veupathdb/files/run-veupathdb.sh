#!/bin/sh

# this is so that we can write to stdout and stderr
# if we have cron run this as the biothings user, we can't have stdout and stderr
(
    cd /home/biothings || exit
    export PYTHONPATH="."
    # only one instance will be running at a time
    flock --verbose --nonblock /crawler.lock \
      su -c '/home/biothings/venv/bin/python3.10 ndjson.py' biothings
) > /proc/1/fd/1 2> /proc/1/fd/2