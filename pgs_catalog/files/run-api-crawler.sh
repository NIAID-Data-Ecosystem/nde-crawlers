#!/bin/sh

(
    cd /home/biothings || exit
    export PYTHONPATH="."
    flock --verbose --nonblock /crawler.lock \
      su -c '/home/biothings/venv/bin/python3.10 ndjson.py' biothings ||
      flock --verbose /crawler.lock -c 'echo "failed to acquire lock, skipping crawling process"'
) > /proc/1/fd/1 2> /proc/1/fd/2
