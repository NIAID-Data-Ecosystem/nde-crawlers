#!/bin/bash

cd /home/biothings/nde-hub
# replace PID 1 process with Python so it will receive the SIGINT when
# people run `docker stop nde-hub`
PYTHONPATH=/home/biothings/nde-hub exec /home/biothings/venv/bin/python /home/biothings/nde-hub/bin/hub.py