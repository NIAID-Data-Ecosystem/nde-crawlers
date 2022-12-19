#!/bin/sh

env >> /etc/environment

/home/biothings/run-dockstore.sh 

# start cron in the foreground (-f) (replacing the current process) set log level (-L) to 15
exec /usr/sbin/cron -f -L 15
