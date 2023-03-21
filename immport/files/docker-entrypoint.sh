#!/bin/sh

env >> /etc/environment
# start cron in the foreground (replacing the current process)
exec /usr/sbin/cron -f -L 15