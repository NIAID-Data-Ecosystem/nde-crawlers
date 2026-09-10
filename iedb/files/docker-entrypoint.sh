#!/bin/sh

env >> /etc/environment

exec /usr/sbin/cron -f -L 15
