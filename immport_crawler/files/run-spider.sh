#!/bin/sh

# this is so that we can write to stdout and stderr
# if we have cron run this as the biothings user, we can't have stdout and stderr
(
    cd /home/biothings || exit
    export PYTHONPATH="."
    export SCRAPY_SETTINGS_MODULE="settings"
    su -c '/home/biothings/venv/bin/scrapy runspider spider.py' biothings
) > /proc/1/fd/1 2> /proc/1/fd/2
