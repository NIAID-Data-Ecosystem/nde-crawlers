#!/bin/sh

# this is so that we can write to stdout and stderr
# if we have cron run this as the biothings user, we can't have stdout and stderr
(
    cd /home/biothings || exit
    
    # https://github.com/john-kurkowski/tldextract/blob/master/tldextract/cache.py#L66
    # pass optional flag for caching
    while getopts "s:" opt; do
      case $opt in
        s)
          SPIDER_NAME=$OPTARG
          export TLDEXTRACT_CACHE="/cache/$SPIDER_NAME/.cache"
          ;;
      esac
    done
    
    export PYTHONPATH="."
    export SCRAPY_SETTINGS_MODULE="settings"
    # only one instance will be running at a time
    flock --verbose --nonblock /crawler.lock \
      su -c '/home/biothings/venv/bin/scrapy runspider spider.py' biothings ||
      flock --verbose /crawler.lock -c 'echo "failed to acquire lock, skipping crawling process"'
) > /proc/1/fd/1 2> /proc/1/fd/2
