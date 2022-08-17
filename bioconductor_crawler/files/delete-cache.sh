#!/bin/sh

# Deletes the scrapy cache of the spider
# https://github.com/scrapy/scrapy/issues/3204
(
    cd /home/biothings || exit
    # I have no clue if there is a better way of getting the spider name other than this way
    # we can manually write in the spider name if no other solution is found
    # NAME=$(grep 'name' /home/biothings/spider.py | cut -d "'" -f 2)
    # temp solution to pass spider name as a parameter
    SPIDER_NAME=$1
    DIR=/cache/$SPIDER_NAME

    if [ -d $DIR ]; then
      flock --verbose --nonblock /crawler.lock \
      su -c "echo 'Removing /cache/$SPIDER_NAME' ; rm -rf $DIR && echo '/cache/$SPIDER_NAME removed'" biothings
    fi
) > /proc/1/fd/1 2> /proc/1/fd/2
