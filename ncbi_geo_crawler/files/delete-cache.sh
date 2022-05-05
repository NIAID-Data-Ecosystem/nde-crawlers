#!/bin/sh

# Deletes the scrapy cache of the spider
# https://github.com/scrapy/scrapy/issues/3204
(
    cd /home/biothings || exit
    # I have no clue if there is a better way of getting the spider name other than this way
    # we can manually write in the spider name if no other solution is found
    # NAME=$(grep 'name' /home/biothings/spider.py | cut -d "'" -f 2)
    # temp solution to pass as a parameter
    DIR=/cache/$1

    if [ -d $DIR ]; then
      flock --verbose --nonblock /crawler.lock \
      -c "echo 'Removing /cache/$1' ; rm -rf $DIR && echo '/cache/$1 removed'"
    fi
) > /proc/1/fd/1 2> /proc/1/fd/2