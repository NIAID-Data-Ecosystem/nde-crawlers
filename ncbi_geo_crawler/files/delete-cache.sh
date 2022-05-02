#!/bin/sh

# Deletes the scrapy cache of the spider
# https://github.com/scrapy/scrapy/issues/3204
# I have no clue if there is a better way of getting the spider name other than this way
NAME=$(grep 'name' spider.py | cut -d "'" -f 2)
DIR=/home/jalin/cache/$NAME/*

if [ -d $DIR ]; then
  echo "Removing cache"
  rm -rf $DIR
fi