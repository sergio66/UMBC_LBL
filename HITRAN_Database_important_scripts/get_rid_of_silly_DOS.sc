#!/bin/sh
cat $1 | tr -d '\r' > ~/.tmp_file.txt
mv ~/.tmp_file.txt $1
