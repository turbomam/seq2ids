#!/usr/bin/env bash

# check if value for dir_name has been passed
# if not write to sqlite_dump
if [[ "$4" != "" ]]; then
    DIR_NAME="$4"
else
    DIR_NAME="sqlite_dump"
    mkdir -p $DIR_NAME
fi

# loop over all comma separated table names
for i in $(echo $3 | tr "," "\n")
do
    echo "Exporting Table: $i"
    echo "=========================================="
    
    db-to-sqlite postgresql://$1@localhost:$2/felix $DIR_NAME/$i.db --table $i --progress
done
