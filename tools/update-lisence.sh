#! /bin/bash

BASENAME=`basename $0`
for file in `find . -type d -name '.git' -prune -o -type f -name ${BASENAME} -prune -o -type f -print | xargs grep -l 'Copyright.*SALMON'`
do
  start_year=`git log --pretty='%ai' ${file} | tail -n 1 | cut -c 1-4`
  echo "${start_year}-`date +%Y`: ${file}"
  sed -i -e "s/Copyright.*SALMON/Copyright ${start_year}-`date +%Y` SALMON/g" ${file}
done
