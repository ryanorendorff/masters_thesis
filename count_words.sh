#!/bin/bash

echo "Compiling"
touch tmp_error.txt
pdflatex -halt-on-error short.tex 1>tmp_error.txt

if [ "$?" -ne "0" ]; then
  echo -e "Error occurred while compiling, spitting out previous error log."
  echo -e "----------------------------------------------------------------\n"
  tail -20 tmp_error.txt
  exit 1
fi

WORDS=`pdftotext short.pdf - | egrep -E '\w\w\w+' | iconv -f ISO-8859-15 -t UTF-8 | wc -w | sed -e "s/ //g"`
echo "Number of words: $WORDS"

#rm tmp_error.txt
