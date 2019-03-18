#!/bin/bash
printf 'assembling...'

# assembly header and data
cat $1.sbgh.* > $1.header

printf "from struct import pack\nfrom os import path\nf=open('"$1".size','wb')\nf.write(pack('>q', path.getsize('"$1".header')))\nf.close()\nquit()" > wd.py
python wd.py


cat $1.sbg.* > $1.data
cat $1.size $1.header $1.data > $1.spk

rm $1.sbg.* $1.sbgh.* wd.py $1.size $1.header $1.data

printf 'done\n'
