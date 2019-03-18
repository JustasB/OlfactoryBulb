grep light1_peak run*/parameters.hoc > $1
grep breath_peak run*/parameters.hoc >> $1
grep net_type run*/parameters.hoc >> $1
grep "n =" run*/num_of_columns.hoc  >> $1
cat $1 | sed "s/run_/run /g" > tmp$1
cat tmp$1 | sed "s/\/p/ p/g" > tmp2$1
cat tmp2$1 | sed "s/\/n/ n/g" > tmp3$1
cat tmp3$1 | sort -n -k 2 > sorted_$1
rm tmp$1 tmp2$1 tmp3$1
mv sorted_$1 $1
