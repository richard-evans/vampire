#/bin/bash

out_file=$1

for nthreads in {1..128}
do
    echo -n "timing $nthreads threads..."

    sed -i "s/gpu:num-threads=.*/gpu:num-threads=$nthreads/" input

    t=$({ time ./vampire &>/dev/null; } |& grep "real" | awk '{ print $2 }' | sed -e 's/0m//' -e 's/s//')

    echo $t

    if [ ! -z "$out_file"]
    then
        echo "$nthreads $t" >> $out_file
    fi
done
