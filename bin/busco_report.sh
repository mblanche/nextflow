#/usr/bin/env bash


for f in "$@"; do
    sp=$(echo short_summary.specific.eukaryota_odb10.test.txt | awk -F'.' '{print $(NF-1)}')
    if [ -z "$i" ]; then
	grep -A8 "Results" $f | tail -n+4| awk '{$1=$1"\t";print $0}'|cut -d$'\t' -f2|awk 'BEGIN{print "Sample"}{$1=$1;print $0}'|paste -s -d, -
	i=1
    fi
    grep -A8 "Results" $f | tail -n+4| awk '{$1=$1"\t";print $0}'|cut -d$'\t' -f1|awk -v sp=${sp} 'BEGIN{print sp}{$1=$1;print $0}'|paste -s -d, -
done
