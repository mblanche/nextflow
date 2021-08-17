#!/usr/bin/env bash

bs=$1;shift
token=$1;shift
host=$1;shift

declare -A R

res=$(bs biosample content --api-server ${host} --access-token ${token} -n ${bs} -f csv)

if [ -z ${res+x} ]; then
    echo "No biosample name $bs found"
    exit 1
fi

for l in $res; do
    [ $l == "Id,FilePath" ] && continue
    id=$(echo $l|cut -d, -f1)
    fq=$(echo $l|cut -d, -f2)


    if [[ $fq == *'_R1_'* || $fq == *'_R2_'* ]]; then
    	R[$fq]=$(echo $id ${R[$fq]}|tr ' ' ,)
    else
    	echo "found a file $fq not R1 or R2"
    	exit 1
    fi
done

for k in "${!R[@]}";do
    f=( $(echo ${R[$k]}|tr , "\n") )
    if [[ ${#f[@]} -gt 1 ]]; then
	for id in ${f[@]}; do
	    fqTS=$( date -d "$(bs file get --api-server ${host} --access-token ${token} -i $id --template='{{.DateCreated}}'|awk '{print $1,$2,$3}')" +%s )
	    if [ -z ${youngID+x} ];then
		youngID=$id
		youngTS=$fqTS
	    elif [ $fqTS -ge $youngTS ]; then
		youngID=$id
		youngTS=$fqTS
	    fi
	done
    else
	youngID=${f[0]}
    fi
    echo "$youngID,$k"
    unset youngID youngTS
done
