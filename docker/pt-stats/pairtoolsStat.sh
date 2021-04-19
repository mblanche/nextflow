

rep1=()
rep2=()

i=0
for f in "$@"; do
    ((++i))
    
    exp=$(basename $f)
    lib=${exp%_*}
    
    header=("lib")
    value=($lib)
    frac=($lib)
    
    while read -r line; do
	f1=$(echo $line|cut -d' ' -f1)
	f2=$(echo $line|cut -d' ' -f2)
	case $f1 in
	    total | total_* | cis | trans | cis_*)
		if [[ $i -eq 1 ]]; then
		    header+=($f1)
		fi
		value+=($f2)
		frac+=($(printf "%.1f%%" "$((1000 *  $f2 / ${value[1]}))e-1"))
		;;
	esac
    done < $f
    
    if [[ $i -eq 1 ]]; then
	rep1+=($(IFS=,; echo "${header[*]}"))
	rep2+=($(IFS=,; echo "${header[*]}"))
    fi
    rep1+=($(IFS=,;echo "${value[*]}"))
    rep2+=($(IFS=,;echo "${frac[*]}"))
done

for l in ${rep1[@]};do
    echo $l
done
echo
for l in ${rep2[@]};do
    echo $l
done

