

VERSION='1.0'
repo='dovetailg'
[ ! -z "$VERSION" ] && v=":$VERSION"

for d in $(find ${PWD}  -maxdepth 1 -type d);do
    dir=$(basename $d)
    img="${repo}/${dir}${v}"
    
    if [ -f "$d/Dockerfile" ]; then
	echo Building $dir
	echo "--->>>"
	docker build -t  $img $d \
	    && docker push $img
    else
	echo $dir HAS NO DOCKERFILE, Check
    fi
done
