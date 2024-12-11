folders="
phasespace_surface/
"

echo "RUN TESTS"

source 00-set_path.sh

for folder in ${folders}
do
	echo "RUN ${folder}"
	pushd $folder
	make testfiles
	popd
	echo
done
