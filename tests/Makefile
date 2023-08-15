testfiles:
	bash make_all_testfiles.sh

test:
	pytest

all: clean testfiles
	pytest
	
clean:
	rm -f */outfile*
	rm -rf */__pycache__/
	rm -rf *.pytest_cache/
	rm -f */contcar*
	rm -f */*.dat
	rm -f */*.log

