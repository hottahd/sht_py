.PHONY: make
make:
	cd sht; make

install:
	cd sht; make clean; make
	pip install .

clean:
	cd sht; make clean