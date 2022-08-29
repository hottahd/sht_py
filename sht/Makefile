FC = gfortran
FFLAGS = -shared -fPIC -O3 -fopenmp -march=native
SRC_FILES := $(wildcard ./*.f90)
LIB_FILES := $(notdir $(SRC_FILES:.f90=.so))

UPDATE = update

$(UPDATE) : $(LIB_FILES)
	touch $(UPDATE)

%.so: %.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(LIB_FILES)
