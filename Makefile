FC = gfortran
# FFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -g -fcheck=all
FFLAGS = -Wall -O3
FLFLAGS =

export FC
export FFLAGS
export FLFLAGS

all:
	$(MAKE) -f Makefile --directory=src
	$(MAKE) -f Makefile --directory=drivers

test:
	$(MAKE) -f Makefile --directory=drivers test

clean:
	$(MAKE) -f Makefile clean --directory=src
	$(MAKE) -f Makefile clean --directory=drivers