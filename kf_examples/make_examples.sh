# gfortran -L../src/ -Wall -O3 d1mach.f rkf45.f random_normal.f90 toms675.f90 kf_example1.f90 -o kf_example1 -ltoms675 -llapack -lblas
# gfortran -L../src/ -Wall -O3 -fcheck=bounds toms675.f90 -o toms675_test -ltoms675 -llapack -lblas
gfortran -L../src/ -Wall -O3 random_normal.f90 toms675.f90 pasta_filter.f90 -o pasta_filter -ltoms675 -llapack -lblas
