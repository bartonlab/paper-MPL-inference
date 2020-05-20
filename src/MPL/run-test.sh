# g++ src/main.cpp src/inf.cpp src/io.cpp -O3 -march=native -lgslcblas -lgsl -o bin/mpl
g++ src/main.cpp src/inf-binary.cpp src/io.cpp -O3 -march=native -lgslcblas -lgsl -o bin/mpl-binary
./bin/mpl-binary -d test -i example_0_T400_ns1000_dt1.dat -o example_0_T400_ns1000_dt1_MPL.dat -g 1e3 -N 1e3 -mu 1e-3
