### HowTo: C++

1. git clone <project@url>

2. build Boost

```sh
cd thirdparty/boost_1_82_0
./bootstrap prefix=.
./b2 install --with-random
```

3. compile project

```sh
g++ -o quench InitialConditions.cpp InitialConditions.h KleinGordonSol.cpp KleinGordonSol.h MonteCarlo.cpp MonteCarlo.h PARAMETERS.h main.h main.cpp -Ithirdparty/eigen-3.4.0 -Ithirdparty/boost_1_82_0/include/ -Lthirdparty/boost_1_82_0/lib -std=c++11 -lboost_random
```

4. set LD_LIBRARY_PATH

```sh
export LD_LIBRARY_PATH="thirdparty/boost_1_82_0/lib"
```
