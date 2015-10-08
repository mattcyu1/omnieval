# omnieval

Build
--------------
Remapping from one representation to another:
```
cd remap
mkdir build
cd build
cmake ../
make
```
Comparing two representations:
```
cd compsph
mkdir build
cd build
cmake ../
make
```

Usage
--------------

Example 1: conversion from equirectangular (2048x1024) to equalarea (1024x512):
```
./remap -i rect -o eqar -m 1024 -b 2048 -n 512 -v 1024 rect.yuv eqar.yuv
```

Example 2: comparing an equirectangular (2048x1024) and an equalarea (1024x512) representation using using equal weighting:
```
./compsph -i rect -o eqar -m 1024 -b 2048 -n 512 -v 1024 rect.yuv eqar.yuv sphere_655362.txt
```

Example 3: comparing an equirectangular (2048x1024) and an equalarea (1024x512) representation using latitude weighted points:
```
./compsph -i rect -o eqar -m 1024 -b 2048 -n 512 -v 1024 -w latweights.txt rect.yuv eqar.yuv sphere_655362.txt
```
