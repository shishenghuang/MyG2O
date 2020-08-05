
# MyG2O
**Authors:** [Shi-Sheng Huang] please contact with me via shishenghuang.net@gmail.com

G2O is a popular math library for solving factor graph based Bundle Adjustment optimization problems in visual SLAM. However, the original [g2o](https://github.com/RainerKuemmerle/g2o) library doesn't contain factors like line-, plane- reprojection factors. In this library, the line-, plane- reprojection factors are integrated, you can try it for free. See details in the doc/factors.pdf


# Prerequisites
We have tested the library in **16.04**, but it should be easy to compile in other platforms. A powerful computer (e.g. i7) will ensure real-time performance and provide more stable and accurate results.

## C++11 or C++0x Compiler
We use the new thread and chrono functionalities of C++11.

## Eigen3
Required Download and install instructions can be found at: http://eigen.tuxfamily.org. **Required at least 3.1.0**.


# Compile
```
mkdir build
cd build
cmake ..
make -j4
```
