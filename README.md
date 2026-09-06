# ADMM_Legalization

# Dependency
- GCC
- Tested on GCC 8/9/10/11
- CUDA
  - Tested on 11.8
- Qt
  - Tested on Qt5
- Eigen3
- Flex
- Bison
- TCL 8.6

# How to run
```
./SkyLine ${command}.tcl
```

If you want to know how to write a command tcl file, please see template files in the test directory.


# How to reproduce the results in the publication

Please note that you should place each benchmark suite in the test directory as ispd2024, iccad2017, ispd2015.

```
python run_ispd24.py
python run_iccad17.py
python run_ispd15.py
```
