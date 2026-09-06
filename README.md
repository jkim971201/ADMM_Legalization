# ADMM_Legalization
An annonymous repository of our IEEE TCAD submission,
"An Open-Source Negotiation-based Legalization Method with Alternating Direction Method of Multipliers"


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

If you want to know how to write a command tcl file, 
please see template files in the test directory.

# How to reproduce the results in the publication

Please note that you should place each benchmark suite in the test directory as ispd2024, iccad2017, ispd2015.

```
cd test
python run_ispd24.py
python run_iccad17.py
python run_ispd15.py
```

# ISPD 2024 Benchmark Suite Global Placement Results
[Link](https://www.dropbox.com/scl/fo/k65yvm0thy5zcqnaar0e6/AL4sOiMd02TvRYVKVXLsz1E?rlkey=gkvy6813d9b6z3pkhfya3xkax&st=j55zfvji&dl=0)
