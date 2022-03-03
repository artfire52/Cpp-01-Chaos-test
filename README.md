# Cpp-01-Chaos-test
This is an implementatino of the 01 chaos test as presented in [**The 0-1 Test for Chaos: A review**](https://www.maths.usyd.edu.au/u/gottwald/preprints/testforchaos_MPI.pdf). An other [**article**](https://hal.archives-ouvertes.fr/hal-02388470/document) could be usefull to understand the implementation. A logistic map generator is used to provide test file in order to apply the choas test.

## Generate test file
You can generate test file with:
```py
logistic_map.py filename mu N x0

```
We remind you that the logistic map is defined as:
<img src="https://latex.codecogs.com/gif.image?\dpi{80}&space;\bg_white&space;x_{n&plus;1}&space;=&space;\mu&space;x_n(1&space;-&space;x_n),&space;~x_0\in&space;[0;&space;1]." title="\bg_white x_{n+1} = \mu x_n(1 - x_n), ~x_0\in [0; 1]." />
N is the number of *x* to generate.

example for chaotic behaviour:

```py
logistic_map.py chaotic 3.7 5000 0.5

```
example for regular behaviour:

```py
logistic_map.py regular 2.8 5000 0.5

```

# Understanding the code

The code may not be really clear about the regression method. A pdf file computaiton_of_k.pdf is provided to help.
