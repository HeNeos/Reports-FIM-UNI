## Finite Element Analysis with Python, C++ and MATLAB 

A simple finite element analysis code for some 2D and 3D elasticity problems. The main code was written in Python and then ported to C++ and MATLAB. Jupyter Notebooks contains the Python code and provides an easy way to define your model.

The main idea of this code is to provide an easy and fast code to solve some problems in the MC516 course, so I only implemented some elements:

- Bar
- Truss
- Beam
- Frame

## Solver

The most important part of the code is the solver. Many naive implementations use the linear solver provided by some numerical libraries. For example ```np.linalg.solve()``` in Python or the ```\``` operation in MATLAB. But this is only a general way of solving a linear system and doesn't exploit the main feature in a linear FEA problem.

### Strassen's Algorithm

My first optimization is to use the Strassen Algorithm for matrix multiplication. This algorithm reduces the complexity of a naive matrix multiplication:

```c++
matrix A, B;
for(int i=0; i<A.size(); i++){
    for(int j=0; j<A[0].size(); j++){
        for(int k=0; k<B[0].size(); k++) //add values
    }
}
```
Time complexity in a naive multiplication is ![n3](https://render.githubusercontent.com/render/math?math=\mathcal{O}(n^{3})), a big complexity!

The Strassen algorithm provides a small optimization, the idea behind of this algorithm is a programming paradigm called *Divide & Conquer*. **Divide** the matrix into *symmetric* sub-matrix and multiply them, if you list the operations you can see that there are 8 of them, but 1 of them can be deduced by the other 7, so that the final complexity is reduced from ![nlog28](https://render.githubusercontent.com/render/math?math=\mathcal{O}(n^{\log_{2}8})) to ![approximation](https://render.githubusercontent.com/render/math?math=\mathcal{O}(n^{\log_{2}7})%20\approx%20\mathcal{O}(n^{2.807355})).

You can read more here:

https://en.wikipedia.org/wiki/Strassen_algorithm

### Conjugate Gradient Method

A really good solver for Finite Element Analysis is the Conjugate Gradient Method. This solver is really fast when the matrix is symmetric and many of its elements are zero.

You can read more here:

https://en.wikipedia.org/wiki/Conjugate_gradient_method

### Sparse Matrix

A nice data structure to optimizate the Conjugate Gradient Method is Sparse Matrix, this data structure reduces the complexity of space and the time to store a big matrix with many zeros. Numpy implemented this structure natively, but C++ needs an external library.

You can read more here:

https://en.wikipedia.org/wiki/Sparse_matrix

## Instructions

### Python

You can clone my repository and open the .ipynb you need. You can also download the notebook and open it with Google Colab. In some notebooks I provide a direct link to do this.

### C++

You need a library for Sparse Matrix data structure. You can find the library at https://github.com/uestla/Sparse-Matrix.
You can notice that I used ```bits/stdc++.h``` as a library, this is a very well known library used in competitive programming, but you can change this line if you don't have this library and add some as:

- ```#include <iostream>```
- ```#include <set>```
- ```#include <map>```
- ```#include <vector>```
- ```#include <algorithm>```
- ```#include <math.h>```
- ```#include <iterator>```

### MATLAB

I'm not a very good coder in MATLAB, so my code probably works slowly in it, so I strongly recommend to use C++ or Python code. But if you want to use my MATLAB version only copy the code an run with your values.
