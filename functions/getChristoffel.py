""" 
MIT License

Copyright (c) 2023 Arun Kulathingal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

""" 

# reference code: https://github.com/kul-arun/christoffel-symbols-calculator/ 

import numpy as np
import sympy as sp
from itertools import product

def getChristoffel(
    t: sp.Symbol, x: sp.Symbol, y: sp.Symbol, z: sp.Symbol, metric: sp.Matrix
) -> np.ndarray:

  metric_inv = metric.inv()
  print("start computation")
  n=4
  X = [t, x, y, z]
  # computing the symbols using the metric equation
  # Create array to store the computed christoffel symbols.
  christoffel_symbols = np.zeros(shape=n, dtype='object')
  simple = False
  for i in range(n):
      dummy_matrix = sp.Matrix.zeros(n,n)
      for (j,k,l) in product(range(n), repeat=3):
          dummy_matrix[j,k] += (
              sp.Rational(1/2)*metric_inv[i,l] * (sp.diff(metric[l,j],X[k])
              +sp.diff(metric[l,k],X[j]) - sp.diff(metric[j,k],X[l]))
          )
          print(f"done connection j: {j} k: {k} l: {l}")
      christoffel_symbols[i] = sp.simplify(dummy_matrix) if simple else dummy_matrix
  C = christoffel_symbols

  return C
