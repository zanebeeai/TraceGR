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