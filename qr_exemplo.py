import numpy as np 
A = np.array([[0.25, 0.21875, 0.15625], [0.459556, 0.315944, 0.2585], [0.459556, 0.315944, 0.315944], [0.095976, 0.083979, 0.035991], [0.095976, 0.083979, 0.083979]]) 
q, r =  np.linalg.qr(A) 
print("Decomposition of matrix:") 
print( "q=\n", q, "\nr=\n", r)
print( "A=\n", np.dot(q,r))
