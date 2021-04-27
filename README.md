# Vector Transport Free LBFGS on Riemannian Manifolds

The code for vector transport free LBFGS quasi-Newton's optimization on the Symmetric Positive Definite Riemannian manifolds.
This project reduces vector transport and Riemannian metric to identity and Euclidean inner product, respectively. 

## Additional notes (related codes)

The code of cautious LBFGS in manopt toolbox and its cautious RLBFGS solver: 
- https://github.com/NicolasBoumal/manopt/blob/master/manopt/solvers/bfgs/rlbfgs.m

The code of Mixest toolbox (it includes its RLBFGS solver): 
- http://visionlab.ut.ac.ir/resources/riemmix.zip
- http://visionlab.ut.ac.ir/resources.html

The code of another version of Mixest toolbox and its RLBFGS solver:
- https://github.com/utvisionlab/mixest
- https://github.com/utvisionlab/mixest/tree/master/mixest/auxiliary/manopt_solvers/lbfgs
- http://visionlab.ut.ac.ir/mixest/

The code of ROPTLIB toolbox and its RBFGS and RLBFGS solvers:
- https://github.com/whuang08/ROPTLIB
- https://github.com/whuang08/ROPTLIB/tree/master/Solvers
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/RBFGS.cpp
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/RBFGS.h
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/LRBFGS.h
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/LRBFGS.cpp
