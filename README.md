# Vector Transport Free LBFGS on Riemannian Manifolds

The code for vector transport free LBFGS quasi-Newton's optimization on the Symmetric Positive Definite Riemannian manifolds.
This project reduces vector transport and Riemannian metric to identity and Euclidean inner product, respectively. 

This is the code for the following paper:
- Reza Godaz, Benyamin Ghojogh, Reshad Hosseini, Reza Monsefi, Fakhri Karray, Mark Crowley, "Vector Transport Free Riemannian LBFGS for Optimization on Symmetric Positive Definite Matrix Manifolds", arXiv preprint arXiv:2108.11019, 2021. 
- Link of paper: https://arxiv.org/abs/2108.11019

## Additional notes (related codes)

The code of manopt toolbox and its cautious RLBFGS solver: 
- https://github.com/NicolasBoumal/manopt
- https://github.com/NicolasBoumal/manopt/blob/master/manopt/solvers/bfgs/rlbfgs.m

The code of Mixest toolbox (it includes its RLBFGS solver): 
- http://visionlab.ut.ac.ir/resources.html
- http://visionlab.ut.ac.ir/resources/riemmix.zip

The code of another version of Mixest toolbox and its RLBFGS solver:
- http://visionlab.ut.ac.ir/mixest/
- https://github.com/utvisionlab/mixest
- https://github.com/utvisionlab/mixest/tree/master/mixest/auxiliary/manopt_solvers/lbfgs

The code of ROPTLIB toolbox and its RBFGS and RLBFGS solvers:
- https://github.com/whuang08/ROPTLIB
- https://github.com/whuang08/ROPTLIB/tree/master/Solvers
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/RBFGS.cpp
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/RBFGS.h
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/LRBFGS.h
- https://github.com/whuang08/ROPTLIB/blob/master/Solvers/LRBFGS.cpp
