import time
from copy import deepcopy
import autograd.numpy as np

from pymanopt.solvers.linesearch import LineSearchBackTracking
from pymanopt.solvers.solver import Solver


class LBFGS_WithWolfeConditions(Solver):
    """
    BFGS (quasi-Newton's method) algorithm based on
    solvers/bfgs/rlbfgs.m from the manopt MATLAB package.

    Based on paper: Huang, Wen, P-A. Absil, and Kyle A. Gallivan. 
    "A Riemannian BFGS method for nonconvex optimization problems." 
    In Numerical Mathematics and Advanced Applications ENUMATH 2015, 
    pp. 627-634. Springer, Cham, 2016.
    """

    def __init__(self, linesearch=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if linesearch is None:
            self._linesearch = LineSearchBackTracking()
        else:
            self._linesearch = linesearch
        self.linesearch = None
        self.max_store_depth = 10

    # Function to solve optimisation problem using BFGS.
    def solve(self, problem, x=None, reuselinesearch=False):
        """
        Perform optimization using gradient descent with linesearch.
        This method first computes the gradient (derivative) of obj
        w.r.t. arg, and then optimizes by moving in the direction of
        BFGS (which is the opposite direction to the gradient).
        Arguments:
            - problem
                Pymanopt problem setup using the Problem class, this must
                have a .manifold attribute specifying the manifold to optimize
                over, as well as a cost and enough information to compute
                the gradient of that cost.
            - x=None
                Optional parameter. Starting point on the manifold. If none
                then a starting point will be randomly generated.
            - reuselinesearch=False
                Whether to reuse the previous linesearch object. Allows to
                use information from a previous solve run.
        Returns:
            - x
                Local minimum of obj, or if algorithm terminated before
                convergence x will be the point at which it terminated.
        """
        man = problem.manifold
        verbosity = problem.verbosity
        objective = problem.cost
        gradient = problem.grad

        if not reuselinesearch or self.linesearch is None:
            self.linesearch = deepcopy(self._linesearch)
        linesearch = self.linesearch

        # If no starting point is specified, generate one at random.
        if x is None:
            x = man.rand()

        # Initialize iteration counter and timer
        iter = -1
        time0 = time.time()

        # Initialize inverse of Hessian:
        dimensionality_of_manifold = len(x)
        grad = gradient(x)
        Hessian_inverse = (1 / (man.inner(x, grad, grad) ** 0.5)) * np.eye(dimensionality_of_manifold)

        if verbosity >= 2:
            print(" iter\t\t   cost val\t    grad. norm")

        self._start_optlog(extraiterfields=['gradnorm'],
                           solverparams={'linesearcher': linesearch})

        # variables:
        self.x = [x]
        self.s = [None]
        self.y = [None]

        while True:
            # Calculate new cost, grad and gradnorm
            cost = objective(x)
            grad = gradient(x)
            gradnorm = man.norm(x, grad)
            iter = iter + 1

            if verbosity >= 2:
                print("%5d\t%+.16e\t%.8e" % (iter, cost, gradnorm))

            if self._logverbosity >= 2:
                self._append_optlog(iter, x, cost, gradnorm=gradnorm)

            # Descent direction calculation:
            desc_dir = self._obtain_descent_direction(p=-grad, iteration=iter%self.max_store_depth, Hessian_inverse=Hessian_inverse, man=man)

            # Perform line-search (with Armijo condition):
            stepsize, x = linesearch.search(objective, man, x, desc_dir, cost, -gradnorm**2)

            # Calculate some variables:
            x_tPlus1 = x
            x_t = self.x[-1]
            s_tPlus1 = man.transp(x_t, x_tPlus1, stepsize * desc_dir)  #--> the step = newx - x = stepsize * d
            # s_tPlus1 = s_tPlus1 / man.norm(x, s_tPlus1)  # Computation of the BFGS step is invariant under scaling of sk and yk by a common factor. For numerical reasons, we scale sk and yk so that sk is a unit norm vector.
            y_tPlus1 = gradient(x_tPlus1) - man.transp(x_t, x_tPlus1, gradient(x_t))

            # Update inverse of Hessian:
            inner_s_y = man.inner(x, s_tPlus1, y_tPlus1)
            inner_y_y = man.norm(x, y_tPlus1) ** 2
            Hessian_inverse = (inner_s_y / inner_y_y) * np.eye(dimensionality_of_manifold)

            # Store the variables:
            self.x = self._store_variable(list_=self.x, variable=x_tPlus1)
            self.s = self._store_variable(list_=self.s, variable=s_tPlus1)
            self.y = self._store_variable(list_=self.y, variable=y_tPlus1)

            stop_reason = self._check_stopping_criterion(
                time0, stepsize=stepsize, gradnorm=gradnorm, iter=iter)

            if stop_reason:
                if verbosity >= 1:
                    print(stop_reason)
                    print('')
                break

        if self._logverbosity <= 0:
            return x
        else:
            self._stop_optlog(x, objective(x), stop_reason, time0,
                              stepsize=stepsize, gradnorm=gradnorm,
                              iter=iter)
            return x, self._optlog

    def _obtain_descent_direction(self, p, iteration, Hessian_inverse, man, n_recursions=0):
        ### base of recursion:
        max_n_recursions = self.max_store_depth  
        if iteration == 0 or n_recursions >= max_n_recursions:
            desc_dir = Hessian_inverse @ p 
            return desc_dir
        ### body of recursion:
        p_tilde = p - (( man.inner(self.x[iteration], self.s[iteration], p) / man.inner(self.x[iteration], self.y[iteration], self.s[iteration]) ) * self.y[iteration])
        # p_tilde = vector_transport_adjoint on p_tilde  #---> TODO
        temp_ = self._obtain_descent_direction(p=p_tilde, iteration=iteration-1, Hessian_inverse=Hessian_inverse, man=man, n_recursions=n_recursions+1)
        p_hat = man.transp(self.x[iteration-1], self.x[iteration], temp_)
        rho = 1 / man.inner(self.x[iteration], self.y[iteration], self.s[iteration])
        return p_hat - ( rho * man.inner(self.x[iteration], self.y[iteration], p_hat) * self.s[iteration] ) + ( rho * man.inner(self.x[iteration], self.s[iteration], self.s[iteration]) * p)

    def _store_variable(self, list_, variable):
        list_.append(variable)
        if len(list_) > self.max_store_depth:
            list_[:self.max_store_depth] = list_[1:]
            list_.pop()
        return list_