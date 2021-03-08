class LineSearchBackTracking:
    """
    Back-tracking line-search based on linesearch.m in the manopt MATLAB
    package.
    """

    def __init__(self, contraction_factor=.5, optimism=2,
                 suff_decr=1e-4, maxiter=25, initial_stepsize=1):
        self.contraction_factor = contraction_factor
        self.optimism = optimism
        self.suff_decr = suff_decr
        self.maxiter = maxiter
        self.initial_stepsize = initial_stepsize

        self._oldf0 = None

    def search(self, objective, manifold, x, d, f0, df0):
        """
        Function to perform backtracking line-search.
        Arguments:
            - objective
                objective function to optimise
            - manifold
                manifold to optimise over
            - x
                starting point on the manifold
            - d
                tangent vector at x (descent direction)
            - f0
                objective function value at the current point
            - df0
                directional derivative at x along d
        Returns:
            - stepsize
                norm of the vector retracted to reach newx from x
            - newx
                next iterate suggested by the line-search
        """
        # Compute the norm of the search direction
        norm_d = manifold.norm(x, d)

        if self._oldf0 is not None:
            # Pick initial step size based on where we were last time.
            alpha = 2 * (f0 - self._oldf0) / df0
            # Look a little further
            alpha *= self.optimism
        else:
            alpha = self.initial_stepsize / norm_d
        alpha = float(alpha)

        # Make the chosen step and compute the cost there.
        newx = manifold.retr(x, alpha * d)
        newf = objective(newx)
        step_count = 1

        # Backtrack while the Armijo criterion is not satisfied
        while (newf > f0 + self.suff_decr * alpha * df0 and
               step_count <= self.maxiter):

            # Reduce the step size
            alpha = self.contraction_factor * alpha

            # and look closer down the line
            newx = manifold.retr(x, alpha * d)
            newf = objective(newx)

            step_count = step_count + 1

        # If we got here without obtaining a decrease, we reject the step.
        if newf > f0:
            alpha = 0
            newx = x

        stepsize = alpha * norm_d

        self._oldf0 = f0

        return stepsize, newx


class LineSearchAdaptive:
    '''
    Adaptive line-search
    '''

    def __init__(self, contraction_factor=.5, suff_decr=.5, maxiter=10,
                 initial_stepsize=1):
        self._contraction_factor = contraction_factor
        self._suff_decr = suff_decr
        self._maxiter = maxiter
        self._initial_stepsize = initial_stepsize
        self._oldalpha = None

    def search(self, objective, man, x, d, f0, df0):
        norm_d = man.norm(x, d)

        if self._oldalpha is not None:
            alpha = self._oldalpha
        else:
            alpha = self._initial_stepsize / norm_d
        alpha = float(alpha)

        newx = man.retr(x, alpha * d)
        newf = objective(newx)
        cost_evaluations = 1

        while (newf > f0 + self._suff_decr * alpha * df0 and
               cost_evaluations <= self._maxiter):
            # Reduce the step size.
            alpha *= self._contraction_factor

            # Look closer down the line.
            newx = man.retr(x, alpha * d)
            newf = objective(newx)

            cost_evaluations += 1

        if newf > f0:
            alpha = 0
            newx = x

        stepsize = alpha * norm_d

        # Store a suggestion for what the next initial step size trial should
        # be. On average we intend to do only one extra cost evaluation. Notice
        # how the suggestion is not about stepsize but about alpha. This is the
        # reason why this line search is not invariant under rescaling of the
        # search direction d.

        # If things go reasonably well, try to keep pace.
        if cost_evaluations == 2:
            self._oldalpha = alpha
        # If things went very well or we backtracked a lot (meaning the step
        # size is probably quite small), speed up.
        else:
            self._oldalpha = 2 * alpha

        return stepsize, newx


class LineSearchLBFGS_CautiousUpdate:
    """
    Back-tracking line-search based on solvers/bfgs/rlbfgs.m in the manopt MATLAB
    package.

    Based on paper: Huang, Wen, P-A. Absil, and Kyle A. Gallivan. 
    "A Riemannian BFGS method for nonconvex optimization problems." 
    In Numerical Mathematics and Advanced Applications ENUMATH 2015, 
    pp. 627-634. Springer, Cham, 2016.
    """

    def __init__(self, contraction_factor=.5, optimism=2,
                 suff_decr=1e-4, maxiter=25, initial_stepsize=1):
        self.contraction_factor = contraction_factor
        self.optimism = optimism
        self.suff_decr = suff_decr
        self.maxiter = maxiter
        self.initial_stepsize = initial_stepsize
        self._oldf0 = None

    def search(self, objective, manifold, x, d, f0, gradient):
        """
        Function to perform backtracking line-search.
        Arguments:
            - objective
                objective function to optimise
            - manifold
                manifold to optimise over
            - x
                starting point on the manifold
            - d
                - Hessian matrix * tangent vector at x (descent direction)
            - f0
                objective function value at the current point
            - gradient
                gradient function defined on the manifold
        Returns:
            - stepsize
                norm of the vector retracted to reach newx from x
            - newx
                next iterate suggested by the line-search
        """
        # Compute the norm of the search direction
        grad = gradient(x)
        search_dir = -grad
        norm_search_dir = manifold.norm(x, search_dir)

        if self._oldf0 is not None:
            # Pick initial step size based on where we were last time.
            gradnorm = manifold.norm(x, grad)
            df0 = -gradnorm**2
            alpha = 2 * (f0 - self._oldf0) / df0
            # Look a little further
            alpha *= self.optimism
        else:
            alpha = self.initial_stepsize / norm_search_dir
        alpha = float(alpha)

        # Make the chosen step and compute the cost there.
        newx = manifold.retr(x, alpha * d)
        newf = objective(newx)
        step_count = 1

        # Backtrack while the Armijo criterions are not satisfied
        t = -1  #--> iteration index
        while ((not Armijo_condition) and step_count <= self.maxiter):
            t += 1
            p = -gradient(x)
            if t == 0:
                zeta_t = Hessian_inverse @ p
            else:
                self.obtain_descent_direction(p, iteration)

            # Armijo condition:
            Armijo_condition = (newf <= f0 + self.suff_decr * alpha * df0)

            # Cautious update of Hessian approximation (denoted by B) with BFGS:
            sk = manifold.transp(x, newx, alpha * d)  #--> the step = newx - x = alpha * d
            # sk = sk / manifold.norm(x, sk)  # Computation of the BFGS step is invariant under scaling of sk and yk by a common factor. For numerical reasons, we scale sk and yk so that sk is a unit norm vector.
            yk = gradient(newx) - manifold.transp(x, newx, gradient(x))
            inner_sk_yk = manifold.inner(x, sk, yk)
            inner_sk_sk = manifold.norm(x, sk)**2    
            rho = 1 / inner_sk_yk
            cautious_condition = ((inner_sk_yk/inner_sk_sk) >= self._monotonic_function(gradnorm))
            if cautious_condition:
                B_tilde = manifold.transp(x, newx, B)
                # B = B_tilde - 
                # B = B + 
                pass

            # Reduce the step size
            alpha = self.contraction_factor * alpha

            # and look closer down the line
            newx = manifold.retr(x, alpha * d)
            newf = objective(newx)

            step_count = step_count + 1

        # If we got here without obtaining a decrease, we reject the step.
        if newf > f0:
            alpha = 0
            newx = x

        stepsize = alpha * norm_d

        self._oldf0 = f0

        return stepsize, newx

    def _monotonic_function(self, t):
        # In rlbfgs.m, it is written that f(t) = 1e-4*t is better than f(t) = t for bfgs
        return 1e-4 * t

    def obtain_descent_direction(p, iteration):
        pass

    
        

