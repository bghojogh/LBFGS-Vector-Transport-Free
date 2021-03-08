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


class LineSearch_WolfeConditions:
    """
    Line search with both Wolfe conditions (Armijo and curvature conditions)
    """

    def __init__(self, contraction_factor=.5, optimism=2,
                 suff_decr=1e-4, maxiter=25, initial_stepsize=1):
        self.contraction_factor = contraction_factor
        self.optimism = optimism
        self.suff_decr = suff_decr
        self.maxiter = maxiter
        self.initial_stepsize = initial_stepsize

        self._oldf0 = None

    def search(self, objective, manifold, x, d, f0, df0, gradient_function):
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
        Armijo_condition, curvature_condition = False, False
        while ((not Armijo_condition) and (not curvature_condition) and
               step_count <= self.maxiter):

            Armijo_condition = (newf <= f0 + self.suff_decr * alpha * df0)
            c2 = 0.9
            # weak Wolfe condition:
            # curvature_condition = (gradient_function(newx) @ manifold.transp(x, newx, d) >= c2 * df0)
            # strong Wolfe condition:
            curvature_condition = (abs(gradient_function(newx) @ manifold.transp(x, newx, d)) <= c2 * abs(df0))  

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

    
        

