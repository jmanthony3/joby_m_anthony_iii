#################################
## Preamble
# import necessary modules/tools

import math
import numpy as np
from scipy.integrate import quad
import sympy as sp
from sympy import simplify, solve
import sys
from types import FunctionType
#   #   #   #   #   #   #   #   #

#################################
## Universal Variables/Methods
# common error messages
must_be_expression = 'I am sorry. The input function must be an expression.'
must_be_collection = 'I am sorry. The input function must be a collection.'
opposite_signs = 'Initial guesses must yield opposite signs.'
solution_not_found = 'Solution could not be found with initial guess or tolerance.'
func_func = 'Input expression used.'
bad_matrix_A = "Warning: Matrix, A of 'Ac = g' is not strictly, diagonally dominant. Solution may be inaccurate."
# string outputs of polynomials
sym_x = sp.Symbol('x')

# common functions
def make_array(X, f):
    print(func_func)
    g = []
    for xi in X: g.append(f(xi))
    return g

def diagonality(A):
    """
    Warnings
    --------
    bad_matrix_A : string
        Printed string to say that matrix, A is not strictly, diagonally dominant. This could lead to poor solution of 'Ac = g'.'
    """
    i, j, n, diags, long = 0, 0, len(A) - 1 , [], []
    for k in range(-A.shape[0] + 1, A.shape[1]):
        b = []
        for a in A.diagonal(k): b.append(abs(a))
        j = n - i
        if i == j: long.append(b)
        else: diags.append(b)
        i += 1
    if not np.sum(long) >= np.sum(diags): print(bad_matrix_A)
#   #   #   #   #   #   #   #   #


#################################
## Classes
# categories of techniques

class test:                     # test class
    def test():                 # test function
        """Was the module loaded correctly?

        Returns
        -------
        success : string
            Prints a message of successful function call.
        """
        success = 'Test complete.'
        sys.exit(success)

class iterative_techniques:     # solving equation(s)
    """Finding solutions to equation(s).

    Attributes
    ----------
    single_variable : function
        Iterative techniques performed on functions of one variable.
    
    multi_variable : function
        Iterative techniques performed on functions/systems of equations of more than one variable.
    """
    class single_variable:      # implicitly find root of equation
        """Techniques to find solution to single-variable equations.

        Methods
        -------
        bisection():
            Bracketed root-finding technique.
        
        false_position():
            Bracketed fixed-point technique.
        
        fixed_point():
            Fixed point technique.
        
        max_iterations():
            Determines maximum number of iterations before exiting procedure.
        
        newton_raphson():
            Fixed point technique.
        
        secant_method():
            Fixed point technique.
        """
        def bisection(f, a, b, tol):
            """Given f(x) in [`a`,`b`] find x within tolerance, `tol`.
            Root-finding method: f(x) = 0.

            Parameters
            ----------
            f : expression
                Input function.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            tol : float
                Signed, specified power of tolerance until satisfying method.
            
            Returns
            -------
            P : list
                Aggregate collection of evaluated points, `p`.
            
            ERROR : list
                Propogation of `error` through method.
            
            I : list
                Running collection of iterations through method.

            Raises
            ------
            solution_not_found : string
                If initial guess or tolerance were badly defined.
            
            opposite_signs : string
                If initial guesses did not evaluate to have opposite signs.
            
            must_be_expression : string
                If input `f` was of array, list, tuple, etcetera...

            Notes
            -----
            Relying on the Intermediate Value Theorem, this is a bracketed, root-finding method. Generates a sequence {p_n}^{inf}_{n=1} to approximate a zero of f(x), `p` and converges by O(1 / (2**N)).

            Examples
            --------
            If  f(x) = x**3 + 4*x**2 = 10
            
            =>  f(x) = x**3 + 4*x**2 - 10 = 0
            """
            a, b, tol = float(a), float(b), float(tol)
            i, error = 0, tol*10        # initialize
            P, ERROR, I = [], [], []    # initialize lists
            # calculate if expression
            if isinstance(f,(FunctionType, sp.Expr)):
                # check if f(a) and f(b) are opposite signs
                if f(a)*f(b) < 0:
                    # exit by whichever condition is TRUE first
                    if error >= tol or \
                        i <= max_iterations(a, b, tol, 'bisection'):
                        x = (b - a)/2
                        p = a + x       # new value, p
                        P.append(p)
                        # adjust next bounds
                        if f(a)*f(p) > 0: a = p
                        else: b = p
                        error = abs(x)  # error of new value, p
                        ERROR.append(error); I.append(i)
                        i += 1          # iterate to i + 1
                    else: print(solution_not_found)
                # abort if f(a) is not opposite f(b)
                else: sys.exit(opposite_signs)
            # abort if not expression
            else: sys.exit(must_be_expression)
            return P, ERROR, I

        def false_position(f, k, a, b, p0, p1, tol):
            """Given f(x) and initial guesses, `p0` and `p1` in [`a`,`b`] find x within tolerance, `tol`.
            
            Root-finding problem: f(x) = 0. 
            
            !!! Use lowest k !!!

            Parameters
            ----------
            f : expression
                Input function.
            
            k : float
                Absolute maximum slope of `f`.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            p0 : float
                First initial guess.
            
            p1 : float
                Second initial guess.
            
            tol : float
                Specified tolerance until satisfying method.
            
            Returns
            -------
            P : list
                Aggregate collection of evaluated points, `p`.
            
            ERROR : list
                Propogation of `error` through method.
            
            I : list
                Running collection of iterations through method.

            Raises
            ------
            solution_not_found : string
                If initial guess or tolerance were badly defined.
            
            opposite_signs : string
                If initial guesses did not evaluate to have opposite signs.
            
            must_be_expression : string
                If input `f` was of array, list, tuple, etcetera...

            Notes
            -----
            Check that |g'(x)| <= (leading coefficient of g'(x)) for all x in [`a`,`b`].

            Theorem:
            1) Existence of a fixed-point:
                If g in C[`a`,`b`] and g(x) in C[`a`,`b`] for all x in [`a`,`b`], then function, g has a fixed point in [`a`,`b`].
            
            2) Uniqueness of a fixed point:
                If g'(x) exists on [`a`,`b`] and a positive constant, `k` < 1 exist with {|g'(x)| <= k  |  x in (`a`,`b`)}, then there is exactly one fixed-point, `p` in [`a`,`b`].

            Converges by O(linear) if g'(p) != 0, and O(quadratic) if g'(p) = 0 and g''(p) < M, where M = g''(xi) that is the error function.

            Examples 
            --------
            If  g(x) = x**2 - 2

            Then    p = g(p) = p**2 - 2
            
            =>  p**2 - p - 2 = 0
            """
            k, a, b, p0, p1, tol = float(k), float(a), float(b), float(p0), float(p1), float(tol)
            i, error = 0, tol*10        # initialize
            P, ERROR, I = [], [], []    # initialize lists
            # calculate if expression
            if isinstance(f,(FunctionType, sp.Expr)):
                # check if f(a) and f(b) are opposites signs
                if f(p0)*f(p1) < 0:
                    # exit by whichever condition is TRUE first
                    if error >= tol or \
                        i <= max_iterations(a, b, tol, 'false position', k, p0):
                        q0, q1 = f(p0), f(p1)
                        # new value, p
                        p = p1 - q1*(p1 - p0)/(q1 - q0)
                        P.append(p)
                        # error of new value, p
                        error = abs(p - p0)
                        ERROR.append(error); I.append(i)
                        # adjust next bounds
                        if f(p)*q1 < 0: p0 = p1
                        p1 = p
                        i += 1          # iterate to i + 1
                    else: print(solution_not_found)
                # abort if f(a) is not opposite f(b)
                else: sys.exit(opposite_signs)
            # abort if not expression
            else: sys.exit(must_be_expression)
            return P, ERROR, I

        def fixed_point(f, k, a, b, p0, tol):
            """Given f(x) and initial guess, `p0` in [`a`,`b`] find x within tolerance, `tol`.
            
            Root-finding problem: f(x) = 0. 
            
            !!! Use lowest k !!!

            Parameters
            ----------
            
            f : expression
                Input function.
            
            k : float
                Absolute maximum slope of `f`.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            p0 : float
                Initial guess.
            
            tol : float
                Specified tolerance until satisfying method.
            
            Returns
            -------
            
            P : list
                Aggregate collection of evaluated points, `p`.
            
            ERROR : list
                Propogation of `error` through method.
            
            I : list
                Running collection of iterations through method.

            Raises
            ------
            solution_not_found : string
                If initial guess or tolerance were badly defined.
            
            must_be_expression : string
                If input `f` was of array, list, tuple, etcetera...

            Notes
            -----
            Check that |g'(x)| <= (leading coefficient of g'(x)) for all x in [`a`,`b`].

            Theorem:
            1) Existence of a fixed-point:
                If g in C[`a`,`b`] and g(x) in C[`a`,`b`] for all x in [`a`,`b`], then function, g has a fixed point in [`a`,`b`].
            
            2) Uniqueness of a fixed point:
                If g'(x) exists on [`a`,`b`] and a positive constant, `k` < 1 exist with {|g'(x)| <= k  |  x in (`a`,`b`)}, then there is exactly one fixed-point, `p` in [`a`,`b`].

            Converges by O(linear) if g'(p) != 0, and O(quadratic) if g'(p) = 0 and g''(p) < M, where M = g''(xi) that is the error function.

            Examples 
            --------
            If  g(x) = x**2 - 2

            Then    p = g(p) = p**2 - 2
            
            =>  p**2 - p - 2 = 0
            """
            k, a, b, p0, tol = float(k), float(a), float(b), float(p0), float(tol)
            i, error = 0, tol*10        # initialize
            P, ERROR, I = [], [], []    # initialize lists
            # calculate if expression
            if isinstance(f,(FunctionType, sp.Expr)):
                # exit by whichever condition is TRUE first
                if error >= tol or \
                    i <= max_iterations(a, b, tol, 'fixed point', k, p0):
                    p = f(p0)           # new value, p
                    P.append(p)
                    error = abs(p - p0) # error of new value, p
                    ERROR.append(error); I.append(i)
                    p0 = p              # set future previous value
                    i += 1              # iterate to i + 1
                else: print(solution_not_found)
            # abort if not expression
            else: sys.exit(must_be_expression)
            return P, ERROR, I

        def max_iterations(a, b, tol, method, k=0, p0=0):
            """Find greatest integer for maximum iterations for tolerance.

            Parameters
            ----------
            
            f : expression
                Input function.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            tol : float
                Specified tolerance until satisfying method.
            
            method : string
                Selection of iterative method for iterations are needed.
            
            k : float
                Maximum possible slope of input function.
            
            p0 : float
                Initial guess for function solution.
            
            Returns
            -------
            N_max : integer
                Maximum number of iterations required for specified tolerance.

            Raises
            ------
            bad_method : string
                Prescribed method is not an available option.

            Notes
            -----
            Will round away from zero to higher integers.

            Examples
            --------
            If `method == 'bisection'` & a=1, b=2, and tol=-3, then:
            `N_max` >= -log(`tol`/(`b` - `a`))/log(2)

            `N_max` >= -log((10**(-3)/(2 - 1))/log(2)

            `N_max` >= 9.96

            `N_max` = 10

            Else, if a=1, b=2, tol=-3, p0=1.5, nd k=0.9, then:
            `N_max` >= log(`tol`/max('p0' - `a`, `b` - `p0`))/log(k)
            
            `N_max` >= log(10**(-3)/max(1.5 - 1, 2 - 1.5))/log(0.9)
            
            `N_max` >= log(10**(-3)/0.5)/log(0.9)
            
            `N_max` >= 58.98
            
            `N_max` >= 59
            """
            bad_method = 'I am sorry. The desired method must be: bisection, fixed point, newton raphson, secant method, or false position.'
            if method == 'bisection':
                N_max = math.ceil(-math.log(tol/(b - a))/math.log(2))
            elif method == ('fixed point', 'newton raphson', 'secant method', 'false position'):
                N_max = math.ceil(math.log(tol/max(p0 - a, b - p0))/math.log(k))
            else: sys.exit(bad_method)
            return N_max

        def newton_raphson(f, x, k, a, b, p0, tol):
            """Given f(x) and initial guess, `p0` in [`a`,`b`], find x within tolerance, `tol`.
            
            Root-finding problem: f(x) = 0. 
            
            !!! Use lowest k !!!

            Parameters
            ----------
            
            f : expression
                Input function.
            
            x : symbol
                Respected variable in derivative.
            
            k : float
                Absolute maximum slope of `f`.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            p0 : float
                Initial guess.
            
            tol : float
                Specified tolerance until satisfying method.
            
            Returns
            -------
            
            P : list
                Aggregate collection of evaluated points, `p`.
            
            ERROR : list
                Propogation of `error` through method.
            
            I : list
                Running collection of iterations through method.

            Raises
            ------
            solution_not_found : string
                If initial guess or tolerance were badly defined.
            
            must_be_expression : string
                If input `f` was of array, list, tuple, etcetera...

            Warnings
            --------
            f'(x) != 0.
            
            Not root-bracketed.

            Initial guess must be close to real solution; else, will converge to different root or oscillate (if symmetric).

            Notes
            -----
            Check that |g'(x)| <= (leading coefficient of g'(x)) for all x in [`a`,`b`].

            Technique based on first Taylor polynomial expansion of `f` about `p0` and evaluated at x = p. |p - p0| is assumed small; therefore, 2nd order Taylor term, the error, is small.

            Newton-Raphson has quickest convergence rate.

            This method can be viewed as fixed-point iteration.

            Theorem:
            1) Existence of a fixed-point:
                If g in C[`a`,`b`] and g(x) in C[`a`,`b`] for all x in [`a`,`b`], then function, g has a fixed point in [`a`,`b`].
            
            2) Uniqueness of a fixed point:
                If g'(x) exists on [`a`,`b`] and a positive constant, `k` < 1 exist with {|g'(x)| <= k  |  x in (`a`,`b`)}, then there is exactly one fixed-point, `p` in [`a`,`b`].

            Converges by O(linear) if g'(p) != 0, and O(quadratic) if g'(p) = 0 and g''(p) < M, where M = g''(xi) that is the error function.

            Examples 
            --------
            If  g(x) = x**2 - 2

            Then    p = g(p) = p**2 - 2
            
            =>  p**2 - p - 2 = 0
            """
            k, a, b, p0, tol = float(k), float(a), float(b), float(p0), float(tol)
            i, error = 0, tol*10        # initialize
            P, ERROR, I = [], [], []    # initialize lists
            # calculate if expression
            if isinstance(f,(FunctionType, sp.Expr)):
                # determine form of derivative
                df = sp.lambdify(x, sp.diff(f))
                # exit by whichever condition is TRUE first
                if error >= tol or \
                    i <= max_iterations(a, b, tol, 'newton raphson', k, p0):
                    fp0 = f(p0)
                    dfp0 = df(p0)
                    p = p0 - (fp0/dfp0) # new value, p
                    P.append(p)
                    error = abs(p - p0) # error of new value, p
                    ERROR.append(error); I.append(i)
                    p0 = p              # set future previous value
                    i += 1              # iterate to i + 1
                else: print(solution_not_found)
            # abort if not expression
            else: sys.exit(must_be_expression)
            return P, ERROR, I

        def secant_method(f, k, a, b, p0, p1, tol):
            """Given f(x) and initial guesses, `p0` and `p1` in [`a`,`b`], find x within tolerance, `tol`.
            Root-finding problem: f(x) = 0. 
            
            !!! Use lowest k !!!

            Parameters
            ----------
            
            f : expression
                Input function.
            
            k : float
                Absolute maximum slope of `f`.
            
            a : float
                Left-hand bound of interval.
            
            b : float
                Right-hand bound of interval.
            
            p0 : float
                First initial guess.
            
            p1 : float
                Second initial guess.
            
            tol : float
                Specified tolerance until satisfying method.
            
            Returns
            -------
            P : list
                Aggregate collection of evaluated points, `p`.
            
            ERROR : list
                Propogation of `error` through method.
            
            I : list
                Running collection of iterations through method.

            Raises
            ------
            solution_not_found : string
                If initial guess or tolerance were badly defined.
            
            opposite_signs : string
                If initial guesses did not evaluate to have opposite signs.
            
            must_be_expression : string
                If input `f` was of array, list, tuple, etcetera...

            Warnings
            --------
            Not root-bracketed.

            Notes
            -----
            Bypasses need to calculate derivative (as in Newton-Raphson).

            Check that |g'(x)| <= (leading coefficient of g'(x)) for all x in [`a`,`b`].

            Theorem:
            1) Existence of a fixed-point:
                If g in C[`a`,`b`] and g(x) in C[`a`,`b`] for all x in [`a`,`b`], then function, g has a fixed point in [`a`,`b`].
            
            2) Uniqueness of a fixed point:
                If g'(x) exists on [`a`,`b`] and a positive constant, `k` < 1 exist with {|g'(x)| <= k  |  x in (`a`,`b`)}, then there is exactly one fixed-point, `p` in [`a`,`b`].

            Converges by O(linear) if g'(p) != 0, and O(quadratic) if g'(p) = 0 and g''(p) < M, where M = g''(xi) that is the error function.

            Examples 
            --------
            If  g(x) = x**2 - 2

            Then    p = g(p) = p**2 - 2
            
            =>  p**2 - p - 2 = 0
            """
            k, a, b, p0, p1, tol = float(k), float(a), float(b), float(p0), float(p1), float(tol)
            i, error = 0, tol*10        # initialize
            P, ERROR, I = [], [], []    # initialize lists
            # calculate if expression
            if isinstance(f,(FunctionType, sp.Expr)):
                # check if f(a) and f(b) are opposite signs
                if f(p0)*f(p1) < 0:
                    # exit by whichever condition is TRUE first
                    if error >= tol or \
                        i <= max_iterations(a, b, tol, 'secant method', k, p0):
                        q0, q1 = f(p0), f(p1)
                        # new value, p
                        p = p1 - q1*(p1 - p0)/(q1 - q0)
                        P.append(p)
                        # error of new value
                        error = abs(p - p0)
                        ERROR.append(error); I.append(i)
                        p0, p1 = p1, p  # set future previous values
                        i += 1          # iterate to i + 1
                    else: print(solution_not_found)
                # abort if f(a) is not opposite f(b)
                else: sys.exit(opposite_signs)
            # abort if not expression
            else: sys.exit(must_be_expression)
            return P, ERROR, I

    class multi_variable:       # implicitly solve system of equations
        """Techniques to find solutions of multi-variate equations or systems of equations.

        Methods
        -------
        l_infinity_norm():
            Returns l_infinity norm between two vectors.
        
        l_2_norm():
            Returns l_2 norm between two vectors.
        
        jacobi():
            [x]_(k) = ( D^(-1)*(L + U) ) * [x]_(k - 1) + ( D^(-1) ) * [b]
        
        gauss_seidel():
            [x]_(k) = ( (D - L)^(-1) * U ) * [x]_(k - 1) + ( (D - L)^(-1) )*[b]
        
        successive_over_relaxation():
            [x]_(k) = ( (D - wL)^(-1) * ((1 - w)*D + w*U) ) * [x]_(k - 1) + w*( (D - w*L)^(-1) )*[b]
        """
        def l_infinity_norm(x, x0=0):
            """Maximum difference between absolute sum of i'th rows.

            Parameters
            ----------
            x : array or vector
                Newly approximated guess.
            
            x0 : array or vector
                Previously approximated guess.

            Returns
            -------
            `np.amax(norm_i)`
                Scalar value.

            Notes
            -----
            Best thought as "actual" distance between vectors.

            Also calculates infinity norm of matrix(ces).

            Examples
            --------
            [x0] = (1, 1, 1)^(t)

            [x] = (1.2001, 0.99991, 0.92538)^(t)

            ||x0 - x|| = max{|1 - 1.2001|, |1 - 0.99991|, |1 - 0.92538|}

            ||x0 - x|| = 0.2001
            """
            # initialize loop
            i, norm_i = 0, np.zeros_like(x)
            while i < len(x):
                j = 0
                while j < len(x[0]):
                    # evaluate and store norm, ||.||
                    if x0 != 0: norm_i[i] += abs(x[i][j] - x0[i][j])
                    else: norm_i[i] += abs(x[i][j])
                    j += 1      # iterate to j + 1 column
                i += 1          # iterate to i + 1 row
            # return the l_infinity norm
            return np.amax(norm_i)

        def l_2_norm(x, x0=0):
            """Square root of sum of differences squared along i'th row.

            Parameters
            ----------
            x : array
                Newly approximated guess.
            
            x0 : array
                Previously approximated guess.

            Returns
            -------
            `math.sqrt(np.sum(norm_i))`
                Scalar value.

            Examples
            --------
            [x0] = (1, 1, 1)^(t)

            [x] = (1.2001, 0.99991, 0.92538)^(t)

            ||x0 - x|| = sqrt[ (1 - 1.2001)^2 \
                + (1 - 0.99991)^2 + (1 - 0.92538)^2 ]

            ||x0 - x|| = 0.21356
            """
            # initialize loop
            i, norm_i = 0, np.zeros_like(x)
            while i < len(x):
                j = 0
                while j <= len(x[0]):
                    # evaluate and store norm, ||.||
                    if x0 != 0: norm_i[i] += (x[i] - x0[i])**2
                    else: norm_i[i] += x[i]**2
                    j += 1      # iterate to j + 1 column
                i += 1          # iterate to i + 1
            # return the l_2 norm
            return math.sqrt(np.sum(norm_i))

        def eigen_values(A):
            r_sym = sp.Symbol('r')
            identity_matrix = np.identity(A)*r_sym
            lambda_matrix = 

        def vector_converge(A, x0, b, N, tol, type, method, w=0):
            """Given [`A`]*[`x`] = [`b`], use `method` and `type` to find [x].

            Parameters
            ----------
            A : matrix
                Characteristic matrix.
            
            x0 : vector
                Dimensions of system of equations.
            
            b : vector
                Input vector.
            
            N : int
                Maximum number of iterations.
            
            tol : float
                Power of desired constraint for final solution.
            
            type : string
                Prescription of desired norm.
            
            method : string
                Actual technique.
            
            w : float
                Relaxation parameter.
            
            Returns
            -------
            X_matrix : array
                Finally evaluated solution.
            
            NORM : list
                Aggregate of yielded norms.
            
            K : list
                Running collection of iterations through method.

            Raises
            ------
            bad_matrix : string
                If [`A`] is not square.
            
            bad_x0 : string
                If {`x0`} is neither n x 1 nor 1 x n array.
            
            bad_b : string
                If {`b`} is neither n x 1 nor 1 x n array.
            
            bad_N : string
                If iterations constraints is not an integer.

            solution_not_found : string
                If initial guess or tolerance were badly defined.

            Notes
            -----
            Unless stated, `w = 0`.
            """
            def jacobi(x0):         
                while i < n:
                    j, y = 0, 0.
                    while j < n:
                        if j != i:
                            y += A[i][j]*x0[j]
                            j += 1
                    xi[i] = (-y + b[i])/A[i][i]
                    i += 1
                return xi
            def gauss_seidel(x0):
                while i < n:
                    j, y1, y2 = 0, 0., 0.
                    while j < i-1:
                        y1 += A[i][j]*xi[j]
                        j += 1
                    j = i + 1
                    while j < n:
                        y2 += A[i][j]*x0[j]
                        j += 1
                    xi[i] = (-y1 - y2 + b[i])/A[i][i]
                    i += 1
                return xi
            def successive_over_relaxation(x0):
                while i < n:
                    gauss_seidel(x0)
                    xi[i] = (1 - w)*x0[i] + w*gauss_seidel(x0)
                    i += 1
                return xi
            bad_matrix = 'Characteristic matrix, A must be square!'
            bad_x0 = 'Systems vector, x0 must be n x 1 or 1 x n array!'
            bad_b = 'Input vector, b must be n x 1 or 1 x n array!'
            bad_N = "Maximum iterations, N must be an integer greater than zero."
            if len(A) != len(A[0]): sys.exit(bad_matrix)
            if len(x0[0]) > 1 and len(x0[1]) > 1: sys.exit(bad_x0)
            if len(b[0]) > 1 and len(b[1]) > 1: sys.exit(bad_b)
            if N <= 0 or not isinstance(N, int): sys.exit(bad_N)
            tol = float(tol)
            n = len(x0)
            k, x0, b, norm = 0, np.reshape(x0,(n,1)), np.reshape(b,(n,1)), tol*10
            xi = np.zeros_like(x0)
            X, NORM, K = [], [], [] 
            X.append(x0); K.append(k)
            if norm > tol or k <= N:
                i = 0
                if method == 'jacobi': xi = jacobi(x0)
                if method == 'gauss_seidel': xi = gauss_seidel(x0)
                if method == 'successive_over_relaxation': xi = successive_over_relaxation(x0)
                if type == 'l_infinity': norm = l_infinity_norm(xi, x0)
                if type == 'l_2': norm = l_2_norm(xi, x0)
                X.append(xi); NORM.append(norm); K.append(k)
                x0 = xi
                k += 1
            else: sys.exit(solution_not_found)
            m, n = len(X[0]), len(X)
            X_matrix, j = np.zeros((m,n)), 0
            while j < n:
                i = 0
                while i < m:
                    X_matrix[i][j] = float(X[j][i])
                    i += 1
                j += 1
            return X_matrix, NORM, K

class interpolations:           # use data set to build polynomial
    """Finding solutions to equation(s).

    Attributes
    ----------
    cubic_spline : function
        Constructs a cubic spline polynomial through some data.
    
    hermite : function
        Oscullating curve that matches data set and its trend.

    lagrange : function
        Constructs a Lagrangian polynomial through some data.
    
    linear_least_squares : function
        Yields a polynomial of desired degree with minimal error to data.
    
    newton_difference : function
        Polynomial built direclty from data set.
    """
    def cubic_spline(X, f, condition, fp=0):
        """Given a domain and range, construct a spline polynomial within interval by some condition.

        Parameters
        ----------
        X : array
            Input domain.
        
        f : array or expression
            Desired/Found range of interest.
        
        condition : string
            Method by which to construct spline polynomial.
        
        fp : array or expression
            Derivative at each point in `f`.
        
        Returns
        -------
        Y : array
            Finally evaluated solution.
        
        splines_j : list
            Aggregate of splines on each interval.
        
        spline : string
            Totally constructed spline polynomial.

        Raises
        ------
        bad_X : string
            If {`X`} is neither n x 1 nor 1 x n array.
        
        bad_f : string
            If `f` is not an expression or function and is not an n x 1 or 1 x n array.
        
        bad_data : string
            If {`X`} and {`f`} are of unequal length.
        
        bad_fp : string
            If `fp` is not an expression or function and is not an n x 1 or 1 x n array.
        
        missing_fp : string
            Output message that derivative data or expression is missing.

        bad_condition : string
            If indicated condition was neither `'clamped'` nor `'natural'`.

        Warnings
        --------
        `fp` will be calculated if not specified.

        See Also
        --------
        num_diff_and_int.endpoint() : Relies on another technique to find derivatives at endpoints if not explicitly provided by data, `fp` nor an expression.

        diagonality() : Determines whether input matrix is strictly, diagonally dominant.

        Notes
        -----
        Method uses many, low-ordered polynomials to fit larger data sets. This minimizes computational load, which conversely greatly increases for larger data sets that yield high-ordered polynomials.

        General form: 
        S(x) = sum_(j=0)^(n){aj + bj(x - xj) + cj(x - xj)^2 + dj(x - xj)^3}

        Clamped splines fit the constructed polynomial to the given data and its der
        ivatives at either endpoint.

        If selected `condition` is `'natural'`, then `fp = 0`, because derivative is assumed to be straight line outside of data set.

        Definitions of cubic spline conditions:
        a) S(x) is a cubic polynomial, Sj(x) on sub-interval [x_(j), x_(j + 1)] for each j = 0, 1, ..., n - 1;

        b) Sj(x_(j)) = f(x_(j)) and Sj(x_(j + 1)) = f(x_(j + 1)) for each j = 0, 1, ..., n - 1;

        c) S_(j + 1)(x_(j + 1)) = Sj(x_(j + 1)) for each j = 0, 1, ..., n - 2;

        d) S_(j + 1)'(x_(j + 1)) = Sj'(x_(j + 1)) for each j = 0, 1, ..., n - 2;

        e) One of the following conditions is satisfied:
            1) S''(x0) = S''(xn) = 0                ->  `'natural'`
            
            2) S'(x0) = f'(x0) and S'(xn) = f'(xn)  ->  `'clamped'`
        """
        def clamped():
            Y, YP = f, fp
            # STEP 1:   build list, h_i
            i, H = 0, np.zeros(n)
            while i < n:
                H[i] = X[i+1] - X[i]
                i += 1
            # STEP 2:   define alpha list endpoints
            A, AP, ALPHA = Y, YP, np.zeros(m)
            ALPHA[0] = 3*(A[1] - A[0])/H[0] - 3*AP[0]
            ALPHA[n] = 3*AP[n] - 3*(A[n] - A[n-1])/H[n-1]
            # STEP 3:   build list, alpha_i
            i = 1
            while i <= n-1:
                ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
                i += 1
            # Algorithm 6.7 to solve tridiagonal
            # STEP 4:   define l, mu, and z first points
            L, MU, Z, C = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)
            L[0], MU[0] = 2*H[0], 0.5
            Z[0] = ALPHA[0]/L[0]
            # STEP 5:   build lists l, mu, and z
            i = 1
            while i <= n-1:
                L[i] = 2*(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
                MU[i] = H[i]/L[i]
                Z[i] = (ALPHA[i] - H[i-1]*Z[i-1])/L[i]
                i += 1
            # STEP 6:   define l, z, and c endpoints
            L[n] = H[n-1]*(2-MU[i-1])
            Z[n] = (ALPHA[n] - H[n-1]*Z[n-1])/L[n]
            C[n] = Z[n]
            # STEP 7:   build lists c, b, and d
            B, D, i, j = np.zeros(n), np.zeros(n), 1, 0
            while i <= n:
                j = n-i
                C[j] = Z[j] - MU[j]*C[j+1]
                B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2*C[j])/3
                D[j] = (C[j+1] - C[j])/(3*H[j])
                i += 1
            return Y, A, B, C, D
        def natural():
            Y = np.zeros(m)
            if isinstance(f,(FunctionType, sp.Expr)):
                Y = make_array(X, f)
            else: Y = f
            # STEP 1:   build list, h_i
            H, i = np.zeros(n), 0
            while i < n:
                H[i] = X[i+1] - X[i]
                i += 1
            # STEP 2:   build list, alpha_i
            A, ALPHA = Y, np.zeros(m)
            i = 1
            while i <= n-1:
                ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
                i += 1
            # Algorithm 6.7 to solve tridiagonal
            # STEP 3:   define l, mu, and z first points
            L, MU, Z, C = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)
            L[0], MU[0], Z[0] = 1, 0, 0
            # STEP 4:   build lists l, mu, and z
            i = 1
            while i <= n-1:
                L[i] = 2*(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
                MU[i] = H[i]/L[i]
                Z[i] = (ALPHA[i] - H[i-1]*Z[i-1])/L[i]
                i += 1
            # STEP 5:   define l, z, and c endpoints
            L[n], Z[n], C[n] = 1, 0, 0
            # STEP 6:   build lists c, b, and d
            B, D, i, j = np.zeros(n), np.zeros(n), 1, 0
            while i <= n:
                j = n-i
                C[j] = Z[j] - MU[j]*C[j+1]
                B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2*C[j])/3
                D[j] = (C[j+1] - C[j])/(3*H[j])
                i += 1
            return Y, A, B, C, D
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_f = 'Input range was neither function nor expression and not an n x 1 or 1 x n array.'
        bad_data = 'Arrays X and f must be of equal length.'
        bad_fp = 'Derivative range was neither function nor expression and not an n x 1 or 1 x n array.'
        bad_fp_data = 'Arrays X, f, and fp must be of equal length.'
        missing_fp = 'Missing derivative data or expression.'
        bad_condition = "Desired condition was not understood. Expected 'clamped' or 'natural'."
        if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
        if not isinstance(f, (FunctionType, sp.Expr)):
            if len(f[0]) > 1 and len(f[1]) > 1: sys.exit(bad_f)
            if len(X) != len(f): sys.exit(bad_data)
        if fp != 0:
            if not isinstance(fp, (FunctionType, sp.Expr)):
                if len(fp[0]) > 1 and len(fp[1]) > 1: sys.exit(bad_fp)
                if len(X) != len(fp): sys.exit(bad_fp_data)
            if isinstance(fp, (FunctionType, sp.Expr)): fp = make_array(X, fp)
        if fp == 0:
            if isinstance(f,(FunctionType, sp.Expr)):
                g = sp.diff(f, sym_x)
                fp = make_array(X, g)
            if not isinstance(f,(FunctionType, sp.Expr)):
                fp = []
                if len(f) > 3:
                    fp.append(num_diff_and_int.endpoint(X, f, X[1]-X[0], 'three', 'left'))
                    fp.append(num_diff_and_int.endpoint(X, f, X[-2]-X[-1], 'three', 'right'))
                if len(f) > 5:
                    fp.append(num_diff_and_int.endpoint(X, f, X[1]-X[0], 'five', 'left'))
                    fp.append(num_diff_and_int.endpoint(X, f, X[-2]-X[-1], 'five', 'right'))
            else: sys.exit(missing_fp)
        m = len(X)
        n = m - 1
        if condition == 'clamped': Y, A, B, C, D = clamped()
        if condition == 'natural': Y, A, B, C, D = natural()
        else: sys.exit(bad_condition)
        j, splines_j = 0, []
        while j <= n-1:
            xj, aj, bj, cj, dj = X[j], A[j], B[j], C[j], D[j]
            sj = aj + bj*(sym_x - xj) + cj*(sym_x - xj)**2 + dj*(sym_x - xj)**3
            splines_j.append(sj)
            j += 1
        spline = sum(splines_j)
        return Y, splines_j, spline

    def hermite(X, FX, FP=0):
        """Given a domain and range, construct a Hermetic polynomial.

        Parameters
        ----------
        X : array
            Input domain.
        
        FX : array
            Desired/Found range of interest.
        
        FP : array or expression
            Derivative at each point in `FX`.
        
        Returns
        -------
        polynomial : string
            Totally constructed Hermetic polynomial.

        Raises
        ------
        bad_X : string
            If {`X`} is neither n x 1 nor 1 x n array.
        
        bad_FX : string
            If {`FX`} is neither n x 1 nor 1 x n array.
        
        bad_data : string
            If {`X`} and {`FX`} are of unequal length.
        
        bad_FP : string
            If `FP` is not an expression or function and is not an n x 1 or 1 x n array.

        bad_FP_data : string
            If {`X`}, {`FX`}, or {`FP`} are of unequal lengths.
        
        missing_FP : string
            If `FP = 0` and `FX` is not an expression, then missing derivative data or expression.

        Warnings
        --------
        `FP` calculated if not specified.

        Slow computation time for larger data sets.

        See Also
        --------
        make_array() : Prints string that expression was used to make array.
        
        Notes
        -----
        Oscullating curve incorporates Taylor and Lagrangian polynomials to kiss the data and match each data point's derivatives. Which fits the curve to the shape of the data and its trend.
        """
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_FX = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X and FX must be of equal length.'
        bad_FP = 'Derivative range was neither function nor expression and not an n x 1 or 1 x n array.'
        bad_FP_data = 'Arrays X, FX, and FP must be of equal length.'
        missing_FP = 'Missing derivative data or expression.'
        if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
        if not isinstance(FX, (FunctionType, sp.Expr)):
            if len(FX[0]) > 1 and len(FX[1]) > 1: sys.exit(bad_FX)
            if len(X) != len(FX): sys.exit(bad_data)
        if isinstance(FX,(FunctionType, sp.Expr)): FX = make_array(X, FX)
        if FP != 0:
            if not isinstance(FP, (FunctionType, sp.Expr)):
                if len(FP[0]) > 1 and len(FP[1]) > 1: sys.exit(bad_FP)
                if len(X) != len(FP): sys.exit(bad_FP_data)
            if isinstance(FP,(FunctionType, sp.Expr)): FP = make_array(X, FP)
        if FP == 0:
            if isinstance(FX,(FunctionType, sp.Expr)):
                g = sp.diff(FX, sym_x)
                FP = make_array(X, g)
            else: print(missing_FP)
        n = len(X)
        i, Q, Z = 0, np.zeros((2*n+1,2*n+1)), np.zeros((n+1,1))
        while i <= n:
            Z[2*i], Z[2*i + 1] = X[i], X[i]
            Q[2*i][0], Q[2*i + 1][0] = FX[i], FX[i]
            Q[2*i + 1][1] = FP[i]
            if i != 0: Q[2*i][1] = (Q[2*i][0] - Q[2*i - 1][0]) \
                / (Z[2*i] - Z[2*i - 1])
            i += 1
        i = 2
        while i <= 2*n + 1:
            j = 2
            while j <= i:
                Q[i][j] = (Q[i][j - 1] - Q[i - 1][j - 1]) \
                / (Z[i] - Z[i - j])
                j += 1
            i += 1
        i, x, terms = 0, 1, []
        while i <= n:
            j, xi = 2*i, (sym_x - X[i])
            qjj, qj1 = Q[j][j], Q[j + 1][j + 1]
            terms.append(qjj*x)
            x = x*xi
            terms.append(qj1*x)
            x = x*xi
            i += 1
        polynomial = simplify(sum(terms))
        return polynomial
    
    def lagrange(X, Y):
        """Given a domain and range, construct a Lagrangian polynomial.

        Parameters
        ----------
        X : array
            Input domain.
        
        Y : array or expression
            Desired/Found range of interest.
        
        Returns
        -------
        yn : list
            Aggregate of Lagrangian terms.
        
        sum(yn) : string
            Totally constructed Lagrangian polynomial.
        
        bound : list
            Propogation of error through construction.
        
        sum(bound)
            Total error.

        Raises
        ------
        bad_X : string
            If {`X`} is neither n x 1 nor 1 x n array.
        
        bad_Y : string
            If {`Y`} is neither n x 1 nor 1 x n array.
        
        bad_data : string
            If {`X`} and {`Y`} are of unequal length.

        Warnings
        --------
        Polynomial will quickly begin to oscillate for larger data sets.

        See Also
        --------
        make_array() : Prints string that expression was used to make array.

        Notes
        --------
        Finds a polynomial of degree n-1.

        Polynomial is of the following form:
        P(x) = f(x0)L_(n,0)(x) + ... + f(xn)L_(n,n)(x), where

        L_(n,k) = prod_(i=0, i!=k)^(n) (x - xi)/(xk - xi)

        Examples
        --------
        A Lagrange polynomial between (2,4) and (5,1) would be found as follows:
        L_(0)(x) = (x - 5)/(2 - 5) = -(x - 5)/3

        L_(1)(x) = (x - 2)/(5 - 2) = (x - 2)/3

        =>  P(x)    = (4)*(-(x - 5)/3) + (1)*((x - 2)/3)
                    = -x + 6
        """
        def term(xk, yk):
            num, den, L_k = [], [], []
            for xl in X:
                if xl != xk:
                    num.append(sym_x-xl)
                    den.append(xk-xl)
            L_k = (np.divide(np.prod(num), np.prod(den)))
            return L_k * yk
        def error(n, xi):
            i, roots, g, xi_error = 0, [], [], []
            while i <= n:
                root = X[i]
                roots.append(sym_x - root)
                g = np.prod(roots)
                k = 0
                while k <= n:
                    xi = sp.diff(xi, sym_x)
                    k += 1
                dxi = np.abs(xi.evalf(subs={sym_x: root})/(math.factorial(k)))
                xi_error.append(np.abs(dxi))
                xi_err = np.max(xi_error)
                g_prime = sp.diff(g, sym_x)
                r = solve(g_prime)
                if i == 0:
                    x = g_prime
                    gx = g.evalf(subs={sym_x: x})
                if i == 1:
                    x = r[0]
                    gx = g.evalf(subs={sym_x: x})
                else:
                    R = []
                    for s in r:
                        if not isinstance(s, complex):
                            R.append(g.evalf(subs={sym_x: s}))
                    gx = np.amax(np.abs(R))
                i += 1
            return np.abs(xi_err*gx)
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_Y = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X and Y must be of equal length.'
        if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
        if not isinstance(Y,(FunctionType, sp.Expr)):
            if len(Y[0]) > 1 and len(Y[1]) > 1: sys.exit(bad_Y)
            if len(X) != len(Y): sys.exit(bad_data)
        if isinstance(Y,(FunctionType, sp.Expr)): Y = make_array(X, Y)
        yn, bound = [], []
        for xk in X:
            k = X.index(xk)
            yn.append(term(xk, Y[k]))
            bound.append(error(k, sum(yn)))
        return yn, sum(yn), bound, sum(bound)

    def linear_least_squares(X_i, Y_i, n):
        """Given a domain and range, construct some polynomial.

        Parameters
        ----------
        X_i : array
            Input domain.
        
        Y_i : array or expression
            Desired/Found range of interest.
        
        n : int
            Degree of polynomial.
        
        Returns
        -------
        polynomial : string
            Totally constructed Lagrangian polynomial.
        
        E : float
            Total error.

        Raises
        ------
        bad_X : string
            If {`X_i`} is neither n x 1 nor 1 x n array.
        
        bad_Y : string
            If {`Y_i`} is neither n x 1 nor 1 x n array.
        
        bad_data : string
            If {`X_i`} and {`Y_i`} are of unequal length.
        
        bad_n : string
            If prescribed `n` is not an integer or is zero.
        """
        def poly(X):
            terms, k = [], 0
            for x in X:
                terms.append(x*(sym_x**k))
                k += 1
            p = simplify(sum(terms))
            err, i = 0, 0
            for x_i in X_i:
                px = p.subs(sym_x, x_i)
                err += (Y_i[i] - px)**2
                i += 1
            return p, err
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_Y = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X_i and Y_i must be of equal length.'
        bad_n = 'Degree of polynomial must be integer and non-zero.'
        if len(X_i[0]) > 1 and len(X_i[1]) > 1: sys.exit(bad_X)
        if len(Y_i[0]) > 1 and len(Y_i[1]) > 1: sys.exit(bad_Y)
        if len(X_i) != len(Y_i): sys.exit(bad_data)
        if not isinstance(n,(int)) or n == 0: sys.exit(bad_n)
        m = len(X_i)
        A, x = np.zeros((n+1, n+1)), np.zeros((n+1,1))
        b, i = np.zeros_like(x), 0
        while i <= n:
            j = 0
            while j <= n:
                a_ij, k = 0, 0
                while k < m:
                    a_ij += (X_i[k])**(i + j)
                    k += 1
                A[i][j] = a_ij
                j += 1
            b_i, k = 0, 0
            while k < m:
                b_i += Y_i[k]*(X_i[k]**(i))
                k += 1
            b[i] = b_i
            i += 1
        x = np.transpose(np.linalg.solve(A, b))
        X, terms, k = x[0], [], 0
        for x in X:
            terms.append(x*(sym_x**k))
            k += 1
        polynomial = simplify(sum(terms))
        E, i = 0, 0
        for x_i in X_i:
            E += (Y_i[i] - polynomial.subs(sym_x, x_i))**2
            i += 1
        return polynomial, E

    def newton_difference(X, FX, x0, direction):
        """Given a domain and range, construct some polynomial by Newton's Divided Difference.

        Parameters
        ----------
        X : array
            Input domain.
        
        FX : array or expression
            Desired/Found range of interest.

        x0 : float
            Point about which polynomial is evaluated.
        
        direction : string
            `'forward'` or `'backward'` construction.
        
        Returns
        -------
        polynomial : string
            Totally constructed polynomial.
        
        px : float
            Evaluation of `polynomial` at `x`.

        Raises
        ------
        bad_X : string
            If {`X_i`} is neither n x 1 nor 1 x n array.
        
        bad_FX : string
            If {`FX`} is neither n x 1 nor 1 x n array.
        
        bad_data : string
            If {`X`} and {`FX`} are of unequal length.
        
        bad_direction : string
            If `direction` is neither `'forward'` nor `'backward'`.

        See Also
        --------
        make_array() : Prints string that expression was used to make array.

        Notes
        -----
        Polynomials best made with even spacing in `X`; although, this is not completely necessary.
        """
        def fterm(i, j):
            fij = (fxn[i][j] - fxn[i-1][j])/(fxn[i][0] - fxn[i-j][0])
            return fij
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_FX = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X and FX must be of equal length.'
        bad_direction = "Supplied direction was not understood. Please specify 'forward' or 'backward'."
        if not isinstance(FX,(FunctionType, sp.Expr)):
            if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
            if len(FX[0]) > 1 and len(FX[1]) > 1: sys.exit(bad_FX)
            if len(X) != len(FX): sys.exit(bad_data)
        if isinstance(FX,(FunctionType, sp.Expr)): FX = make_array(X, FX)
        if direction != 'forward' or 'backward': sys.exit(bad_direction)
        x0 = float(x0)
        m = len(X)
        n = m + 1
        fxn, coeff, term, poly = np.zeros((m,n)), [], [], []
        m, n = m - 1, n - 1     # change m and n from length to index
        j, fxn[:m,0], fxn[:m,1] = 1, X, FX
        while j <= n:
            i = 1
            while i <= m:
                fk = fterm(i, j)
                fxn[i][j+1] = fk
                if direction == 'forward' and i == j:
                    coeff.append(fk)
                if direction == 'backward' and i == m - 1:
                    coeff.append(fk)
                i += 1
            j += 1
        for c in coeff:
            k = coeff.index(c)
            term.append(sym_x - X[k])
            poly.append(c*np.prod(term))
        if direction == 'forward': polynomial = simplify(sum(poly) + FX[0])
        if direction == 'backward': polynomial = simplify(sum(poly) + FX[m])
        px = polynomial.subs(sym_x, x0)
        return polynomial, px

class num_diff_and_int:         # computational differentiation/integration
    """Differentiate and integrate some function.

    Attributes
    ----------
    composite_simpson : function
        Iterative techniques performed on functions of one variable.
    
    composite_trapz : function
        Iterative techniques performed on functions/systems of equations of more than one variable.
    
    endpoint : function
        Yields the derivative of a data set's endpoints.
    
    gaussian_quadrature : function
        Optimization of Simpson's Rule with Legendre polynomial.
    
    integrate : function
        Short-hand function that actually uses Gaussian Quadrature.

    midpoint : function
        Yields the derivative at some point within the data set.

    richard_extrapolation : 
    """

    def composite_simpson(f, h, X=0, a=0, b=0):
        """Find the integral of a function within some interval, using Simpson's Rule.

        Parameters
        ----------
        f : expression
            Polynomial equation that defines graphical curve.
        
        h : float
            Step-size through interval.
        
        X : list
            Domain over which `f` is evaluated.
        
        a : float
            Left-hand bound of interval.
        
        b : float
            Right-hand bound of interval.
        
        Returns
        -------
        XJ : list
            Values of domain at which `f` was analyzed.
        
        YJ : list
            Evaluations of `f` from domain.
        
        F : float
            Total area under curve, `f`.

        Raises
        ------
        bad_X : string
            If {`X_i`} is neither n x 1 nor 1 x n array.
        
        bad_f : string
            If {`f`} is not an expression.

        Warns
        -----
        func_func : string
            Evaluate input expression for Newton difference approximation.
        
        Notes
        -----
        `X = 0` if not a list nor n x 1 or 1 x n array.

        Unless specified and if `X` is defined, `a` and `b` will be the minimum and maximum, respectively, of `X`.

        Theorem:
        Let f be in C4[a,b], n be even, h = (b-a)/n, and xj = a + jh for j = 0, 1, ..., n. There exists a mu in (a,b) for which the quadrature for n sub-intervals can be written with its error term as:
        int_(a)^(b)f(x)dx = h[f(a) + 2*[sum_(j=1)^(n/2 - 1){f(x_(2j))}] + 4*[sum_(j=1)^(n/2){f(x_(2j-1))}] + f(b)]/3 - (b-a)*(h^4)f''''(mu)/180.

        Where: (b-a)*(h^4)f''''(mu)/180 -> O(h^4)
        """
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_f = 'Input range must be expression, not list or tuple.'
        func_func = 'Input expression used.'
        if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
        if not isinstance(f,(FunctionType, sp.Expr)):
            if len(f[0]) > 1 and len(f[1]) > 1: sys.exit(bad_X)
            else: sys.exit(bad_f)
        if isinstance(f,(FunctionType, sp.Expr)): print(func_func)
        if X != 0: a, b = min(X), max(X)
        h, a, b = float(h), float(a), float(b)
        n = math.ceil((b-a)/h)
        XJ1, XJ2, XJ, = [], [], []
        YJ1, YJ2, YJ, = [], [], []
        XJ.append(a); YJ.append(f(a))
        j, z1 = 1, 0
        while j <= (n/2)-1:
            xj = a + 2*j*h
            yj = f(xj)
            XJ1.append(xj); YJ1.append(yj)
            z1 += yj
            j += 1
        k, z2 = 1, 0
        while k <= n/2:
            xj = a + (2*k - 1)*h
            yj = f(xj)
            XJ2.append(xj); YJ2.append(yj)
            z2 += yj
            k += 1
        l = 0
        while l < len(XJ1):
            XJ.append(XJ2[l]); YJ.append(YJ2[l])
            XJ.append(XJ1[l]); YJ.append(YJ1[l])
            l += 1
        XJ.append(XJ2[l]); YJ.append(YJ2[l])
        XJ.append(b); YJ.append(f(b))
        F = h/3*(f(a) + 2*z1 + 4*z2 + f(b))
        return XJ, YJ, F

    def composite_trapezoidal(f, h, X=0, a=0, b=0):
        """Find the integral of a function within some interval, using Trapezoidal Rule.

        Parameters
        ----------
        f : expression
            Polynomial equation that defines graphical curve.
        
        h : float
            Step-size through interval.
        
        X : list
            Domain over which `f` is evaluated.
        
        a : float
            Left-hand bound of interval.
        
        b : float
            Right-hand bound of interval.
        
        Returns
        -------
        XJ : list
            Values of domain at which `f` was analyzed.
        
        YJ : list
            Evaluations of `f` from domain.
        
        F : float
            Total area under curve, `f`.

        Raises
        ------
        bad_X : string
            If {`X_i`} is neither n x 1 nor 1 x n array.
        
        bad_f : string
            If {`f`} is not an expression.

        Warns
        -----
        func_func : string
            Evaluate input expression for Newton difference approximation.
        
        Notes
        -----
        `X = 0` if not a list nor n x 1 or 1 x n array.

        Unless specified and if `X` is defined, `a` and `b` will be the minimum and maximum, respectively, of `X`.

        Theorem:
        Let f be in C2[a,b], h = (b-a)/n, and xj = a + jh for j = 0, 1, ..., n. There exists a mu in (a,b) for which the quadrature for n sub-intervals can be written with its error term as:
        int_(a)^(b)f(x)dx = h[f(a) + 2*[sum_(j=1)^(n - 1){f(xj)}] + f(b)]/2 - (b-a)*(h^2)f''(mu)/12.

        Where: (b-a)*(h^2)f''(mu)/12 -> O(h^2)
        """
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_f = 'Input range must be expression, not list or tuple.'
        func_func = 'Input expression used.'
        if not isinstance(f,(FunctionType, sp.Expr)):
            if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
            else: sys.exit(bad_f)
        if isinstance(f,(FunctionType, sp.Expr)): print(func_func)
        if X != 0: a, b = min(X), max(X)
        h, a, b = float(h), float(a), float(b)
        XJ, YJ = [], []
        XJ.append(a); YJ.append(f(a))
        j, n, z = 1, math.ceil((b-a)/h), 0
        while j <= n-1:
            x_j = a + j*h
            XJ.append(x_j)
            y_j = f(x_j)
            YJ.append(y_j)
            z += y_j
            j += 1
        XJ.append(b); YJ.append(f(b))
        F = h/2*(f(a) + 2*z + f(b))
        return XJ, YJ, F

    def endpoint(X, Y, h, point_type, which_end):
        """Find the derivative at an endpoint of data set.

        Parameters
        ----------
        X : list
            Domain of collected data.
        
        Y : array or expression
            Range of collected data.
        
        h : float
            Step-size through interval.
        
        point_type : string
            Determines if 3 or 5 pt. method is used.
        
        which_end : string
            Dictates whether evaluated point is left or right most data point.
        
        Returns
        -------
        dY : float
            Evaluated derivative at point.

        Raises
        ------
        bad_X : string
            If {`X`} is neither n x 1 nor 1 x n array.
        
        bad_Y : string
            If {`Y`} is not an expression.
        
        bad_data : string
            If `X` and `Y` are of unequal length.

        See Also
        --------
        make_array() : Prints string that expression was used to make array.

        Notes
        -----
        5 point is more accurate than 3 point; however, round-off error increases.
        """
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_Y = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X_i and Y_i must be of equal length.'
        if not isinstance(Y,(FunctionType, sp.Expr)):
            if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
            if len(Y[0]) > 1 and len(Y[1]) > 1: sys.exit(bad_Y)
            if len(X) != len(Y): sys.exit(bad_data)
        if isinstance(Y,(FunctionType, sp.Expr)): Y = make_array(X, Y)
        h, dY = float(h), 0
        if which_end == 'left':
            i = 0
            if point_type == 'three':
                dY = (-3*Y[i] + 4*Y[i+1] - Y[i+2])/(2*h)
            if point_type == 'five':
                dY = (-25*Y[i] + 48*Y[i+1] \
                    - 36*Y[i+2] + 16*Y[i+3] \
                        - 3*Y[i+4])/(12*h)
        if which_end == 'right':
            i = -1
            if point_type == 'three':
                dY = (-3*Y[i] + 4*Y[i-1] - Y[i-2])/(2*h)
            if point_type == 'five':
                dY = (-25*Y[i] + 48*Y[i-1] \
                    - 36*Y[i-2] + 16*Y[i-3] \
                        - 3*Y[i-4])/(12*h)
        return dY

    def gaussian_legendre(function, a, b):
        return quad(function, a, b)

    def integrate(function, a, b):
        return quad(function, a, b)

    def midpoint(X, Y, h, point_type, i):
        """Find derivative information at some point within data set.

        Parameters
        ----------
        X : list
            Domain of collected data.
        
        Y : array or expression
            Range of collected data.
        
        h : float
            Step-size through interval.
        
        point_type : string
            Determines if 3 or 5 pt. method is used.

        i : int
            Index at which point is to be evaluated.
        
        Returns
        -------
        dY : float
            Evaluated derivative at point.

        Raises
        ------
        bad_X : string
            If {`X`} is neither n x 1 nor 1 x n array.
        
        bad_Y : string
            If {`Y`} is not an expression.
        
        bad_data : string
            If `X` and `Y` are of unequal length.
        
        bad_i : string
            `i` must be an integer and non-zero for indexing.
        
        bad_type : string
            If `point_type` was not an acceptable option.

        See Also
        --------
        make_array() : Prints string that expression was used to make array.

        Notes
        -----
        5 point is more accurate than 3 point; however, round-off error increases.
        """
        bad_X = 'Input domain was neither an n x 1 nor a 1 x n array.'
        bad_Y = 'Input range was neither an n x 1 nor a 1 x n array.'
        bad_data = 'Arrays X_i and Y_i must be of equal length.'
        bad_i = 'Index must be an integer.'
        bad_type = "I am sorry. The selected type was not understood. Please select: 'three', 'five', or '2nd_derivative'."
        if not isinstance(Y,(FunctionType, sp.Expr)):
            if len(X[0]) > 1 and len(X[1]) > 1: sys.exit(bad_X)
            if len(Y[0]) > 1 and len(Y[1]) > 1: sys.exit(bad_Y)
            if len(X) != len(Y): sys.exit(bad_data)
        if isinstance(Y,(FunctionType, sp.Expr)): Y = make_array(X, Y)
        if not isinstance(i,int): sys.exit(bad_i)
        h, dY = float(h), 0
        if point_type == 'three':
            dY = (Y[i+1] - Y[i-1])/(2*h)
        if point_type == 'five':
            dY = (Y[i-2] - 8*Y[i-1] \
                + 8*Y[i+1] - Y[i+2])/(12*h)
        if point_type == '2nd_derivative':
            dY = (Y[i-1] - 2*Y[i] + Y[i+1])/(h**2)
        else: sys.exit(bad_type)
        return dY

    def richard_extrapolation(function, x0, h, order, direction):
        """Results in higher-accuracy of derivative at point in function with lower-order formulas to minimize round-off error and increase O(h) of truncation error.

        Parameters
        ----------
        function : expression
            Polynomial over which derivative must be calculated.
        
        x0 : float
            Point about which extrapolation centers
        
        h : float
            Step-size through interval.

        order : int
            Order for rate of convergence.
        
        direction : string
            `'forward'` or `'backward'` construction.
        
        Returns
        -------
        polynomial : string
            Totally constructed polynomial.
        
        px : float
            Evaluation of `polynomial` at `x`.

        Raises
        ------
        bad_function : string
            If `function` is not an expression.
        
        bad_order : string
            `order` must be an integer and non-zero.
        
        bad_direction : string
            If `direction` is neither `'forward'` nor `'backward'`.
        
        Warns
        -----
        func_func : string
            Evaluate input expression for Newton difference approximation.
        """
        bad_function = 'Function must be expression.'
        func_func = 'Input expression used.'
        bad_order = 'Expected integer.'
        bad_direction = "Supplied direction was not understood. Please specify 'forward' or 'backward'."
        if not isinstance(function,(FunctionType, sp.Expr)): 
            sys.exit(bad_function)
        if isinstance(function,(FunctionType, sp.Expr)):  print(func_func)
        if not isinstance(order,int): sys.exit(bad_order)
        if direction != 'forward' or 'backward': sys.exit(bad_direction)
        def f(h):
            x = x0 + h
            return x, function(x)
        x0, h = float(x0), float(h)
        i, X, FX = 0, [], []
        while i < order:
            dx = h / (2**order) * (2**i)
            x_i, fx_i = f(dx)
            X.append(x_i); FX.append(fx_i)
            i += 1
        m = len(X)
        n = m + 1
        return interpolations.newton_difference(X, FX, x0, direction)
#   #   #   #   #   #   #   #   #


#################################
## End of Code
test.test()     # 'Test complete.'
#   #   #   #   #   #   #   #   #