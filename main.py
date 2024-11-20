import random
import calcutil as calc

# Méthode de gradient à pas fixe
def fixedStepGradientMethod (r,h,l):

    X_k = [r,h,l]

    # Initialisation arbitraire du vecteur direction
    dk = [10,10,10]
    iterations = 1
    while (calc.norm(dk) > 0.5 and iterations < 10000):

        # Calcul du gradient de la fonction lagrangienne
        dk = calc.gradientLagrangien(X_k)

        # Mise à jour de dk avec le pas fixe
        dk = calc.scalarMultiplication(0.0001, dk)
        
        # Mise à jour de X_k, avec un pas fixe de 0.0001
        X_k = calc.subtractVectors(X_k, dk)
        
        print(iterations, ' | Pas fixe Lagrangien: ', calc.lagrangien(X_k))
        iterations += 1
    return X_k
        
# Méthode de gradient à pas optimal, utilisant le Wolfe step
def optimalStepGradientMethod(r,h,l):

    alpha = 0
    X_k = [r,h,l]

    # Initialisation arbitraire du vecteur direction
    dk = [10,10,10]

    iterations = 1
    while (calc.norm(dk) > 0.8 and iterations < 100000):

        dk = calc.gradientLagrangien(X_k)
        
        # Mise à jour de dk, avec un signe négatif
        dk = calc.scalarMultiplication(-1, dk)

        # Calcul du pas optimal avec la méthode de Wolfe et dk negatif
        alpha = calc.wolfeStep(X_k,dk)
        print('OptimalStep a: ', alpha)

        # Mise à jour de dk avec le pas optimal
        dk = calc.scalarMultiplication(alpha, dk)
        
        # Mise à jour de X_k, NOTE: dk est negatif
        X_k = calc.addVectors(X_k, dk)

        print(iterations, ' | Pas Optimal Lagrangien: ', calc.lagrangien(X_k))
        iterations += 1
    return X_k

# Méthode de Newton
def newtonMethod(r, h, l):
    X_k = [r, h, l]

    # Initialisation arbitraire du vecteur direction
    dk = [10, 10, 10]
    iterations = 1
    
    while calc.norm(dk) > 0.5 and iterations < 50:

        grad = calc.gradientLagrangien(X_k)
        hess = calc.hessianLagrangien(X_k)
        hess_inv = calc.inverseMatrix(hess)
        
        # Calcul de la direction de Newton
        dk = calc.matrixVectorMultiplication(hess_inv, grad)

        # Mise à jour de X_k
        X_k = calc.subtractVectors(X_k, dk)
            
        print(iterations, ' | Newton Lagrangien: ', calc.lagrangien(X_k))
        iterations += 1
    print(iterations)
    return X_k

def bfgs(r, h, l):

    # Initialisation of initial guess
    X_k = [r,h,l] 
    tolerance = 1e-5

    # Initialisation of Hessian approximation, Identity matrix
    H = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    
    for iterations in range(1, 10000):
        grad = calc.gradientLagrangien(X_k)
        
        # Check convergence
        grad_norm = calc.scalarProduct(grad, grad) ** 0.5
        if grad_norm < tolerance:
            return X_k
        
        # Compute search direction
        dk = calc.scalarMultiplication(-1, calc.matrixVectorMultiplication(H, grad))
        
        # Line search to find optimal step size
        alpha = calc.wolfeStep(X_k, dk)
        
        # Update x
        X_k1 = calc.addVectors(X_k, calc.scalarMultiplication(alpha, dk))
        s = calc.subtractVectors(X_k1, X_k) # Calculation of s by the difference of X_k1 and X_k
        X_k = X_k1
        
        # Update gradient and compute y
        grad_k1 = calc.gradientLagrangien(X_k1)
        y = calc.subtractVectors(grad_k1, grad)
        
        # Update Hessian approximation (BFGS formula)
        if calc.scalarProduct(y, s) > 0:
            rho = 1.0 / calc.scalarProduct(y, s)
            sy_outer = calc.outerProduct(s, y)
            ss_outer = calc.outerProduct(s, s)
            H = calc.subtractMatrices(calc.subtractMatrices(H, calc.scalarMatrixMultiplication(rho, calc.multiplyMatrices(sy_outer, H))), calc.scalarMatrixMultiplication(rho, calc.multiplyMatrices(H, sy_outer)))
            H = calc.addMatrices(H, calc.scalarMatrixMultiplication(rho, ss_outer))
        else:
            H = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        
        print(iterations, ' | BFGS Lagrangien: ', calc.lagrangien(X_k))
    
    # If max iterations reached, return current estimate
    return X_k, iterations

def main ():
    r = 1#random.randint(0, 100)
    h = 1#random.randint(0, 100)
    l = 20#random.randint(0, 100)
    # print('Initial values: ', r, h, l)
    # print('Fixed step gradient method: ', fixedStepGradientMethod(r,h,l))
    # print('Optimal step gradient method: ', optimalStepGradientMethod(r,h,l))
    # print('Newton method: ', newtonMethod(r,h,l))
    print('BFGS method: ', bfgs(r,h,l))
main()