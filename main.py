import random
import calcutil as calc

# Méthode de gradient à pas fixe
def fixedStepGradientMethod (r,h,l):

    X_k = [r,h,l]

    # Initialisation arbitraire du vecteur direction
    dk = [10,10,10]
    iterations = 1
    while (calc.norm(dk) > 0.5 and iterations < 1000):

        # Calcul du gradient de la fonction lagrangienne
        dk = calc.gradientLagrangien(X_k)
    
        # Mise à jour de dk avec le pas fixe
        dk = calc.scalarMultiplication(0.0001, dk)
        
        # Mise à jour de X_k, avec un pas fixe de 0.0001
        X_k = calc.subtractVectors(X_k, dk)
        
        print(iterations, ' | Pas fixe Lagrangien: ', calc.lagrangien(X_k))
        iterations += 1
    return X_k, iterations
        
# Méthode de gradient à pas optimal, utilisant le Wolfe step
def optimalStepGradientMethod(r,h,l):

    alpha = 0
    X_k = [r,h,l]

    # Initialisation arbitraire du vecteur direction
    dk = [10,10,10]

    iterations = 1
    while (calc.norm(dk) > 0.8 and iterations < 1000):

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
    return X_k, iterations

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
    return X_k, iterations

def bfgs(r, h, l):

    # Initialisation of initial guess
    X_k = [r,h,l] 
    tolerance = 5

    # Initialisation of Hessian approximation, Identity matrix
    H = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    
    for iterations in range(1, 100):
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

        if abs(calc.lagrangien(X_k1) - calc.lagrangien(X_k)) < 0.00001 and iterations > 10:
            break
        
        s = calc.subtractVectors(X_k1, X_k) # Calculation of s by the difference of X_k1 and X_k
        X_k = X_k1

        
        # Update gradient and compute y
        grad_k1 = calc.gradientLagrangien(X_k1)
        y = calc.subtractVectors(grad_k1, grad)
        
        # Update Hessian approximation (BFGS formula)
        if calc.scalarProduct(y, s) != 0:
            rho = 1 / calc.scalarProduct(y, s)
            sy_outer = calc.outerProduct(s, y)
            ss_outer = calc.outerProduct(s, s)
            H = calc.subtractMatrices(calc.subtractMatrices(H, calc.scalarMatrixMultiplication(rho, calc.multiplyMatrices(sy_outer, H))), calc.scalarMatrixMultiplication(rho, calc.multiplyMatrices(H, sy_outer)))
            H = calc.addMatrices(H, calc.scalarMatrixMultiplication(rho, ss_outer))
        else:
            H = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Reinitialiser la matrice H si le produit scalaire est 0
        
        print(iterations, ' | BFGS Lagrangien: ', calc.lagrangien(X_k))
    
    # If max iterations reached, return current estimate
    return X_k, iterations



def main ():
    r = random.randint(0, 100)
    h = random.randint(0, 100)
    l = random.randint(0, 100)

    fixedStepResults, i1 = fixedStepGradientMethod(r,h,l)
    optimalStepResults, i2 = optimalStepGradientMethod(r,h,l)
    newtonResults, i3 = newtonMethod(r,h,l)
    bfgsResults, i4 = bfgs(r,h,l)

    print('Fixed step gradient method: ', fixedStepResults, ' | Iterations: ', i1)
    print('Optimal step gradient method: ', optimalStepResults, ' | Iterations: ', i2)
    print('Newton method: ', newtonResults, ' | Iterations: ', i3)
    print('BFGS method: ', bfgsResults, ' | Iterations: ', i4)
main()