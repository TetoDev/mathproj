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
        dk = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        
        # Mise à jour de X_k, avec un pas fixe de 0.0001
        X_k = [
            X_k[0] - 0.0001 * dk[0],
            X_k[1] - 0.0001 * dk[1],
            X_k[2] - 0.0001 * dk[2]
            ]
        
        print(iterations, ' | Pas fixe Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
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

        dk = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        for i in range(len(dk)):
            dk[i] = -dk[i]

        # Calcul du pas optimal avec la méthode de Wolfe
        alpha = calc.wolfeStep(X_k,dk)
        print('OptimalStep a: ', alpha)
        
        # Mise à jour de X_k
        X_k = [
            X_k[0] + alpha * dk[0],
            X_k[1] + alpha * dk[1],
            X_k[2] + alpha * dk[2]
            ]

        print(iterations, ' | Pas Optimal Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
        iterations += 1
    return X_k

# Méthode de Newton
def newtonMethod(r, h, l):
    X_k = [r, h, l]

    # Initialisation arbitraire du vecteur direction
    dk = [10, 10, 10]
    iterations = 1
    
    while calc.norm(dk) > 0.5 and iterations < 50:

        grad = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        hess = calc.hessianLagrangien(X_k[0], X_k[1], X_k[2])
        hess_inv = calc.inverseMatrix(hess)
        
        # Calcul de la direction de Newton
        dk = [
            hess_inv[0][0] * grad[0] + hess_inv[0][1] * grad[1] + hess_inv[0][2] * grad[2],
            hess_inv[1][0] * grad[0] + hess_inv[1][1] * grad[1] + hess_inv[1][2] * grad[2],
            hess_inv[2][0] * grad[0] + hess_inv[2][1] * grad[1] + hess_inv[2][2] * grad[2]
        ]
            
        X_k = [
            X_k[0] - dk[0],
            X_k[1] - dk[1],
            X_k[2] - dk[2]
        ]
            
        print(iterations, ' | Newton Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
        iterations += 1
    print(iterations)
    return X_k

def bfgs(r, h, l):
    # Initial guess
    X_k = [r, h, l]
    iterations = 0
    max_iterations = 1000
    tolerance = 1e-5

    # Initial Hessian approximation (identity matrix)
    H_k = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    while iterations < max_iterations:
            
            
            grad = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])

            # Compute search direction
            dk = [-sum(H_k[i][j] * grad[j] for j in range(3)) for i in range(3)]

            if calc.norm(dk) < tolerance:
                break
    
            alpha = calc.wolfeStep(X_k,dk)
    
            # Update X_k
            X_k1 = [
                X_k[0] + alpha*dk[0],
                X_k[1] + alpha*dk[1],
                X_k[2] + alpha*dk[2],
            ]
    
            # Compute new gradient
            grad1 = calc.gradientLagrangien(X_k1[0], X_k1[1], X_k1[2])
    
            # Compute s_k and y_k
            s_k = [
                X_k1[0] - X_k[0],
                X_k1[1] - X_k[1],
                X_k1[2] - X_k[2]
            ]
            y_k = [
                grad1[0] - grad[0],
                grad1[1] - grad[1],
                grad1[1] - grad[2]
            ]
    
            # Update Hessian approximation using BFGS formula
            correction = calc.scalarProduct(s_k, y_k)
            if correction != 0:
                rho_k = 1.0 / correction
            else:
                rho_k = 0.0

            outer_sk_yk = [[s_k[i] * y_k[j] for j in range(3)] for i in range(3)]
            outer_sk_sk = [[s_k[i] * s_k[j] for j in range(3)] for i in range(3)]
            H_k = [[H_k[i][j] - rho_k * sum(H_k[i][k] * outer_sk_yk[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            H_k = [[H_k[i][j] - rho_k * sum(outer_sk_yk[i][k] * H_k[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            H_k = [[H_k[i][j] + rho_k * outer_sk_sk[i][j] for j in range(3)] for i in range(3)]

    
            # Update X_k and gradient
            X_k = X_k1
            grad = grad1
    
            print(iterations, ' | BFGS Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
            iterations += 1
            
    return X_k

def main ():
    r = 1#random.randint(0, 100)
    h = 1#random.randint(0, 100)
    l = -20#random.randint(0, 100)
    # print('Initial values: ', r, h, l)
    # print('Fixed step gradient method: ', fixedStepGradientMethod(r,h,l))
    # print('Optimal step gradient method: ', optimalStepGradientMethod(r,h,l))
    # print('Newton method: ', newtonMethod(r,h,l))
    print('BFGS method: ', bfgs(r,h,l))
main()