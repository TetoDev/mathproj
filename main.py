import random
import calcutil as calc

usedBFGS = False

def fixedStepGradientMethod (r,h,l):

    X_k = [r,h,l]

    dk = [10,10,10]
    iterations = 1
    while (calc.norm(dk) > 0.5 and iterations < 10000):

        dk = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        
        X_k = [
            X_k[0] - 0.0001 * dk[0],
            X_k[1] - 0.0001 * dk[1],
            X_k[2] - 0.0001 * dk[2]
            ]
        print(iterations, ' | Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
        iterations += 1
    return X_k


def wolfeStep(X_k, dk):

    lagrange = calc.lagrangien(X_k[0],X_k[1],X_k[2]) # X_k
    gradient = calc.gradientLagrangien(X_k[0],X_k[1],X_k[2]) # X_k

    psk = calc.scalarProduct(dk,gradient)

    c1 = 0.01
    c2 = 0.99

    a_min = 0
    a_max = 100
    a =1.0E-06
    condition1 = 0
    condition2 = 0
    iterations = 1 
    while (((condition1 + condition2) < 2) and (iterations < 100)):

        X_k1 = [
            X_k[0] + a*dk[0],
            X_k[1] + a*dk[1],
            X_k[2] + a*dk[2],
        ]

        lagrange_X_k1 = calc.lagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        gradient_X_k1 = calc.gradientLagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        psk1 = calc.scalarProduct(gradient_X_k1,dk)

        iterations += 1

        if (lagrange_X_k1 > (lagrange + c1*a*psk)):
            condition1 = 0

            a_max = a

            a= (a_min + a_max)/2
        else:
            condition1 = 1

        if (-psk1 > -c2*psk) :
            condition2 = 0
            
            if a_max < 100:
                a_min = a

                a = (a_min+a_max)/2
            else :
                a = 2*a
        else :
            condition2 = 1

    return a
        

def optimalStepGradientMethod(r,h,l):

    alpha = 0
    X_k = [r,h,l]
    dk = [10,10,10]

    
    iterations = 1
    while (calc.norm(dk) > 0.5 and iterations < 100000):

        dk = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        for i in range(len(dk)):
            dk[i] = -dk[i]

        alpha = wolfeStep(X_k,dk)
        print('OptimalStep a: ', alpha)
        
        
        X_k = [
            X_k[0] + alpha * dk[0],
            X_k[1] + alpha * dk[1],
            X_k[2] + alpha * dk[2]
            ]

        print(iterations, ' | Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
        iterations += 1
    return X_k

def newtonMethod(r, h, l):
    X_k = [r, h, l]
    dk = [10, 10, 10]
    iterations = 1
    
    while calc.norm(dk) > 0.5 and iterations < 50:

        grad = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        hess = calc.hessianLagrangien(X_k[0], X_k[1], X_k[2])
        hess_inv = calc.inverseMatrix(hess)
            
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

def quasiNewtonMethod(r, h, l):
        global usedBFGS
    # Initial guess
        X_k = [r, h, l]
        iterations = 0
        max_iterations = 1000
    
        # Initial Hessian approximation (identity matrix)
        H_k = [[1,0,0],[0,1,0],[0,0,1]]
    
        while iterations < max_iterations:

            # Actual hessian of lagrange for current X_k
            hess_inv = calc.inverseMatrix(calc.hessianLagrangien(X_k[0], X_k[1], X_k[2]))
    
            grad = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])

            # Compute search direction
            dk = [-sum(H_k[i][j] * grad[j] for j in range(3)) for i in range(3)]

            # Compute newton direction
            newton_dk = [
                hess_inv[0][0] * grad[0] + hess_inv[0][1] * grad[1] + hess_inv[0][2] * grad[2],
                hess_inv[1][0] * grad[0] + hess_inv[1][1] * grad[1] + hess_inv[1][2] * grad[2],
                hess_inv[2][0] * grad[0] + hess_inv[2][1] * grad[1] + hess_inv[2][2] * grad[2]
            ]

            # Check if newton direction is a descent direction
            condition = calc.scalarProduct(newton_dk, grad) >= 0
            if condition:
                break
            usedBFGS = True
            
            # Find alpha using Wolfe conditions, switch to newton if values are too big
            try:
                alpha = wolfeStep(X_k,dk)
            except:
                print('No alpha found')
                break

            # Set Xk+1
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
            estimation = calc.scalarProduct(s_k, y_k)
            if estimation != 0:
                rho_k = 1.0 / calc.scalarProduct(y_k, s_k)
            else:
                rho_k = 0

            V = [[1 - rho_k * s_k[i] * y_k[j] for j in range(3)] for i in range(3)]
            H_k = [[sum(V[i][k] * H_k[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            H_k = [[H_k[i][j] + rho_k * s_k[i] * s_k[j] for j in range(3)] for i in range(3)]

    
            # Update X_k and gradient
            X_k = X_k1
            grad = grad1
    
            print(iterations, ' | BFGS Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
            iterations += 1
        
        
    
        print(iterations)
        return newtonMethod(X_k[0], X_k[1], X_k[2])
    



def main ():
    global usedBFGS
    #solution = fixedStepGradientMethod(r, h, l)
    #solution = optimalStepGradientMethod(r, h, l)
    # solution = newtonMethod(r, h, l)
    while not usedBFGS:
        r= random.random()*100
        h= random.random()*100
        l= random.random()*100
        print(r,h,l)
        hessian = calc.hessianLagrangien(r,h,l)
        print('Hessian')
        for i in range(len(hessian)):
            print(hessian[i])
        print()
        print('Inverse Hessian')
        inverse = calc.inverseMatrix(hessian)
        for i in range(len(inverse)):
            print(inverse[i])
        solution = quasiNewtonMethod(r, h, l)
        print(solution)
main()