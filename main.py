import math, random
import calcutil as calc

# Volume contrainte minimale
V0 = 10


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
        print (iterations,': ',dk)
        iterations += 1
    print(iterations)


def wolfeStep(X_k, dk):

    lagrange = calc.lagrangien(X_k[0],X_k[1],X_k[1]) # X_k
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
        
        print(iterations,'. a: ', a)


        X_k1[0] = [
            X_k[0] + a*dk[0],
            X_k[1] + a*dk[1],
            X_k[2] + a*dk[2],
        ]

        lagrange_X_k1 = calc.lagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        gradient_X_k1 = calc.gradientLagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        psk1 = calc.scalarProduct(gradient_X_k1,dk)

        iterations += 1

        #print(lagrange_X_k1, ' > ', lagrange + c1*a*psk)
        if (lagrange_X_k1 > (lagrange + c1*a*psk)):
            condition1 = 0

            a_max = a

            a= (a_min + a_max)/2
        else:
            condition1 = 1

        #print(-psk1 , ' > ', -c2*psk)
        if (-psk1 > -c2*psk) :
            condition2 = 0
            
            if a_max < 100:
                a_min = a

                a = (a_min+a_max)/2
            else :
                a = 2*a
        else :
            condition2 = 1
        
        print(condition1, condition2)

    
    return a
        

def optimalStepGradientMethod(r,h,l):

    alpha = 0
    X_k = [r,h,l]
    dk = [10,10,10]

    
    iterations = 1
    while (calc.norm(dk) > 0.5 and iterations < 100):

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

        print (iterations,' Lagrangien: ', calc.lagrangien(X_k[0], X_k[1], X_k[2]))
        iterations += 1

def newtonMethod(r, h, l):
    pass

def quasiNewtonMethod(r, h, l):
    X_k = [r, h, l]
    dk = [10, 10, 10]
    iterations = 1
    
    while calc.norm(dk) > 0.5 and iterations < 10000:

        grad = calc.gradientLagrangien(X_k[0], X_k[1], X_k[2])
        hess = calc.hessianLagrangien(X_k[0], X_k[1], X_k[2])
        hess_inv = calc.inverseMatrix(hess)
            
        dk = [
            # I think this is wrong
            hess_inv[0][0] * grad[0] + hess_inv[0][1] * grad[1] + hess_inv[0][2] * grad[2],
            hess_inv[1][0] * grad[0] + hess_inv[1][1] * grad[1] + hess_inv[1][2] * grad[2],
            hess_inv[2][0] * grad[0] + hess_inv[2][1] * grad[1] + hess_inv[2][2] * grad[2]
        ]
            
        X_k = [
            X_k[0] - dk[0],
            X_k[1] - dk[1],
            X_k[2] - dk[2]
        ]
            
        print(iterations, ': ', dk)
        iterations += 1
    print(iterations)


    



def main ():
    r= 1 #random.random()*100
    h= 1 #random.random()*100
    l= -20 #random.random()*100
    #fixedStepGradientMethod(r, h, l)
    optimalStepGradientMethod(r, h, l)
    #quasiNewtonMethod(r, h, l)

main()