import math, random

# Volume contrainte minimale
V0 = 10

def lagrangien(r, h, l):
    return (math.pi**2) * (r**2) * ((r**2) + (h**2)) + l*((math.pi/3)* (r**2) * h - V0)

def norm(x):
    return math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

def gradientLagrangien(r, h ,l ):
    x_component = (math.pi**2) * (4* (r**3) + 2* (h**2) *r) + l*((2*math.pi*r*h)/3)
    y_component = 2*h*(math.pi**2)*(r**2) + l*(math.pi/3)*(r**2)
    z_component = (math.pi/3)*(r**2)*h - V0
    return [x_component, y_component, z_component]

def hessianLagrangien(r, h, l):
    # Second partial derivative with respect to r
    d2f_dr2 = 1
    # Second partial derivative with respect to h
    d2f_dh2 = 1
    # Second partial derivative with respect to l (constant term)
    d2f_dl2 = 1
    # Mixed partial derivative with respect to r and h
    d2f_drh = 2
    # Mixed partial derivative with respect to r and l
    d2f_drl = 2
    # Mixed partial derivative with respect to h and l
    d2f_dhl = 5
        
    return [
        [d2f_dr2, d2f_drh, d2f_drl],
        [d2f_drh, d2f_dh2, d2f_dhl],
        [d2f_drl, d2f_dhl, d2f_dl2]
    ]

def inverseMatrix(matrix):
    inverse = [[matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1], matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2], matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]],
            [matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2], matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0], matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]],
            [matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0], matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1], matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]]]
            
    det = determinant(matrix)

    for i in range(len(inverse)):
        for j in range(len(inverse[i])):
            inverse[i][j] = inverse[i][j] / det
    
    return inverse

def determinant(matrix):
    # Calculating the determinant of a 3x3 matrix
    return (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
            matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
            matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))

def fixedStepGradientMethod (r,h,l):

    X_k = [r,h,l]

    dk = [10,10,10]
    iterations = 1
    while (norm(dk) > 0.5 and iterations < 10000):

        dk = gradientLagrangien(X_k[0], X_k[1], X_k[2])
        
        X_k = [
            X_k[0] - 0.001 * dk[0],
            X_k[1] - 0.001 * dk[1],
            X_k[2] - 0.001 * dk[2]
            ]
        print (iterations,': ',dk)
        iterations += 1
    print(iterations)

def scalarProduct(v1,v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def wolfeStep(X_k, dk):

    lagrange = lagrangien(X_k[0],X_k[1],X_k[1]) # X_k
    gradient = gradientLagrangien(X_k[0],X_k[1],X_k[2]) # X_k

    psk = scalarProduct(dk,gradient)

    c1 = 0.1
    c2 = 0.99

    a_min = 0
    a_max = 100
    a = (a_min + a_max)/2

    tooBig = True
    tooSmall = True
    iterations = 1 
    while ((tooSmall or tooBig) and (iterations < 100)):
        
        #print(iterations,'. a: ', a)


        X_k1 = [
            X_k[0] + a*dk[0],
            X_k[1] + a*dk[1],
            X_k[2] + a*dk[2],
        ]

        lagrange_X_k1 = lagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        gradient_X_k1 = gradientLagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        psk1 = scalarProduct(gradient_X_k1,dk)

        iterations += 1

        print(lagrange_X_k1, ' > ', lagrange + c1*a*psk)
        if (lagrange_X_k1 > (lagrange + c1*a*psk)):
            tooBig = True

            a_max = a

            a= (a_min + a_max)/2
        else:
            tooBig = False

        if (-psk1 > -c2*psk) :
            tooSmall = True
            
            a_min = a

            a= (a_min+a_max)/2
        else :
            tooSmall = False

    
    return a
        
        


def optimalStepGradientMethod(r,h,l):

    alpha = 0
    X_k = [r,h,l]
    dk = [10,10,10]

    
    iterations = 1
    while (norm(dk) > 0.5 and iterations < 100):

        dk = gradientLagrangien(X_k[0], X_k[1], X_k[2])
        alpha = wolfeStep(X_k,dk)
        print('OptimalStep a: ', alpha)
        
        
        X_k = [
            X_k[0] - alpha * dk[0],
            X_k[1] - alpha * dk[1],
            X_k[2] - alpha * dk[2]
            ]
        print (iterations,': ',dk)
        iterations += 1
    print(iterations)

def newtonMethod(r, h, l):
    pass

def quasiNewtonMethod(r, h, l):
    X_k = [r, h, l]
    dk = [10, 10, 10]
    iterations = 1
    
    while norm(dk) > 0.5 and iterations < 10000:

        grad = gradientLagrangien(X_k[0], X_k[1], X_k[2])
        hess = hessianLagrangien(X_k[0], X_k[1], X_k[2])
        hess_inv = inverseMatrix(hess)
            
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