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

def optimalStepGradientMethod(r,h,l):

    alpha = 0.001

    X_k = [r,h,l]

    dk = [10,10,10]
    iterations = 1
    while (norm(dk) > 0.5 and iterations < 10000):

        dk = gradientLagrangien(X_k[0], X_k[1], X_k[2])
        
        X_k = [
            X_k[0] - alpha * dk[0],
            X_k[1] - alpha * dk[1],
            X_k[2] - alpha * dk[2]
            ]
        print (iterations,': ',dk)
        iterations += 1
    print(iterations)
    



def main ():
    r= 1 #random.random()*100
    h= 1 #random.random()*100
    l= -20 #random.random()*100
    fixedStepGradientMethod(r, h, l)
    optimalStepGradientMethod(r, h, l)

main()