import math

# Volume contrainte minimale
V0 = 10

# Calcul de la fonction lagrangienne, son gradient et sa hessienne
def lagrangien(r, h, l):
    return (math.pi**2) * (r**2) * ((r**2) + (h**2)) + l*((math.pi/3)* (r**2) * h - V0)

def gradientLagrangien(r, h ,l ):
    x_component = (math.pi**2) * (4* (r**3) + 2* (h**2) *r) + l*((2*math.pi*r*h)/3)
    y_component = 2*h*(math.pi**2)*(r**2) + l*(math.pi/3)*(r**2)
    z_component = (math.pi/3)*(r**2)*h - V0
    return [x_component, y_component, z_component]

def hessianLagrangien(r, h, l):
    # Second partial derivative with respect to r
    d2f_dr2 = (math.pi**2) * (12 * (r**2) + 2 * (h**2)) + l * (2 * math.pi * h / 3)
    # Second partial derivative with respect to h
    d2f_dh2 = 2 * (math.pi**2) * (r**2)
    # Second partial derivative with respect to l (constant term)
    d2f_dl2 = 0
    # Mixed partial derivative with respect to r and h
    d2f_drh = 4 * h * (math.pi**2) + l * (2*math.pi * r / 3)
    # Mixed partial derivative with respect to r and l
    d2f_drl = (2 * math.pi * r * h / 3)
    # Mixed partial derivative with respect to h and l
    d2f_dhl = (math.pi * r**2 / 3)
        
    return [
        [d2f_dr2, d2f_drh, d2f_drl],
        [d2f_drh, d2f_dh2, d2f_dhl],
        [d2f_drl, d2f_dhl, d2f_dl2]
    ]

# Calculating the determinant of a 3x3 matrix
def determinant(matrix):
    return (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
            matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
            matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))

def inverseMatrix(matrix):
    inverse = [[matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1], matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2], matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]],
            [matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2], matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0], matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]],
            [matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0], matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1], matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]]]
            
    det = determinant(matrix)

    for i in range(len(inverse)):
        for j in range(len(inverse[i])):
            inverse[i][j] = inverse[i][j] / det
    
    return inverse

# Calcul du produit scalaire de deux vecteurs
def scalarProduct(v1,v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

# Calcul de la norme d'un vecteur
def norm(x):
    return math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

# Algorithme de Wolfe pour trouver un pas optimal pour un certain X_k et dk
def wolfeStep(X_k, dk):

    lagrange = lagrangien(X_k[0],X_k[1],X_k[2]) # X_k
    gradient = gradientLagrangien(X_k[0],X_k[1],X_k[2]) # X_k

    psk = scalarProduct(dk,gradient)

    c1 = 0.01
    c2 = 0.99

    # Bornes initiales pour le pas, initialisation avec une valeur petite de a pour accélerer la convergence
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

        lagrange_X_k1 = lagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        gradient_X_k1 = gradientLagrangien(X_k1[0],X_k1[1],X_k1[2]) # X_k1
        psk1 = scalarProduct(gradient_X_k1,dk)

        iterations += 1

        # Conditions de Wolfe
        # Condition verifiant que le pas ne soit pas trop grand
        if (lagrange_X_k1 > (lagrange + c1*a*psk)):
            condition1 = 0

            a_max = a

            a= (a_min + a_max)/2
        else:
            condition1 = 1

        # Condition verifiant que le pas ne soit pas trop petit
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