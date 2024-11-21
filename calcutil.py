import math

# Volume contrainte minimale, 1e-3 m^3 pour 1 litre de volume
V0 = 1e-3

# Calcul de la fonction lagrangienne, son gradient et sa hessienne
def lagrangien(X_k):
    try:
        r = X_k[0]
        h = X_k[1]
        l = X_k[2]
        return (math.pi**2) * (r**2) * ((r**2) + (h**2)) + l*((math.pi/3)* (r**2) * h - V0)
    except:
        return 0

def gradientLagrangien(X_k):
    try:
        r = X_k[0]
        h = X_k[1]
        l = X_k[2]
        x_component = (math.pi**2) * (4* (r**3) + 2* (h**2) *r) + l*((2*math.pi*r*h)/3)
        y_component = 2*h*(math.pi**2)*(r**2) + l*(math.pi/3)*(r**2)
        z_component = (math.pi/3)*(r**2)*h - V0
        return [x_component, y_component, z_component]
    except:
        return [0, 0, 0]

def hessianLagrangien(X_k):
    r = X_k[0]
    h = X_k[1]
    l = X_k[2]
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
    try:
        return math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    except:
        return 0

def scalarMultiplication(scalar, v):
    return [scalar * a for a in v]

def scalarMatrixMultiplication(scalar, matrix):
    return [[scalar * element for element in row] for row in matrix]

def matrixVectorMultiplication(matrix, vector):
    return [scalarProduct(row, vector) for row in matrix]

def addVectors(v1, v2):
    return [a + b for a, b in zip(v1, v2)]

def subtractVectors(v1, v2):
    return [a - b for a, b in zip(v1, v2)]

def outerProduct(v1, v2):
    return [[a * b for b in v2] for a in v1]

def addMatrices(m1, m2):
    return [[a + b for a, b in zip(row1, row2)] for row1, row2 in zip(m1, m2)]

def subtractMatrices(m1, m2):
    return [[a - b for a, b in zip(row1, row2)] for row1, row2 in zip(m1, m2)]

def multiplyMatrices(m1, m2):
    return [[scalarProduct(row, col) for col in zip(*m2)] for row in m1]

def wolfeStep(X_k, dk):
    alpha = 1.0
    # Paramètres de Wolfe
    c1 = 1e-4
    c2 = 0.9
    max_iterations = 100  # Nombre maximal d'itérations pour trouver le pas optimal

    for _ in range(max_iterations):
        g = gradientLagrangien(X_k)
        X_k1 = addVectors(X_k, scalarMultiplication(alpha, dk))
        # Condition pour verifier si le pas n'est pas trop grand
        if lagrangien(X_k1) > lagrangien(X_k) + c1 * alpha * scalarProduct(g, dk):
            alpha *= 0.5
        else:
            g_k1 = gradientLagrangien(X_k1)
            # Condition pour vérifier si le pas n'est pas trop petit
            if scalarProduct(g_k1, dk) >= c2 * scalarProduct(g, dk):
                return alpha
            alpha *= 2.0

    # Renvoyer 0,00001 si le pas optimal n'est pas trouvé
    return 0.00001