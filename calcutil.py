import math

# Volume contrainte minimale
V0 = 1

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

def scalarProduct(v1,v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]