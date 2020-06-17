class Matrix(object):
    def __init__(self, matrix):
        assert len(matrix) > 0 and len(matrix[0]) > 0, 'provided list must be in 2D.'
        self.__matrix = matrix
        self.shape = (len(self.__matrix), len(self.__matrix[0]))

    def __str__(self):
        return str(self.__matrix)

    def __getitem__(self, tup):
        i,j = tup
        return self.__matrix[i][j]

    def __setitem__(self, tup, value):
        i,j = tup
        self.__matrix[i][j] = value

    def __mul__(self, other):
        assert len(self.__matrix[0]) == len(other.__matrix), 'column size of first is not the same as row size of second matrix.'
        
        row = len(self.__matrix)
        col = len(other.__matrix[0])
        ret = [[0 for j in range(col)] for i in range(row)] 

        # O(n^3)
        for i in range(row):
            for j in range(col):
                for k in range(len(self.__matrix[0])):
                    ret[i][j] += self.__matrix[i][k]*other.__matrix[k][j]
        
        return Matrix(ret)
        
def __exchange_row(A, r0):
    if A[r0,r0] != 0:
        return True

    for r1 in range(r0+1,A.shape[0]):
        if A[r1,r0] != 0:
            for i in range(A.shape[1]):
                A[r1,i], A[r0,i] = A[r0,i],A[r1,i]
            return True
    
    return False

def __elimination(A, b, r1, r0):
    if A[r0,r0] == 0:
        return False
    
    mult = A[r1,r0]/A[r0,r0]
    for i in range(r0,A.shape[1]):
        A[r1,i] = A[r1,i]-mult*A[r0,i]

    b[r1,0] = b[r1,0]-mult*b[r0,0]
    return True

def __substitution(A, b):
    ret = [[0] for i in range(b.shape[0])]
    for i in range(b.shape[0]):
        ret[b.shape[0]-1-i][0] = b[b.shape[0]-1-i,0]
        for j in range(b.shape[0]-i,b.shape[0]):
            ret[b.shape[0]-1-i][0] = ret[b.shape[0]-1-i][0]-A[b.shape[0]-1-i,j]*ret[j][0]
        ret[b.shape[0]-1-i][0] /= A[b.shape[0]-1-i,b.shape[0]-1-i]
    
    return Matrix(ret)

def gaussian_elimination(A, b):
    '''
        Ax = b
        given:
            A, b
        return:
            x that satisfies the linear equations
    '''
    assert A.shape[0] == b.shape[0], 'size of A and b are not consistent.'
    
    for r0 in range(A.shape[0]):
        if not __exchange_row(A,r0):
            return False # A is not invertible
        for r1 in range(r0+1,A.shape[0]):
            __elimination(A,b,r1,r0)
    
    return __substitution(A,b)
