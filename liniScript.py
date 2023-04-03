
def gcd(a,b):
    if b == 0:
        return a
    
    if abs(b) > abs(a):
        return gcd(b, a)
    return gcd(b, a % b)

def isPrime(p):
    for i in range(1,p):
        if p % i == 0:
            return False
    return True

def inverse(a, p):
    return InverseRec(a+p, p, 1, 0)
    
def InverseRec(a, b, s, t):
    if b == 0:
        return s
    return InverseRec(b, a%b, t, s - (a//b)*t)

class poly():
    coefs = []
    field = ""
    
    def __init__(self, field, coefs=[0]):
        if type(coefs) == poly:
            coefs = coefs.coefs
        if type(coefs) != list:
            coefs = [coefs]
        self.field = field
        for i in range(len(coefs)):
            self[i] = coefs[i]
            self.deg += 1
        self.deg()
    
    def __len__(self):
        return len(self.coefs)
    
    def __getitem__(self, key):
        return self.coefs[key]
    
    def __setitem__(self, key, item):
        if key >= len(self):
            while len(self) <= key:
                self.coefs.append(scalar(self.field, 0))
        self.coefs[key].setValue(item)
        self.deg()
    
    def __eq__(self, other):
        other = poly(other)
        self.deg()
        other.minimizeDeg()
        return self.coefs == other.coefs
    
    def __neg__(self):
        return poly(self.field, coefs=[-self[i] for i in range(len(self))])
    
    def __add__(self, other):
        result = poly(self.field, coefs=[self[i] + other[i] for i in range(len(self))])
        result.deg()
        return result
     def __sub__(self, other):
        return self + (-other)
    
    def __mul__(self, other):
        
        
    
    def deg(self):
        for i in range(len(self) - 1, 0):
            if self[i] == 0:
                self.coefs.pop()
        return len(self) - 1
    
class scalar():
    field = ""
    val = 0
    def __init__(self, f, v):
        if f[0] == 'Z':
            if not isPrime(int(f[1:])):
                raise Exeption("finite ring is not a field")
        self.field = f
        self.setVal(v)
        
    def __repr__(self):
        self.less()
        sym = {'C' : " + i", 'Q' : " / "}
        if self.field[0] in "RZ":
            return str(self.val)
        if self.field in "CQ":
            return str(self.val[0]) + sym[self.field] + str(self.val[1])
    
    def __eq__(self, other):
        other = scalar(self.field, other)
        other.less()
        self.less()
        return self.val == other.val
    
    def __neg__(self):
        result = scalar(self.field, 0)
        if self.field == "C":
            result.val[0] = -self.val[0]
            result.val[1] = -self.val[1]
        if self.field == "Q":
            result.val[0] = self.val[0]
            result.val[1] = -self.val[1]
        else:
            result.val = -self.val
        result.less()
        return result
    
    def __add__(self, other):
        other = scalar(self.field, other)
        result = scalar(self.field, 0)
        
        if self.field in "CQ":
            a, b, c, d = self.val[0], self.val[1], other.val[0], other.val[1]
            if self.field == "Q":
                result.val = [a*d + b*c, b*d]
            else:
                result.val = [a+c, b+d]
        else:
            result.val = self.val + other.val
        result.less()
        return result
    
    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        other = scalar(self.field, other)
        result = scalar(self.field, 0)
        
        if self.field in "CQ":
            a, b, c, d = self.val[0], self.val[1], other.val[0], other.val[1]
            if self.field == "Q":
                result.val = [a*d - b*c, b*d]
            else:
                result.val = [a-c, b-d]
        else:
            result.val = self.val - other.val
        result.less()
        return result
    
    def __rsub__(self, other):
        return -(self - other)
    
    def __mul__(self, other):
        if type(other) == scalar:
            other = scalar(self.field, other)
            result = scalar(self.field, 0)    
            if self.field in "CQ":
                a, b, c, d = self.val[0], self.val[1], other.val[0], other.val[1]
                if self.field == "Q":
                    result.val = [a*c, b*d]
                else:
                    result.val = [a*c - b*d, a*d + b*c]
            else:
                result.val = self.val * other.val
            result.less()
            return result
    
    def __rmul__(self, other):
        return self * other
    
    def __truediv__(self, other):
        other = scalar(self.field, other)
        result = scalar(self.field, 0)
        
        if self.field in "CQ":
            a, b, c, d = self.val[0], self.val[1], other.val[0], other.val[1]
            if self.field == "Q":
                result.val = [a*d, b*c]
            else:
                result.val = [(a*c + b*d) / (c*c + d*d), (-a*d + b*c) / (c*c + d*d)]
        elif self.field == "R":
            result.val = self.val / other.val
        else:
            self.val * inverse(other.val, int(other.field[1:]))
        result.less()
        return result
    
    def __rtruediv__(self, other):
        I = scalar(self.field, 1)
        return I / (self / other)
    
    def setVal(self, v):
        if type(v) == scalar:
            v = v.val
            
        elif type(v) == int:
            if self.field == "Q":
                v = [v, 1]
            if self.field == "C":
                v = [v, 0]
        
        elif type(v) == float:
            if type.field[0] in "ZQ":
                raise Exeption("cant convert to rational or integer")
            if type.filed == "C":
                v = [v, 0]

        elif type(v) == str:
            if self.field[0] == "Z":
                v = int(v)
            if self.field == "Q":
                try:
                    v = [int(v), 1]
                except:
                    v = v.split(' / ')
                    v[0] = int(v[0])
                    v[1] = int(v[1])
            if self.field == "C":
                v = v.split(' + i')
                v[0] = float(v[0])
                v[1] = float(v[1])
            if  self.field == "R":
                v = float(v)
        self.val = v
        self.less
    

    
    def less(self):
        if self.field == "Q":
            d = gcd(self.val[0], self.val[1])
            self.val[0] //= d
            self.val[1] //= d
        
        if self.field[0] == 'Z':
            self.val %= int(self.field[1:])
    
    class vector():
        val = []
        field = ""
        
        def __init__(self, field, val):
            self.field = field
            for v in val:
            self.val = val
        
        def __len__
        
        def __getitem__(self, key):
            return self.val[key]
        
        def 
    
    
    
class matrix():
    m = 0
    n = 0
    mat =  []
    field = ""
    eigenValues
    
    def __init__(self, m, n, field):
        self.mat = [[scalar(field, 0) for j in range(n)] for i in range(m)]
        self.m = m
        self.n = n
        self.field = field
    
    def __repr__(self):
        cell = 2 + max([max([len(str(self[(i, j)])) for j in range(self.n)]) for i in range(self.m)])
        result = ""
        s = ""
        l = 0
        for i in range(self.m):
            result += '|'
            for j in range(self.n):
                s = str(self[(i, j)])
                l = cell - len(s)
                result += (l // 2) * " "
                result += s
                result += (l - (l // 2)) * " "
            result += "|\n"
        return result
                
                
    
    def __getitem__(self, key):
        return self.mat[key[0]][key[1]]
    
    def __setitem__(self, key, item):
        self.mat[key[0]][key[1]].setVal(item)
    
    def __eq__(self, other):
        return self.mat == other.mat
    
    def __neg__(self):
        result = matrix(self.m, self.n, slef.field)
        for i in range(self.m):
            for j in range(self.n):
                result[(i, j)] = -self[(i, j)]
        return result
    
    def __add__(self, other):
        if self.m != other.m or self.n != other.m:
            raise Exeption("dimentions not right for addition")
        result = matrix(self.m, self.n, self.field)
        for i in range(self.m):
            for j in range(self.n):
                result[(i, j)] = self[(i, j)] + other[(i, j)]
        return result
        
    def __sub__(self, other):
        if self.m != other.m or self.n != other.m:
            raise Exeption("dimentions not right for subtraction")
        result = matrix(self.m, self.n, self.field)
        for i in range(self.m):
            for j in range(self.n):
                result[(i, j)] = self[(i, j)] - other[(i, j)]
        return result
        
    def __mul__(self, other):
        if type(other) != matrix:
            other = scalar(self.field, other)
            B = matrix(self.m, self.m, self.field)
            for i in range(self.m):
                B[(i, i)] = other
        if type(other) == matrix:
            B = other
        result = matrix(self.m, B.n, self.field)
        if self.n != B.m:
            raise Exeption("dimentions not right for multiplication")
        for i in range(self.m):
            for j in range(self.n):
                for k in range(self.n):
                    result[(i, j)] = result[(i, j)] + ( self[(i, k)] * B[(k, j)] )
        return result
    
    def __rmul__(self, other):
        return self * other
    
    def input(self):
        for i in range(self.m):
            for j in range(self.n):
                self[(i, j)] = input("({0}, {1}): ".format(i, j))
    
    def copy(self):
        result = matrix(self.m, self.n, self.field)
        for i in range(self.m):
            for j in range(self.n):
                result[(i, j)] = self[(i, j)]
        return result
        
    def P1(self, i, k):
        temp = scalar(self.field, 0)
        for j in range(self.m):
            temp.setVal(self[(i, j)])
            self[(i, j)] = self[(k, j)]
            self[(k, j)] = temp
    
    def P2(self, i, L):
        for j in range(self.m):
            self[(i, j)] = L * self[(i, j)]
    
    def P3(self, i, k, L):
        for j in range(self.m):
            self[(i, j)] = self[(i, j)] + ( L * self[(k, j)] )
    
    def Derug(self, other):
        if other.m != self.m:
                raise Exeption("hights dont match")
        A = self.copy()
        det = scalar(self.field, 1)
        r = 0
        c = 0
        while r < self.m and c < self.n:
            for i in range(r, self.m):
                if not A[(i, c)] == 0:
                    if r != i:
                        other.P1(r, i)
                        A.P1(r, i)
                        det = -det
                    other.P2(r, 1/A[(r, c)])
                    A.P2(r, 1/A[(r, c)])
                    det = det / A[(r, c)]
                    for k in range(self.m):
                        if k != r:
                            other.P3(k, r, -A[(k, c)])
                            A.P3(k, r, -A[(k, c)])
                    r += 1
                    break
            c += 1
        if r != self.m:
            det = 0
        return r, det, A, other
    
    def rank(self):
        Z = matrix(self.m, self.m, self.field)
        return self.Derug(Z)[0]
        
    def det(self):
        if self.m != self.n:
            raise exeption("matrix is not squere")
        
        Z = matrix(self.m, self.m, self.field)
        return self.Derug(Z)[1]
    
    def echelon(self):
        Z = matrix(self.m, self.m, self.field)
        return self.Derug(Z)[2]
    
    def inverse(self):
        if self.m != self.n:
            raise exeption("matrix is not squere")
        
        I = matrix(self.m, self.m, self.field)
        for i in range(self.m):
            I[(i, i)] = 1
        return self.Derug(I)[3]
        
    