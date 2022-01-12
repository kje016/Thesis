from sage.all import *

var('x')
R = PolynomialRing(GF(2), x)
R.inject_variables()

g = x**24 + x**23 + x**21 + x**20 + x**17 + x**15 + x**13 + x**12 + x**8 + x**4 + x**2 + x + x**0
n = g.degree();
m = (x**2+x+x**0);
if m.degree() < 11:
    m = m*x**(11-m.degree());
k = m.degree()+1;


# More efficient version (from the paper)
D = zero_matrix(k, n);
D[k-1] = [g[n-i] for i in range(1,n)];
t = k-1;
while t >=0:
    for j in range(1, L-1):
        D[t-1, j-1] = Mod(D[t,j] + D[t,0]*g[n-j],2);
        D[t-1,n-1] = D[t, 0]*g[0];
    t = t-1;
print(D)

# CRC bits from polynomial computation
r = (m*x**n).mod(g);
print(Matrix([r[n-1-i] for i in range(n)]));

# CRC bits from matrix computation
m = [m[k-1-i] for i in range(k)];

print(Matrix(GF(2), m)*C)
print(Matrix(GF(2), m)*D)

