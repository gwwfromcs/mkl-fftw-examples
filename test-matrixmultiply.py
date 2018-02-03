import numpy as np

m=10; k=5; n=7;
a = np.zeros((m,k))
b = np.zeros((k,n))

for i in range(m):
    for j in range(k):
        a[i,j] = float(i+1)/float(j+1)

for i in range(k):
    for j in range(n):
        b[i,j] = float(i+1)/float(j+1)

print "a = \n", a
print "b = \n", b
print "c = a * b"
print np.dot(a,b)
print "d = a^T * a"
print np.dot(a.transpose(), a)
print "e = a * a^T"
print np.dot(a, a.transpose())
