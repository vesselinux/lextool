# 
# Original file: ks-compute-2x2x4.sage
# 

# key schedule for standard lexs-2x2x4
BITS = 4
ROWS = 2
COLS = 2
KEYS = 10
NSTATE = ROWS*COLS*BITS
Nk = NSTATE*KEYS
N=Nk

#load "leaks-2x2x4.sage"

P=BooleanPolynomialRing(N,'x',order='lex')

# keys
k = []
for r in range(0,KEYS):
    k.append([])
    for i in range((NSTATE*r)+0,(NSTATE*r)+NSTATE):
        k[r].append(P.gen(i))

# field of fractions
R.<z> = P.fraction_field()[]
#rijndael polynomial for GF16
rp = z^4+z+1

# SubBytes4() step 1 - inverse transformation; calculate inverses mod 16
def inverse4(a):

    #print "a",a
    #a0 = a[0]; a1 = a[1]; a2 = a[2]; a3 = a[3]
    apoly = 0

    # represent a as a polynomial in z
    for i in range(0,BITS):
        apoly+=(a[i])*z^i
        #apoly+=(a[i])*z^(7-i)   # inverse order

    # initial
    inv_apoly = 1
    
    # ((((inv_poly^2)^2)^2)^2) = inv_poly^16
    power = 1
    for i in range(1,BITS):# 4 times lift a to the power of 2 (a^2,a^4,a^8,a^16)
        # lift a to the power of 2
        apoly = (apoly * apoly) % rp
        # multiply a to the accumulated product
        inv_apoly = (inv_apoly * apoly) % rp
    #inv_apoly = apoly^14 % rp

    inv_coeff = inv_apoly.list()

    degree = inv_apoly.degree()     
    assert degree < BITS

    for icoeff in range(degree+1,BITS):
        inv_coeff.append(0)

    #return result
    return inv_coeff

# SubBytes4() step 2 - affine transformation
m4 = matrix([[P(1),P(1),P(1),P(0)],
             [P(0),P(1),P(1),P(1)],
             [P(1),P(0),P(1),P(1)],
             [P(1),P(1),P(0),P(1)]])
v4 = vector ([P(0),P(1),P(1),P(0)]) # 0x6 = 0110 
#v =   vector([1,1,0,0,0,1,1,0]) # 0x63 = b0110'0011 
inv_m4 = m4.inverse()
inv_v4 = m4.inverse() * v4

# SubBytes4() step 1 - inverse transformation; calculate inverses mod 16
def invt4(in_byte):
     out_byte = []
     out_byte = inverse4(in_byte)
     #out_byte = precalculated_inverse4(in_byte) # !
     return out_byte

def afft4(in_byte,m,v):
    out_byte = []
    vector_byte = vector(in_byte)
    #print "lin",m4*vector_byte
    aff_vector_byte = m*vector_byte + v
    aff_byte = aff_vector_byte.list()
    #aff_byte.reverse()          # WARNING!
    for ibit in range(0,BITS):
        out_byte.append(aff_byte[ibit])
    return out_byte

# SubByte4 on one word
def sbox4(a):
    b = invt4(a)
    c = afft4(b,m4,v4)
    return c

# add two BITS-bit strings in GF(2) ie. xor their elements
def add(a,b):

    c = [None]*BITS
    # constant as polygon
    for i in range(0,BITS):
        c[i] = a[i]+b[i]
    return c

# precomputed round constant for the key schedule (computed with rcon4() from lexs-compute-2x2x4.sage)
#rc = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

# generates the round constants for the key schedule
def rcon4():
    rc = []
    # initial value for ro is 0x01 = 0*z^1 + 1*z^0 = 1
    #ro = 1
    x = 1
    rc.append([1,0,0,0])        # 0x1 = [1,0,0,0] = a, so that a[0]=1,a[1]=0,a[2]=0,a[3]=0 
    #rc.append(x)
    for i in range(1,10):       # 10 keys
        rc.append([])
        x = z*x % rp           # calc constant
        #print "x[",i,"]",x
        xcf = x.list()          # get pol coeffs
        d = x.degree()          # get pol degree
        assert d < BITS
        for j in range(d+1,BITS): 
            xcf.append(0)       # padd the zero coeffs with 0-s
        for j in range(0,BITS):
            rc[i].append(xcf[j]) # append constant
    return rc

# key schedule; K is the initial key
def ks(key):
    ek = []                     # extended key

    # generate round constants
    rc = rcon4()

    # copy the initial key as the first element of the extended key
    ek.append([])
    for i in range(0,NSTATE):
        ek[0].append(key[i])

    print "ks: ek",ek

    #for i in range(1,10): # 9 keys
    for i in range(1,2): # 2 keys

        s0 = invt4(ek[i-1][(NSTATE-BITS):NSTATE]) # k_{i-1},1
        t0 = afft4(s0,m4,v4)                      # L(s0) + d
     
        print "sb-output-0[",i,"]",t0
        
        s1 = invt4(ek[i-1][(NSTATE-2*BITS):NSTATE-BITS]) # k_{i-1},0
        t1 = afft4(s1,m4,v4)                             # L(s1) + d

        print "sb-output-1[",i,"]",t1
    
        ek.append([])
        # [ki00, ki10  |ki01  ,ki11  ] | [k_{i+1}00,k_{i+1}10|k_{i+1}01,k_{i+1}11]
        # [BITS, 2*BITS|3*BITS,4*BITS]
    
        # calculate first element: ki00,ki10
        t0 = add(t0,rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])      # L(s0) + d + rc_i + k_{i-1},00
        ki10 = add(t1,ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
                    
        # calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                
    
        ek[i] = ki00 + ki10 + ki01 + ki11
        #ek[i] = ki00 + ki10 + add(ki01,ki00) + add(ki11,ki10) #vincent correction, 20090117

    return ek
# key schedule }

# store the transformed equations in file
def write_eqs_to_file(e):

    n = 32

    s = str(e)
    # replace xi with x[i] and yi with y[i]
    for i in range(0,n+1):
        s1 = 'x'+str(n-i)
        s2 = 'x['+str(n-i)+']'
        s = s.replace(s1,s2)
    # replace x[8:16] with y[0:8]
    #for i in range(N/2,N):
    #    s1 = 'x['+str(i)+']'
    #    s2 = 'y['+str(i-N/2)+']'
    #    s = s.replace(s1,s2)

    fname = 'ks-eqs-2x2x4.txt'  
    print "fname",fname

    f = open(fname,"w")
    f.write("e = ")
    f.write(s)
    f.close()
    return s


#k = key
print "k",k

#k0=ekey[0]
#k0=k[0]
#ek=ks(k0)
ek=ks(k[0])

# write key equations for two keys
eq = []
for i in range(0,NSTATE):
    eq.append(k[1][i] + ek[1][i])

#vincent transformation, 20090117
for i in range(NSTATE/2,NSTATE):
    eq[i] = eq[i] + eq[i-NSTATE/2]

#s = write_eqs_to_file(eq)
#print "\ns",s
#load "ks-2x2x2x4.sage"
load "keysched-eqs-precomputed.sage"

t = keqs(k[0],k[1])
