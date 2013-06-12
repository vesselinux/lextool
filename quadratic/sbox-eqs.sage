# 
# Original name: sbox-compute-2x2x4.sage
# 

# computing equations from the 4-bit input and 4-bit output of a single sbox
# used in lex224-quadratic-main.sage (old name: lexs-compute-2x2x4.sage)
# 
BITS=4
N=2*BITS
P=BooleanPolynomialRing(N,'x',order='lex')
# sbox input
a=[]
for i in range(0,N/2):
    a.append(P.gen(i))
# sbox output
b=[]
for i in range(N/2,N):
    b.append(P.gen(i))

# field of fractions
R.<z> = P.fraction_field()[]
#rijndael polynomial for GF16
rp = z^4+z+1

#load "leaks.sage"

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

# wirte a byte as a polynomial in z
def poly(a):
    apoly = 0
    for i in range(0,BITS):
        apoly+=(a[i])*z^i
    return apoly

def inv4_eqs(a,b):
    ap = poly(a)
    bp = poly(b)
    r = (ap * bp) % rp
    return r    


def sbox4_eqs(a,b):
    # equations list
    e = []
    #print "a",a
    #print "b",b
    # invert the linear transformation of the output
    c = afft4(b,inv_m4,inv_v4)
    #print "c",c
    apoly = poly(a)
    #print "apoly",apoly
    cpoly = poly(c)
    #print "cpoly",cpoly

    ### equation set (1): f1(z)*f2(z) = 1 (mod rp(z)) : 3 equations
    # product mod rijndael polynomial
    p = (apoly*cpoly + 1)%rp 
    #print "p0",p
    cf = p.coeffs()
    cf = p.list()
    #print "cf",cf
    d = p.degree()     
    assert d < BITS
    for i in range(d+1,BITS):
        cf.append(0)
    #print "cf",cf
    #if a != [0,0,0,0]:          # if a=0 then a*b != 1, so we don't have equations for this case!
    #e=e+cf[1:BITS]               # # don't add the first equation 
    e=e+cf               # add also the first equation (!)

    ### equation set (2): {f1(z)}^2*f2(z) = f1(z) (mod rp(z)) : 4 equations
    p = (apoly*apoly*cpoly+apoly)%rp 
    #print "p1",p
    cf = p.coeffs()
    cf = p.list()
    #print "cf",cf
    d = p.degree()     
    assert d < BITS
    for i in range(d+1,BITS):
        cf.append(0)
    #print "cf",cf
    e=e+cf

    ### equation set (3): f1(z)*{f2(z)}^2 = f2(z) (mod rp(z)) : 4 equations
    p = (apoly*cpoly*cpoly+cpoly)%rp 
    #print "p2",p
    cf = p.coeffs()
    cf = p.list()
    #print "cf",cf
    d = p.degree()     
    assert d < BITS
    for i in range(d+1,BITS):
        cf.append(0)
    #print "cf",cf
    e=e+cf

    return e

#
# test vectors for inversion, linear map and sbox for SR(n,1,1,4) according to Matthew-Robshaw
#
# inversion in GF(2^4)
# input:  0 1 2 3 4 5 6 7 8 9 a b c d e f
# output: 0 1 9 e d b 7 6 f 2 c 5 a 4 3 8
#
# GF(2)-linear map in GF(2^4) (afft4() without the addition of constant)
# input:  0 1 2 3 4 5 6 7 8 9 a b c d e f
# output: 0 d b 6 7 a c 1 e 3 5 8 9 4 2 f
#
# look-up table for the entire sbox in GF(2^4)
# input:  0 1 2 3 4 5 6 7 8 9 a b c d e f
# output: 6 b 5 4 2 e 7 a 9 d f c 3 1 0 8
#
def sb4_lut(b):
    if   b==[0, 0, 0, 0]: return [0, 1, 1, 0]#0
    elif b==[1, 0, 0, 0]: return [1, 1, 0, 1]# 1
    elif b==[0, 1, 0, 0]: return [1, 0, 1, 0]# 2
    elif b==[1, 1, 0, 0]: return [0, 0, 1, 0]# 3
    elif b==[0, 0, 1, 0]: return [0, 1, 0, 0]# 4
    elif b==[1, 0, 1, 0]: return [0, 1, 1, 1]# 5
    elif b==[0, 1, 1, 0]: return [1, 1, 1, 0]# 6
    elif b==[1, 1, 1, 0]: return [0, 1, 0, 1]# 7
    elif b==[0, 0, 0, 1]: return [1, 0, 0, 1]# 8
    elif b==[1, 0, 0, 1]: return [1, 0, 1, 1]# 9
    elif b==[0, 1, 0, 1]: return [1, 1, 1, 1]# 0xa
    elif b==[1, 1, 0, 1]: return [0, 0, 1, 1]# 0xb
    elif b==[0, 0, 1, 1]: return [1, 1, 0, 0]# 0xc
    elif b==[1, 0, 1, 1]: return [1, 0, 0, 0]# 0xd
    elif b==[0, 1, 1, 1]: return [0, 0, 0, 0]# 0xe
    elif b==[1, 1, 1, 1]: return [0, 0, 0, 1]# 0xf
    else: print "ERROR sb4_lut() invalid input",b

def bit_to_hex(b):
    if   b==[0, 0, 0, 0]: return 0
    elif b==[1, 0, 0, 0]: return 1
    elif b==[0, 1, 0, 0]: return 2
    elif b==[1, 1, 0, 0]: return 3
    elif b==[0, 0, 1, 0]: return 4
    elif b==[1, 0, 1, 0]: return 5
    elif b==[0, 1, 1, 0]: return 6
    elif b==[1, 1, 1, 0]: return 7
    elif b==[0, 0, 0, 1]: return 8
    elif b==[1, 0, 0, 1]: return 9
    elif b==[0, 1, 0, 1]: return 0xa
    elif b==[1, 1, 0, 1]: return 0xb
    elif b==[0, 0, 1, 1]: return 0xc
    elif b==[1, 0, 1, 1]: return 0xd
    elif b==[0, 1, 1, 1]: return 0xe
    elif b==[1, 1, 1, 1]: return 0xf
    else: print "ERROR bit_to_hex() cannot convert",b

#b = sbox4(a)
def test_sbox4():
    for b3 in range(0,2):
        for b2 in range(0,2):
            for b1 in range(0,2):
                for b0 in range(0,2):
                    b=[b0,b1,b2,b3]
                    print hex(bit_to_hex(b)),hex(bit_to_hex((sbox4(b))))
                    e = sbox4_eqs(b,sbox4(b))
                    print "e",e
                    assert e == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#test_sbox4()

def write_eqs_to_file(e):

    s = str(e)
    # replace xi with x[i] and yi with y[i]
    for i in range(0,N):
        s1 = 'x'+str(i)
        s2 = 'x['+str(i)+']'
        s = s.replace(s1,s2)
    # replace x[8:16] with y[0:8]
    for i in range(N/2,N):
        s1 = 'x['+str(i)+']'
        s2 = 'y['+str(i-N/2)+']'
        s = s.replace(s1,s2)

    fname = 'sbox-eqs-2x2x4.sage'  
    print "fname",fname

    f = open(fname,"w")
    f.write("e = ")
    f.write(s)
    f.close()
    return s

#a[0]=0
#a[1]=1
#a[2]=0
#a[3]=0

#b[0]=1
#b[1]=0
#b[2]=1
#b[3]=0

#print "a",a
#print "b",b

#e = sbox4_eqs(a,b)
#print "e",e

#s = write_eqs_to_file(e)
#print s

#k = [[[1, 1, 1, 0, 0, 0, 1, 1], [1, 0, 0, 0, 1, 1, 1, 0]], [[1, 1, 0, 1, 0, 1, 0, 1], [1, 1, 1, 1, 1, 0, 0, 0]]]
#ek=[[1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0],[0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0]]

#nle1[ 1 ] [1, 0, 0, 0] , [1, 1, 0, 1]
#nle2[ 1 ] [1, 1, 1, 0] , [0, 1, 0, 1]
#a=[1, 0, 0, 0]; b=[1, 1, 0, 1]
#c=[1, 1, 1, 0]; d=[0, 1, 0, 1]

#print "a",a
#a=ek[0][8:12]
#b=ek[1][0:4]
#print "a",a
#t=sbox4(a)
#print "t",t
#print "b",b

#pt = [P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1)] # plaintext (ie. x,IV,etc)
#p=[None]*32;p[31]=P(0);p[30]=P(0);p[29]=P(0);p[28]=P(0);p[27]=P(1);p[26]=P(1);p[25]=P(1);p[24]=P(0);p[23]=P(1);p[22]=P(0);p[21]=P(1);p[20]=P(1);p[19]=P(1);p[18]=P(0);p[17]=P(1);p[16]=P(0);p[15]=P(0);p[14]=P(1);p[13]=P(1);p[12]=P(1);p[11]=P(0);p[10]=P(0);p[9]=P(0);p[8]=P(1);p[7]=P(1);p[6]=P(1);p[5]=P(0);p[4]=P(0);p[3]=P(0);p[2]=P(1);p[1]=P(1);p[0]=P(1)

#l0 = [P(0),P(1),P(0),P(1),P(0),P(0),P(0),P(0),P(0),P(0),P(1),P(1),P(0),P(0),P(1),P(1)]
#l1 = [P(1),P(1),P(0),P(1),P(1),P(1),P(0),P(1),P(0),P(1),P(1),P(0),P(1),P(0),P(0),P(1)]

#x0=pt
#x1=l0
#x2=l1


#sbl=[]
#for i in range(0,len(l4)):
    #print "i",i
    #sbl.append([])
#    li0=sbox4(l4[i][0:4])
#    li1=sbox4(l4[i][4:8])
##    li2=sbox4(l4[i][8:12])
#    li3=sbox4(l4[i][12:16])
#    sbl.append(li0+li1+li2+li3)
#print "sbl",sbl
#assert t==ek[1][0:4]

#e = sbox4_eqs(a,b)
#print "e",e
#print "ek",ek
