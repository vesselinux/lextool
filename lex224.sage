#
# Small-scale variants of LEX built on small-scale variants of AES based on Murphy-Robshaw paper "Small Scale Variants of the AES"
#
# The code here is based on code from SYMAES and more specifically the file fwd.sage - the implementation of AES128 in encryption (forward) mode. SYMAES is available at: http://www.ecrypt.eu.org/tools/symaes
#
# All small-scale variants of LEX work in GF(2^4) (the full-scale of LEX works in GF(2^8))
#
# 
# Original name: lexs.sage
# 

# NOTE: bits are ordered like this: 0x1 = [1,0,0,0] = a, so that a[0]=1,a[1]=0,a[2]=0,a[3]=0 

# number of bits in one word
BITS = 4
# number of rows and columns in the state
ROWS = 2
# number of rows and columns in the state
COLS = 2

# Number of rounds for which we generate equations
ROUNDS=4

# Number of leaked bits after every round
#LEAK = ROWS*COLS*BITS - 4
LEAK = 4

# Number of leaked bits after round zero; normally should be equal to LEAK
# Note: Equations are generate starting from round 1 so that we can use
# the bits from the leak after round zero (the initial leak)
#INITIAL_LEAK = ROWS*COLS*BITS - 4 
INITIAL_LEAK = LEAK

# 1 aes round is composed of one application of (SubBytes+ShiftRows+MixColumns+AddRoundKey)
AES_ROUNDS = 9
# 1 lex round is composed of 10 aes rounds
# LEX_ROUNDS = 1 + divmod(ROUNDS,AES_ROUNDS)[0]
# keys for Lex: 0,1,2,...,9
KEYS = 2
# list of equations
E = []
# number of equations
En = 0

# Number of plaintext variables x
Np = ROWS*COLS*BITS

# the flag is 1 when we DO NOT introduce ciphertext variables y
FLAG_ONLY_LEAKS=0

# the flag is 1 when we want to compute the round keys from the extended key
FLAG_CALCULATE_ROUND_KEYS=0

# Number of key variables
if FLAG_CALCULATE_ROUND_KEYS==0:
    if ROUNDS < AES_ROUNDS:
        Nk = Np*ROUNDS
    else:  
        Nk = Np*AES_ROUNDS# number of key variables : k
elif FLAG_CALCULATE_ROUND_KEYS==1:
    Nk = Np

# Number of ciphertext variables
if FLAG_ONLY_LEAKS==0:
    Nc = Np*ROUNDS# number of ciphertext variables : y
elif FLAG_ONLY_LEAKS==1:
    Nc=0

# Total number of variables
N = Np + Nk + Nc

# Define the ring of Boolean polynomials
#P = PolynomialRing(GF(2), N, 'x',order='lex')
P = BooleanPolynomialRing(N,'x',order='lex') # WARNING!!!!
#P.<x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,x91,x92,x93,x94,x95,x96,x97,x98,x99,x100,x101,x102,x103,x104,x105,x106,x107,x108,x109,x110,x111,x112,x113,x114,x115,x116,x117,x118,x119,x120,x121,x122,x123,x124,x125,x126,x127,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96,k97,k98,k99,k100,k101,k102,k103,k104,k105,k106,k107,k108,k109,k110,k111,k112,k113,k114,k115,k116,k117,k118,k119,k120,k121,k122,k123,k124,k125,k126,k127,x128,x129,x130,x131,x132,x133,x134,x135,x136,x137,x138,x139,x140,x141,x142,x143> = BooleanPolynomialRing(N, order='lex')
# Generate field equations: xi^2-xi=0
FE=[]
for i in range(0,N):
    FE.append(P.gen(i)^2+P.gen(i))
#print "FE",FE

# Generate plaintext variables
x = []
for i in range(0,Np):
    x.append(P.gen(i))

# Generate keys variables
k = []
if FLAG_CALCULATE_ROUND_KEYS==0:
    if ROUNDS < AES_ROUNDS:
        R = ROUNDS
    else:  
        R = AES_ROUNDS
    for r in range(0,R):
        k.append([])
    #    for i in range(ROWS*COLS*BITS*(r+1),ROWS*COLS*BITS*(r+1)+ROWS*COLS*BITS):
        for i in range((Np+Np*r)+0,(Np+Np*r)+Np):
            k[r].append(P.gen(i))

elif FLAG_CALCULATE_ROUND_KEYS==1:
    for i in range(Np+0,Np+Np):
        k.append(P.gen(i))

# Generate ciphertext variables
if FLAG_ONLY_LEAKS==0:
    # ciphertext(s)
    y = []
    #for r in range(0,AES_ROUNDS*LEX_ROUNDS):
    for r in range(0,ROUNDS):
        y.append([])
        for i in range((Np+Nk+Np*r)+0,(Np+Nk+Np*r)+Np):
            y[r].append(P.gen(i))    

R.<z> = P.fraction_field()[]
if BITS==4:
    rp = z^4+z+1#rijndael polynomial for GF16
if BITS==8:
    rp = z^8+z^4+z^3+z+1

 # Calculate the inverse in GF(2^4) modulo 
# the Rijndael polynomial for SR(10, 2,2,4)

# a /in GF16 = GF(2^4) 
def inverse4(a):

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

# 
# Precalculated equations from inverse4(a)
# 
def precalculated_inverse4(a):

    b = [None]*BITS

    if a[0] == 0: a0 = P(0)
    elif a[0] == 1: a0 = P(1)
    else: a0 = a[0]            

    if a[1] == 0: a1 = P(0)
    elif a[1] == 1: a1 = P(1)
    else: a1 = a[1]            

    if a[2] == 0: a2 = P(0)
    elif a[2] == 1: a2 = P(1)
    else: a2 = a[2]            

    if a[3] == 0: a3 = P(0)
    elif a[3] == 1: a3 = P(1)
    else: a3 = a[3]            

    b[0]  = a0*a1*a2 + a0*a2 + a0 + a1*a2*a3 + a1*a2 + a1 + a2 + a3
    b[1]  = a0*a1*a3 + a0*a1 + a0*a2 + a1*a2 + a1*a3 + a3
    b[2]  = a0*a1 + a0*a2*a3 + a0*a2 + a0*a3 + a2 + a3
    b[3]  = a0*a3 + a1*a2*a3 + a1*a3 + a1 + a2*a3 + a2 + a3

    return b

# 
# SubBytes4() step 1 - inverse transformation; calculate inverses mod 16
# Note: a byte is 4 bits in GF(2^4)
# 
def invt4(in_byte):
     out_byte = []
     out_byte = inverse4(in_byte)
     #out_byte = precalculated_inverse4(in_byte) # !
     return out_byte

# Matrix for SubBytes4() step 2 - affine transformation
m4 = matrix([[P(1),P(1),P(1),P(0)],
             [P(0),P(1),P(1),P(1)],
             [P(1),P(0),P(1),P(1)],
             [P(1),P(1),P(0),P(1)]])

# Vector for SubBytes4() step 2 - affine transformation
v4 = vector ([P(0),P(1),P(1),P(0)]) # 0x6 = 0110 

# 
# SubBytes4() step 2 - affine transformation
# 
def afft4(in_byte):
    out_byte = []
    vector_byte = vector(in_byte)
    aff_vector_byte = m4*vector_byte + v4
    aff_byte = aff_vector_byte.list()
    for ibit in range(0,BITS):
        out_byte.append(aff_byte[ibit])
    return out_byte

# 
# SubBytes4 : step 1 + step 2
# 4-bit Rijndael SubBytes transformation: invt4() + afft4()
# 
def sb4(in_state):
    out_state = []
    # cycle through all bytes of the state
    for ibyte in range(0,ROWS*COLS):
        inv_byte = invt4(in_state[BITS*ibyte:(BITS*ibyte+BITS)])
        out_byte = afft4(inv_byte)
        for ibit in range(0,BITS):
            out_state.append(out_byte[ibit])
    return out_state

# 
# ShiftRows4
# left shift constants for each row
# 
C0 = 0; C1 = 1; C2 = 2; C3 = 3
def sr4(in_state):
    out_state = []
    C = [C0,C1,C2,C3]
    for irow in range(0,ROWS):
        for ibyte in range(0,COLS):
            # index of shifted byte
            ishift = irow*COLS + ((ibyte+C[irow]) % COLS) # COLS == Nb
            # copy the bits of the shifted byte
            for ibit in range(0,BITS):
                out_state.append(in_state[ishift*BITS + ibit])
                #print "appended bit ",byte*8 + ibit
    return out_state    

# 
# MixColumn4
# 
def mc4(in_state):

    #print "mc4: in_state",in_state

    # state is 1x1
    mz1 = matrix ([1])
    # state is 2x2
    mz2 = matrix ([[z+1,z],
                   [z,z+1]])
    # state is 4x4 - full-scale
    mz4 = matrix ([[z,z+1,1,1],
                   [1,z,z+1,1],
                   [1,1,z,z+1],
                   [z+1,1,1,z]])

    # Debug:
    # for i in range(0,len(in_state)):
    #     print "\nin_state[",i,"]",in_state[i]

    # Output state
    out_state = []
    # List of polynomials in z
    zpolys = []

    # 
    # represent each byte of the state as a polynomial in z
    # for SR(2,2,4) the order should be: 
    # 
    # zpoly[0] = in_state[3]*z^3  + in_state[2]*z^2  + in_state[1]*z^1  + in_state[0]*z^0
    # zpoly[1] = in_state[7]*z^3  + in_state[6]*z^2  + in_state[5]*z^1  + in_state[4]*z^0
    # zpoly[2] = in_state[11]*z^3 + in_state[10]*z^2 + in_state[9]*z^1  + in_state[8]*z^0
    # zpoly[3] = in_state[15]*z^3 + in_state[14]*z^2 + in_state[13]*z^1 + in_state[12]*z^0
    # 
    for ibyte in range(0,ROWS*COLS):
        ipoly = 0
        for ibit in range(0,BITS):
            #print "ibyte*BITS+ibit",ibyte*BITS+ibit
            ipoly+= (in_state[ibyte*BITS+ibit])*z^ibit
        #print "\nipoly[",ibyte,"]",ipoly
        zpolys.append(ipoly)

    # 
    # put the polynomials of the state in 4x4 matrix; 'sz' stands for "state-polynomial in z"
    # 
    #    sz = matrix([[zpolys[0], zpolys[1], zpolys[2], zpolys[3]],
    #                 [zpolys[4], zpolys[5], zpolys[6], zpolys[7]],
    #                 [zpolys[8], zpolys[9], zpolys[10],zpolys[11]],
    #                 [zpolys[12],zpolys[13],zpolys[14],zpolys[15]]])
    # 
    # SR(r,1,1,4)
    if (ROWS == 1) and (COLS == 1): 
        mz = mz1
        sz = matrix([zpolys[0]])
    # SR(r,2,1,4)
    elif (ROWS == 2) and (COLS == 1): 
        mz = mz2
        sz = matrix([[zpolys[0]],
                     [zpolys[1]]])
    # SR(r,2,2,4)
    elif (ROWS == 2) and (COLS == 2): 
        mz = mz2
        sz = matrix([[zpolys[0], zpolys[1]],
                     [zpolys[2], zpolys[3]]])
    # SR(r,2,4,4)
    elif (ROWS == 2) and (COLS == 4): 
        mz = mz2
        sz = matrix([[zpolys[0], zpolys[1], zpolys[2], zpolys[3]],
                     [zpolys[4], zpolys[5], zpolys[6], zpolys[7]]])
    # SR(r,4,4,4)
    elif ROWS == 4: 
        mz = mz4
    else: 
        print "ERROR invalid rows",ROWS

    # multiply matrices modulo the rijndael polynomial
    mc = (mz*sz) % rp
    # put the coefficients of the polynomials of mc in the output state
    for irow in range(0,ROWS):
        for ibyte in range(0,COLS):
            #print "mc[",irow,",",ibyte,"]=", mc[irow][ibyte]
            mcpoly = mc[irow][ibyte].list() # WARNING!

            # get the degree of the polynomial; if it is less than 7, padd the missing coeffs with zeros
            degree = mc[irow][ibyte].degree()     
            assert degree < BITS
            # pad the "missing" coefficients from degree+1 up to 7 with zeros
            for icoeff in range(degree+1,BITS):
                mcpoly.append(0)

            for ibit in range(0,BITS):
                # get the polynomial representing each byte (the "mcpoly" polynomial)
                # print "mcpoly[",ibit,"]",mcpoly[ibit]
                # put the 8 coefficents of mcpoly into 8 of the bits of the state
                out_state.append(mcpoly[ibit])
    return out_state

# 
# AddRoundKey4
# 
def ark4(in_state,rk):
    out_state = []
    # add round key
    for ibit in range(0,ROWS*COLS*BITS):
        out_bit = in_state[ibit] + rk[ibit]
        out_state.append(out_bit)
    return out_state

# 4-bit SR(r,2,2,4)/Lex round transformation: SubByte4 + ShiftRows4 + MixColumn4 + AddRoundKey4
# y = round(x) : x->SB(x)->sx->SR(sx)->rx->MC(rx)->mx->ARK(mx)->y
def r4(x, rk):

    #print "r4: x",x
    #print "r4: rk",rk
    sx = sb4(x)
    #print "r4: sx",sx
    rx = sr4(sx)
    #print "r4: rx",rx
    mx = mc4(rx)
    #print "r4: mx",mx
    y = ark4(mx,rk)
    #print "r4: y",y
    return y

# 
# r rounds of SR(r,2,2,4)
# 
def aes4(x, k):
    x = ark4(x,k[0])

    # 9 rounds
    for i in range(1,9):        # WARNING! shoud be range(1,10) !!!
        print "aes round#",i
        x = sb4(x)
        x = sr4(x)
        x = mc4(x)
        x = ark4(x,k[i])
    # 1 final round 
    print "aes round#",i+1
    x = sb4(x)
    x = sr4(x)
    y = ark4(x,k[9])

    return y

# test inverse: a = d3 = b0011 = 0x3; a = [1,1,0,0] # b = a^{-1} mod 16 = d14 = 0xE = b1110 #b = inverse4 (a) print "a",a,"b",b
def test_inverse():
    for a0 in range(0,2):
        for a1 in range(0,2):
            for a2 in range(0,2):
                for a3 in range(0,2):
                    if (ROWS==1) and (COLS==1):
                        a = [a3,a2,a1,a0]
                    elif (ROWS==2) and (COLS==1):
                        a = [a3,a2,a1,a0,a3,a2,a1,a0]
                    elif (ROWS==2) and (COLS==4):
                        a = [a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0,a3,a2,a1,a0]
                    else: print "ERROR invalid rows,cols",ROWS,COLS
                    b = inverse4 (a)
                    c = precalculated_inverse4 (a)
                    xd = afft4(b)
                    b = sb4(a)
                    c = sr4(b)
                    d = mc4(c)
                    e = ark4(d,k[0])
                    f = r4(a,k[0])
                    print "a",a,"b",b,"c",c,"d",d,"e",e,"f",f


#
# generates the round constants for the key schedule
#
def rcon4():

    rc = []
    # initial value for ro is 0x01 = 0*z^1 + 1*z^0 = 1
    #ro = 1
    x = 1
    rc.append([1,0,0,0])        # 0x1 = [1,0,0,0] = a, so that a[0]=1,a[1]=0,a[2]=0,a[3]=0 
    #rc.append(x)
    for i in range(1,AES_ROUNDS):
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

# add two BITS-bit strings in GF(2) ie. xor their elements
def add(a,b):

    c = [None]*BITS
    # constant as polygon
    for i in range(0,BITS):
        c[i] = a[i]+b[i]
    return c

# key matrix layout:
# 
# [k00 k01 k02 k03  
#  k10 k11 k12 k13
#  k20 k21 k22 k23
#  k30 k31 k32 k33]
# 
# key vector layout by spec:
# 
# [[k00,k10,k20,k30],[k01,k11,k21,k31],[k02,k12,k22,k32],[k03,k13,k23,k33]]
# 
# my key vector layout (because the plaintext has the same layout):
# 
# [[k00 k01 k02 k03],[k10 k11 k12 k13],[k20 k21 k22 k23],[k30 k31 k32 k33]]
# 
# 
# key schedule; K is the initial key
# 
# NOTE: For more details see paper: Matt Robshaw, "Small-scale Variants of the AES"
# 
def ks(key):
    #k = []                      # key in words by 4
    ek = []                     # extended key

    # generate round constants
    rc = rcon4()

    # copy the initial key as the first element of the extended key
    ek.append([])
    for i in range(0,ROWS*COLS*BITS):
        ek[0].append(key[i])

    if ROUNDS < AES_ROUNDS:
        R = ROUNDS
    else:
        R = AES_ROUNDS

    if (ROWS==1):
        
        for i in range(1,R): # 9 keys

            # ROWS*COLS - the number of 4-bit elements in one round key
            #prev_key = ek[i-1]
            #print "ek[",i-1,"]",ek[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS] # the last 4 bits of the previous key
            #t = ek[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS]
            s0 = invt4(ek[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS])
            #print "s0", s0

            # SR(R,1,1,4)
            if (COLS==1):
                t = afft4(s0)   # L(s0) + d
                t = add(t,rc[i]) # L(s0) + d + rc_i
                ek.append(t)
                
            # SR(R,1,2,4)
            if (COLS==2):

                ek.append([])

                t = afft4(s0)   # L(s0) + d
                t = add(t,rc[i]) # L(s0) + d + rc_i
                ki0 = add(t,ek[i-1][0:BITS]) # L(s0) + d + rc_i + k_{i-1},0
                #print "ki0[",i,"]",ki0
                ki1 = add(ki0, ek[i-1][BITS:2*BITS]) # L(s0) + d + rc_i + k_{i-1},0 + k_{i-1},1
                #print "ki1[",i,"]",ki1
                for j in range(0,BITS):
                    ek[i].append(ki0[j])
                for j in range(0,BITS):
                    ek[i].append(ki1[j])
                #print "ek[",i,"]",ek[i]

    elif (ROWS==2):
        
        for i in range(1,R): # 9 keys

            s0 = invt4(ek[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS]) # k_{i-1},1
            #print "s0[",i,"]",s0
            s1 = invt4(ek[i-1][(ROWS*COLS*BITS-2*BITS):ROWS*COLS*BITS-BITS]) # k_{i-1},0
            #print "s1[",i,"]",s1

            # SR(R,2,1,4)
            if (COLS==1):

                ek.append([])

                t = afft4(s0)   # L(s0) + d
                ki0 = add(t,rc[i]) # L(s0) + d + rc_i
                ki1 = afft4(s1)   # L(s1) + d
                for j in range(0,BITS):
                    ek[i].append(ki0[j])
                for j in range(0,BITS):
                    ek[i].append(ki1[j])

            # SR(R,2,2,4)
            elif (COLS==2):

                ek.append([])
                # [ki00, ki10  |ki01  ,ki11  ] | [k_{i+1}00,k_{i+1}10|k_{i+1}01,k_{i+1}11]
                # [BITS, 2*BITS|3*BITS,4*BITS]

                # calculate first element: ki00,ki10
                t = afft4(s0)   # L(s0) + d
                t = add(t,rc[i]) # L(s0) + d + rc_i
                ki00 = add(t,ek[i-1][0:BITS]) # L(s0) + d + rc_i + k_{i-1},00
                t = afft4(s1)   # L(s1) + d
                ki10 = add(t,ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
                
                # calculate second element: ki01,ki11
                ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
                ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                

                for j in range(0,BITS):
                    ek[i].append(ki00[j])
                for j in range(0,BITS):
                    ek[i].append(ki10[j])
                for j in range(0,BITS):
                    ek[i].append(ki01[j])
                for j in range(0,BITS):
                    ek[i].append(ki11[j])

                #print "k00[",i,"]",ki00
                #print "k10[",i,"]",ki10
                #print "k01[",i,"]",ki01
                #print "k11[",i,"]",ki11

            elif (COLS==4):

                ek.append([])
                # [ki00, ki10  |ki01  ,ki11   |ki02  ,ki12 |ki03  ,ki13] | [k_{i+1}00,k_{i+1}10|k_{i+1}01,k_{i+1}11|k_{i+1}02,k_{i+1}12|k_{i+1}03,k_{i+1}13]
                # [BITS, 2*BITS|3*BITS,4*BITS]

                # calculate first element: ki00,ki10
                t = afft4(s0)   # L(s0) + d
                t = add(t,rc[i]) # L(s0) + d + rc_i
                ki00 = add(t,ek[i-1][0:BITS]) # L(s0) + d + rc_i + k_{i-1},00
                t = afft4(s1)   # L(s1) + d
                ki10 = add(t,ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10                
                # calculate second element: ki01,ki11
                ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
                ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                
                # calculate the third element: ki02,ki12
                ki02 = add(ki01,ek[i-1][4*BITS:5*BITS])
                ki12 = add(ki11,ek[i-1][5*BITS:6*BITS])
                # calculate the third element: ki03,ki13
                ki03 = add(ki02,ek[i-1][6*BITS:7*BITS])
                ki13 = add(ki12,ek[i-1][7*BITS:8*BITS])

                for j in range(0,BITS):
                    ek[i].append(ki00[j])
                for j in range(0,BITS):
                    ek[i].append(ki10[j])
                for j in range(0,BITS):
                    ek[i].append(ki01[j])
                for j in range(0,BITS):
                    ek[i].append(ki11[j])
                for j in range(0,BITS):
                    ek[i].append(ki02[j])
                for j in range(0,BITS):
                    ek[i].append(ki12[j])
                for j in range(0,BITS):
                    ek[i].append(ki03[j])
                for j in range(0,BITS):
                    ek[i].append(ki13[j])

    elif (ROWS==4):
        
        for i in range(1,R): # 9 keys

            #print "ek[",i-1,"]",ek[i-1]
            #print "ROWS*COLS*BITS-BITS",ROWS*COLS*BITS-BITS,"ROWS*COLS*BITS",ROWS*COLS*BITS

            s0 = invt4(ek[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS]) # k_{i-1},3
            s1 = invt4(ek[i-1][(ROWS*COLS*BITS-2*BITS):ROWS*COLS*BITS-BITS]) # k_{i-1},2
            s2 = invt4(ek[i-1][(ROWS*COLS*BITS-3*BITS):ROWS*COLS*BITS-2*BITS]) # k_{i-1},1
            s3 = invt4(ek[i-1][(ROWS*COLS*BITS-4*BITS):ROWS*COLS*BITS-3*BITS]) # k_{i-1},0

            #print "s0[",i,"]",s0
            #print "s1[",i,"]",s1
            #print "s2[",i,"]",s2
            #print "s3[",i,"]",s3

            # SR(R,4,1,4)
            if (COLS==1):

                ek.append([])

                # ki00,ki10,ki20,ki30  
                t = afft4(s0) 
                ki00 = add(t,rc[i]) # L(s0) + d + rc_i
                ki10 = afft4(s1) 
                ki20 = afft4(s2) 
                ki30 = afft4(s3) 

                for j in range(0,BITS):
                    ek[i].append(ki00[j])
                for j in range(0,BITS):
                    ek[i].append(ki10[j])
                for j in range(0,BITS):
                    ek[i].append(ki20[j])
                for j in range(0,BITS):
                    ek[i].append(ki30[j])

            # SR(R,4,2,4)
            elif (COLS==2):

                # [ki00,ki10,ki20,ki30][ki01,ki11,ki21,ki31]
                ek.append([])

                # calculate first element: ki00,ki10,ki20,ki30
                t = afft4(s0) 
                t = add(t,rc[i]) # L(s0) + d + rc_i
                ki00 = add(t,ek[i-1][0:BITS])

                t= afft4(s1) 
                ki10 = add(t,ek[i-1][BITS:2*BITS]) 

                t = afft4(s2) 
                ki20 = add(t,ek[i-1][2*BITS:3*BITS]) 

                t = afft4(s3) 
                ki30 = add(t,ek[i-1][3*BITS:4*BITS]) 

                # calculate second element: ki01,ki11,ki21,ki31
                ki01 = add(ki00,ek[i-1][4*BITS:5*BITS])
                ki11 = add(ki10,ek[i-1][5*BITS:6*BITS])
                ki21 = add(ki20,ek[i-1][6*BITS:7*BITS])
                ki31 = add(ki30,ek[i-1][7*BITS:8*BITS])

                for j in range(0,BITS):
                    ek[i].append(ki00[j])
                for j in range(0,BITS):
                    ek[i].append(ki10[j])
                for j in range(0,BITS):
                    ek[i].append(ki20[j])
                for j in range(0,BITS):
                    ek[i].append(ki30[j])
                for j in range(0,BITS):
                    ek[i].append(ki01[j])
                for j in range(0,BITS):
                    ek[i].append(ki11[j])
                for j in range(0,BITS):
                    ek[i].append(ki21[j])
                for j in range(0,BITS):
                    ek[i].append(ki31[j])

    return ek

# key schedule for one round
def ks_round(k1):
    #k = []                      # k1 in words by 4
    k2 = []                     # extended k1

    # generate round constants
    rc = rcon4()

    # copy the initial k1 as the first element of the extended k1
    k2.append([])
    for i in range(0,ROWS*COLS*BITS):
        k2[0].append(k1[i])

    R = 2
  
    for i in range(1,R): # 9 k1s

        s0 = invt4(k2[i-1][(ROWS*COLS*BITS-BITS):ROWS*COLS*BITS]) # k_{i-1},1
        #print "s0[",i,"]",s0
        s1 = invt4(k2[i-1][(ROWS*COLS*BITS-2*BITS):ROWS*COLS*BITS-BITS]) # k_{i-1},0
        #print "s1[",i,"]",s1

        k2.append([])
        # [ki00, ki10  |ki01  ,ki11  ] | [k_{i+1}00,k_{i+1}10|k_{i+1}01,k_{i+1}11]
        # [BITS, 2*BITS|3*BITS,4*BITS]

        # calculate first element: ki00,ki10
        t = afft4(s0)   # L(s0) + d
        t = add(t,rc[i]) # L(s0) + d + rc_i
        ki00 = add(t,k2[i-1][0:BITS]) # L(s0) + d + rc_i + k_{i-1},00
        t = afft4(s1)   # L(s1) + d
        ki10 = add(t,k2[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
                
        # calculate second element: ki01,ki11
        ki01 = add(ki00,k2[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,k2[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                

        for j in range(0,BITS):
            k2[i].append(ki00[j])
        for j in range(0,BITS):
            k2[i].append(ki10[j])
        for j in range(0,BITS):
            k2[i].append(ki01[j])
        for j in range(0,BITS):
            k2[i].append(ki11[j])

        #print "k00[",i,"]",ki00
        #print "k10[",i,"]",ki10
        #print "k01[",i,"]",ki01
        #print "k11[",i,"]",ki11

    return k2


# calculate round keys for 'rn' rounds starting from initial key k0
def generate_rkeys(k0,rn):
    ke = ks(k0) #xvpv
    ##print ke

    ke = ke[0:rn]  
  
    s = str(ke)
    s = s.replace('x16','k[0]')
    s = s.replace('x17','k[1]')
    s = s.replace('x18','k[2]')
    s = s.replace('x19','k[3]')
    s = s.replace('x20','k[4]')
    s = s.replace('x21','k[5]')
    s = s.replace('x22','k[6]')
    s = s.replace('x23','k[7]')
    s = s.replace('x24','k[8]')
    s = s.replace('x25','k[9]')
    s = s.replace('x26','k[10]')
    s = s.replace('x27','k[11]')
    s = s.replace('x28','k[12]')
    s = s.replace('x29','k[13]')
    s = s.replace('x30','k[14]')
    s = s.replace('x31','k[15]')

    fname = 'ks_' + str(rn) + 'x2x2x4.sage'  
    print "fname",fname

    f = open(fname,"w")
    f.write("ke = ")
    f.write(s)
    f.close()

    return ke

# print explicit equations for one round of LEX-2x2x4 in a file
def generate_eqs(x,kr):

    xo = r4(x,kr)

    s = str(xo)
    s = s.replace('x16','k[0]')
    s = s.replace('x17','k[1]')
    s = s.replace('x18','k[2]')
    s = s.replace('x19','k[3]')
    s = s.replace('x20','k[4]')
    s = s.replace('x21','k[5]')
    s = s.replace('x22','k[6]')
    s = s.replace('x23','k[7]')
    s = s.replace('x24','k[8]')
    s = s.replace('x25','k[9]')
    s = s.replace('x26','k[10]')
    s = s.replace('x27','k[11]')
    s = s.replace('x28','k[12]')
    s = s.replace('x29','k[13]')
    s = s.replace('x30','k[14]')
    s = s.replace('x31','k[15]')
    s = s.replace('x','x[')
    s = s.replace('*',']*')
    s = s.replace(' + ','] + ')
    #s = s.replace(' , ','],')
    s = s.replace(']]',']')
    s = s.replace('k[15]','k[15]]')

    fname = 'lex224-round-eqs.sage'  
    print "fname",fname

    f = open(fname,"w")
    f.write("y = ")
    f.write(s)
    f.close()

# calculate k1 from k0 according to the key schedule
def next_key(k0):
    # this produces all 10 round keys computed from k0, but we only use the 1-st
    #ke = ks(k0)
    #k1 = ke[1]  
    #k0 = [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10], k[11], k[12], k[13], k[14], k[15]]
    k1 = [k0[0] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13]*k0[15] + k0[12]*k0[14]*k0[15] + k0[12]*k0[14] + k0[12]*k0[15] + k0[12] + k0[13]*k0[14]*k0[15] + k0[13]*k0[15] + k0[13] + k0[15], k0[1] + k0[12]*k0[13]*k0[15] + k0[12]*k0[14]*k0[15] + k0[13]*k0[14]*k0[15] + k0[13]*k0[14] + k0[13] + k0[14]*k0[15] + k0[15], k0[2] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13] + k0[12]*k0[14]*k0[15] + k0[12] + k0[13]*k0[14] + k0[13]*k0[15] + k0[14]*k0[15] + k0[14] + k0[15] + 1, k0[3] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13]*k0[15] + k0[12]*k0[13] + k0[12]*k0[15] + k0[12] + k0[14]*k0[15] + k0[15], k0[4] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9]*k0[11] + k0[8]*k0[10]*k0[11] + k0[8]*k0[10] + k0[8]*k0[11] + k0[8] + k0[9]*k0[10]*k0[11] + k0[9]*k0[11] + k0[9] + k0[11], k0[5] + k0[8]*k0[9]*k0[11] + k0[8]*k0[10]*k0[11] + k0[9]*k0[10]*k0[11] + k0[9]*k0[10] + k0[9] + k0[10]*k0[11] + k0[11] + 1, k0[6] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9] + k0[8]*k0[10]*k0[11] + k0[8] + k0[9]*k0[10] + k0[9]*k0[11] + k0[10]*k0[11] + k0[10] + k0[11] + 1, k0[7] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9]*k0[11] + k0[8]*k0[9] + k0[8]*k0[11] + k0[8] + k0[10]*k0[11] + k0[11], k0[0] + k0[8] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13]*k0[15] + k0[12]*k0[14]*k0[15] + k0[12]*k0[14] + k0[12]*k0[15] + k0[12] + k0[13]*k0[14]*k0[15] + k0[13]*k0[15] + k0[13] + k0[15], k0[1] + k0[9] + k0[12]*k0[13]*k0[15] + k0[12]*k0[14]*k0[15] + k0[13]*k0[14]*k0[15] + k0[13]*k0[14] + k0[13] + k0[14]*k0[15] + k0[15], k0[2] + k0[10] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13] + k0[12]*k0[14]*k0[15] + k0[12] + k0[13]*k0[14] + k0[13]*k0[15] + k0[14]*k0[15] + k0[14] + k0[15] + 1, k0[3] + k0[11] + k0[12]*k0[13]*k0[14] + k0[12]*k0[13]*k0[15] + k0[12]*k0[13] + k0[12]*k0[15] + k0[12] + k0[14]*k0[15] + k0[15], k0[4] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9]*k0[11] + k0[8]*k0[10]*k0[11] + k0[8]*k0[10] + k0[8]*k0[11] + k0[8] + k0[9]*k0[10]*k0[11] + k0[9]*k0[11] + k0[9] + k0[11] + k0[12], k0[5] + k0[8]*k0[9]*k0[11] + k0[8]*k0[10]*k0[11] + k0[9]*k0[10]*k0[11] + k0[9]*k0[10] + k0[9] + k0[10]*k0[11] + k0[11] + k0[13] + 1, k0[6] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9] + k0[8]*k0[10]*k0[11] + k0[8] + k0[9]*k0[10] + k0[9]*k0[11] + k0[10]*k0[11] + k0[10] + k0[11] + k0[14] + 1, k0[7] + k0[8]*k0[9]*k0[10] + k0[8]*k0[9]*k0[11] + k0[8]*k0[9] + k0[8]*k0[11] + k0[8] + k0[10]*k0[11] + k0[11] + k0[15]]
    return k1

# produce real leaks for Lex(2,2,4)
# input: x - plaintext, k - key, r - rounds, l - number of bits to leak every round
def lex_leaks(x,k,ln):

    assert BITS==4

    l = []                      # list of leaks

    #print "lex4_leaks: x",x

    # calculate extended key
    ke = ks(k)
    #print "lex4_leaks: ke",ke

    for ri in range (0,ROUNDS): 

        #print "\n round#",ri, "key", ri%KEYS

        l.append([])

        # the round transformation
        x = r4(x,ke[ri%KEYS])
        # print "lex4_leaks x",x

        for i in range(0,ln):
            l[ri].append(x[i])

    return l

# generate actual leaks from a given plaintext, key and leak length
def generate_leaks(pt,key,LEAK):

    #pt  = [1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1] # plaintext (ie. x,IV,etc)
    #pt  = [P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1)] # plaintext (ie. x,IV,etc)
    #key = [1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0] # initial key
    #key = [P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0)] # initial key

    l = lex_leaks(pt,key,LEAK)
    s = str(l)
    f = open("leaks.sage","w")
    f.write("l = ")
    f.write(s)
    f.close()


# START MAIN

assert BITS == 4
# Note: we analyze Lex(2,2,4) with two round keys only
assert KEYS == 2

print "\nLEX: ROUNDS",ROUNDS,"ROWS",ROWS,"COLS",COLS,"BITS",BITS,"KEYS",KEYS
print "\nLEAK:",LEAK,"bits (out of",ROWS*COLS*BITS,") initial leak", INITIAL_LEAK 

# Generate random plaintext
pt  = [P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1))] 

# Generate random key
key = [P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1))] 

print " plaintext", pt
print " key      ", key

# Generate round keys
rkeys = []
rkeys = ks(key)
for i in range(0,2):
    print " rkey[",i,"]",rkeys[i]

# Generate the output after every round
# From it we shall extract the leaks
l = lex_leaks(pt,key,16)
#print "l", l
for i in range(0,len(l)):
    print "rleak[",i,"]",l[i]


# { --- Debug ---
# Assign the values of the first two keys
#k[0] = rkeys[0]
#k[1] = rkeys[1]
#print "k[0]", k[0]
#print "k[1]", k[1]
#a = next_key(k[0])
#b = ks_round(k[0])
#assert a == k[1]
#assert b[0] == k[0]
#assert b[1] == k[1]
# --- Debug --- }

