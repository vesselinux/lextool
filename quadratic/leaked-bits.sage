#
# Compute the full output from ROUNDS rounds of Lex224
# for a given master key (16 bits) and plaintext (16 bits)
# From these outputs later the leaks are extracted
# Note: most the code in this file is copied from lex224.sage; small part is copied from sbox-eqs.sage
# 

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
# SubVytes4() step 2 - affine transformation
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

# produce real leaks for Lex(2,2,4)
# input: x - plaintext, k - key, r - rounds, l - number of bits to leak every round
def lex_leaks(x,k,ln):

    assert BITS==4

    l = []                      # list of leaks

    print "lex4_leaks: x",x

    # calculate extended key
    ke = ks(k)
    print "lex4_leaks: ke",ke

    for ri in range (0,AES_ROUNDS): 

        #print "\n round#",ri, "key", ri%KEYS

        l.append([])

        # the round transformation
        x = r4(x,ke[ri%KEYS])
        # print "lex4_leaks x",x

        # store the leaks as elements of the ring P
        for i in range(0,ln):
            l[ri].append(P(x[i]))

    return l

# SubByte4 on one word
#def sbox4(a):
#    b = invt4(a)
#    c = afft4(b,m4,v4)
#    return c
###

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

# 
# compute the Sbox outputs of the leaks
# (used in the construction of the Quadratic equations)
# 
def sbox_leaks(l4):
    sbl=[]
    for i in range(0,len(l4)):
        #print "i",i
        #sbl.append([])
        li0=sb4_lut(l4[i][0:4])
        li1=sb4_lut(l4[i][4:8])
        li2=sb4_lut(l4[i][8:12])
        li3=sb4_lut(l4[i][12:16])
        sbl.append(li0+li1+li2+li3)
        # print "sbl",sbl
    return sbl
