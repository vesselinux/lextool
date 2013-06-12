# Original file: ks2-compute-2x2x4.sage
#
# computing equations for the key schedule of lex 2x2x4 according to the DeCanniere-Biruykov arrangement of the key variables (see next)
# used lexs-compute-2x2x4.sage, sbox-compute-2x2x4.sage, sbox-solve-2x2x4.sage
#
# the variables of the key equations are: all bits of the initial key (16 bits) plus the inputs and the outputs of the sbox-es during the key schedule. in total we have 16 + 8*16 + 8 = 9*16 + 8 = 152 
#

#
# the key variables are organized as follow:
#
# k[0] = [initial_key[0:8],initial_key[8:16]==sbox-input[0:8]]
# k[1] = [sbox-output[0:8],sbox-input[8:16]]
# k[2] = [sbox-output[0:8],sbox-input[8:16]]
# ...
# k[8] = [sbox-output[0:8],sbox-input[8:16]]
# k[9] = [sbox-output[0:8],sbox-input[8:16]]
#
# in total we have 9*16+8=152 variables 
# every (sbox-input[0:8],sbox-output[0:8]) pair results in 2*11=22 sbox non-linear equations. we have 9 such pairs, thus: 198 equations 
# every (sbox-output[0:8],sbox-input[8:16]) pair results in 8 linear equations. we have 9 such pairs, thus: 72 linear equations
# in total we have: 198 + 72 = 270 equations in 152 variables
#
#load "sbox-eqs-2x2x4.sage"
load "sbox-eqs-precomputed.sage"

# pre-computed equations resulting from the aes-2x2x4 sbox in GF(2^4)
BITS = 4
ROWS = 2
COLS = 2
KEYS = 10
NSTATE = ROWS*COLS*BITS

# key variables
#Nk = NSTATE*(KEYS-1)+2*BITS     # only 8 bits sbox output for the last key ie.e the last key does not have 8 bit sbox input
Nk = NSTATE*KEYS     # only 8 bits sbox output for the last key ie.e the last key does not have 8 bit sbox input

N=Nk

# ring
P=BooleanPolynomialRing(N,'x',order='lex')
#P = PolynomialRing(GF(2), N, 'x',order='lex')
# keys
k = []
for r in range(0,KEYS):
    k.append([])
    k[r].append([])
    for i in range((NSTATE*r)+0,(NSTATE*r)+NSTATE/2):
        k[r][0].append(P.gen(i))
    k[r].append([])
    for i in range((NSTATE*r)+NSTATE/2,(NSTATE*r)+NSTATE):
        k[r][1].append(P.gen(i))

# field of fractions
R.<z> = P.fraction_field()[]
#rijndael polynomial for GF16
rp = z^4+z+1

#load "leaks-2x2x4.sage"
load "leaked-bits-precomputed.sage"

# {equations from the sbox - used to compute cheat vectors of the key schedule---

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

# ---equations from the sbox - used to compute cheat vectors of the key schedule}

# add two BITS-bit strings in GF(2) ie. xor their elements
def add(a,b):

    c = [None]*BITS
    # constant as polygon
    for i in range(0,BITS):
        c[i] = a[i]+b[i]
    return c

# precomputed round constant for the key schedule (computed with rcon4() from lexs-compute-2x2x4.sage)
rc = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

# normal key schedule; K is the initial key
def ks(key,rc):
    ek = []                     # extended key

    ek2 = []

    # generate round constants
    #rc = rcon4()

    # copy the initial key as the first element of the extended key
    ek.append([])
    ek2.append([])
    for i in range(0,NSTATE):
        ek[0].append(key[i])
        ek2[0].append(key[i])

    print "ks: ek",ek
    for i in range(1,10): # 9 keys

        s0 = invt4(ek[i-1][(NSTATE-BITS):NSTATE]) # k_{i-1},1
        t0 = afft4(s0,m4,v4)                      # L(s0) + d
     
        #print "sb-output-0[",i,"]",t0
        
        s1 = invt4(ek[i-1][(NSTATE-2*BITS):NSTATE-BITS]) # k_{i-1},0
        t1 = afft4(s1,m4,v4)                             # L(s1) + d

        #print "sb-output-1[",i,"]",t1
        #print "ek",t0,t1
    
        ek.append([])
        ek2.append([])
        # [ki00, ki10  |ki01  ,ki11  ] | [k_{i+1}00,k_{i+1}10|k_{i+1}01,k_{i+1}11]
        # [BITS, 2*BITS|3*BITS,4*BITS]
    
        # calculate first element: ki00,ki10
        t0 = add(t0,rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])      # L(s0) + d + rc_i + k_{i-1},00
        ki10 = add(t1,ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
                    
        # calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                
    
        #print "ek[",i,"]",ki00,ki10,ki01,ki11
        #print "ek[",i,"]",t0,t1,ki01,ki11

        ek[i] = ki00 + ki10 + ki01 + ki11
        ek2[i] = t0 + t1 + ki01 + ki11
        print "ek[",i,"]",ek[i]

    return ek

# (22+8)*9 equations for 10 keys = 270 equations for the 10-key key schedule
# 22 sbox equations * 9 times + 8 linear equations * 9 times
# building equations from the new key schedule based on the new variable arrangement (only 8 unknown bits)
def ks2_equations(k,rc):
    ek = []                      # extended key
    keq = []

    # copy the initial key as the first element of the extended key
    ek.append(k[0][0]+k[0][1])

    NK = len(k)
    
    for i in range(1,NK):
        ek.append([])

        # (sb-output,sb-input) - non-linear equatios
        nle1 = sb4_eqs(k[i-1][1][BITS:2*BITS],k[i][0][0:BITS])
        nle2 = sb4_eqs(k[i-1][1][0:BITS],k[i][0][BITS:2*BITS])
        
        # calculate first element: ki00,ki10
        t0 = add(k[i][0][0:4],rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])                # L(s0) + d + rc_i + k_{i-1},00
        ki10 = add(k[i][0][4:8],ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
                    
        # calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                

        # (sb-output,sb-input) - 2x4 linear equations
        le1 = []
        le2 = []
        for j in range(0,BITS):
            le1.append(ki01[j] + k[i][1][j])
        for j in range(0,BITS):
            le2.append(ki11[j] + k[i][1][BITS+j])

        keq = keq + nle1 + nle2 + le1 + le2

        # the following line computes the round keys
        ek[i] = ki00 + ki10 + ki01 + ki11
        #ek[i] = k[i][0][0:BITS] + k[i][0][BITS:2*BITS] + ki01 + ki11

    return keq

# new key schedule based on the new variable arrangement (only 8 unknown bits)
# the input key k is 10 round keys according to the DC-B arranegment; the function returns the actual round keys with which encryption is done
# this function is: from-DCB-arrangement--to-normal-arrangement
def ks2_dcb_to_normal(k,rc):
    ek = []                      # extended key
    keq = []
    nml_k=[]                    # real (normal) round keys with which encrytpion is done

    # copy the initial key as the first element of the extended key
    ek.append(k[0][0]+k[0][1])
    #nml_k.append(k[0][0]+k[0][1])
    nml_k.append([])
    nml_k[0].append([])
    nml_k[0][0] = k[0][0]
    nml_k[0].append([])
    nml_k[0][1] = k[0][1]

    NK = len(k)
    #for i in range(1,KEYS):
    for i in range(1,NK):
        ek.append([])

        # NO SBOXES!!! - in dcb the words are already sboxed

        # calculate first element: ki00,ki10
        t0 = add(k[i][0][0:BITS],rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])                # L(s0) + d + rc_i + k_{i-1},00
        ki10 = add(k[i][0][BITS:2*BITS],ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10

        # Calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11 

        # the following line computes the round keys
        #ek[i] = ki00 + ki10 + k[i][1][0:BITS] + k[i][1][BITS:2*BITS]
        ek[i] = ki00 + ki10 + ki01 + ki11

        #print "ek3[",i,"]",ek[i]
        nml_k.append([])
        nml_k[i].append([])
        #nml_k[i][0] = ki00 + ki10
        nml_k[i][0] = ki00 + ki10
        nml_k[i].append([])
        #nml_k[i][1] = ki01 + ki11
        nml_k[i][1] = k[i][1][0:BITS] + k[i][1][BITS:2*BITS]

    return nml_k

# compute the actual values of the extended key-2 (new key schedule based on the new variable arrangement (only 8 unknown bits))
# the input key k is *the* 10 round keys with which encryption is done
# this function is: from-normal-arrangement--to-DCB-arrangement
#def ks2_values(k,rc):
def ks2_normal_to_dcb(k,rc):
    ek = []                      # extended key
    keq = []
    dcb_k=[]                    # key according to DeCaniere-Biruykov (DC-B) arrangement

    # copy the initial key as the first element of the extended key
    ek.append(k[0][0]+k[0][1])
    #dcb_k.append(k[0][0]+k[0][1])
    dcb_k.append([])
    dcb_k[0].append([])
    dcb_k[0][0] = k[0][0]
    dcb_k[0].append([])
    dcb_k[0][1] = k[0][1]

    NK = len(k)
    #for i in range(1,KEYS):
    for i in range(1,NK):
        ek.append([])
        #dcb_k.append([])

        # the first 4 bits of current key are the sbox outputs of the last four bits of the previous key
        s0 = sbox4(k[i-1][1][BITS:2*BITS])
        # the second 4 bits of current key are the sbox outputs of the last but one four bits of the previous key
        s1 = sbox4(k[i-1][1][0:BITS])

        # calculate first element: ki00,ki10
        t0 = add(s0,rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])      # L(s0) + d + rc_i + k_{i-1},00
        ki10 = add(s1,ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10

        # Calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11 

        # the following line computes the round keys
        ek[i] = ki00 + ki10 + ki01 + ki11

        # the following line computes the values of the DeCanniere-Biruykov arrangement
        dcb_k.append([])
        dcb_k[i].append([])
        dcb_k[i][0] = s0 + s1
        dcb_k[i].append([])
        dcb_k[i][1] = ki01 + ki11

        #print "ek_values[",i,"]",ek[i]

    #fname = 'dcb_k.sage'  
    #print "fname",fname
    #f = open(fname,"w") 
    #f.write("actual values of the extended key-2\n")
    #f.write("\n\ndcb_k = ")
    #f.write(str(dcb_k))
    #f.close()
    return dcb_k

# add equations for the initial key and find the solutions to the rest 
def ks2_test1():
    ek,e=ks2(k,rc)
    # cheat: add 8 bits from the real key 
    for i in range(0,NSTATE/2):
        e.append(ekey[0][i] + k[0][0][i])
    for i in range(NSTATE/2,NSTATE):
        e.append(ekey[0][i] + k[0][1][i-NSTATE/2])
    print "\nextended key\n",ek
    print "\nkeqs\n",e
    
    I=ideal(e)
    d=I.dimension()
    print "ideal dimension",d
    G=I.groebner_basis()
    I2=ideal(G)
    V=I2.variety()
    print "V",V


def write_ks2_eqs_to_file(k,rc):

    keq = ks2_equations(k,rc)

    NN = NSTATE*len(k)

    s = str(keq)
    # replace "x12*x16 +" with "x12*k[16] +"
    for i in range(0,NN):
        s1 = 'x'+str(i)+' +'
        s2 = 'k['+str(i)+'] +'
        s = s.replace(s1,s2)
    # replace "x12*x16 +" with "k[12]*x16 +"
    for i in range(0,NN):
        s1 = 'x'+str(i)+'*'
        s2 = 'k['+str(i)+']*'
        s = s.replace(s1,s2)
    # replace "x12," with "+ k[12],"
    for i in range(0,NN):
        s1 = 'x'+str(i)+','
        s2 = 'k['+str(i)+'],'
        s = s.replace(s1,s2)
    # replace x159] with k[159]]
    s1 = 'x'+str(NN-1)+']'
    s2 = 'k['+str(NN-1)+']]'
    s = s.replace(s1,s2)

    fname = 'ks2-eqs-2x2x4-'+str(len(k))+'-keys.sage'  
    print "fname",fname
    f = open(fname,"w")
    f.write("# explicit equations for the key schedule 2 of lex-2x2x4 for [")
    f.write(str(len(k)))
    f.write("] keys \n")
    f.write("# the key variables are input and output of sboxes \n")
    f.write("# the input key x is expected to be in the following format \n")
    f.write("# x = " + str(k) + "\n")
    f.write("def ks2_eqs(k):\n") 
    f.write("    e = ")
    f.write(s)
    f.write("\n")
    f.write("    return e")
    f.close()


def write_dcb_round_keys_to_file(k,rc):

    nml_k = ks2_dcb_to_normal(k,rc)

    NN = NSTATE*len(k)

    s = str(nml_k)
    # replace "x12*x16 +" with "x12*k[16] +"
    for i in range(0,NN):
        s1 = 'x'+str(i)+' +'
        s2 = 'k['+str(i)+'] +'
        s = s.replace(s1,s2)
    # replace "x12," with "k[12],"
    for i in range(0,NN):
        s1 = 'x'+str(i)+','
        s2 = 'k['+str(i)+'],'
        s = s.replace(s1,s2)
    # replace x159] with k[159]]
    for i in range(0,NN):
        s1 = 'x'+str(i)+']'
        s2 = 'k['+str(i)+']]'
        s = s.replace(s1,s2)
    #s1 = 'x'+str(NN-1)+']'
    #s2 = 'k['+str(NN-1)+']]'
    #s = s.replace(s1,s2)

    fname = 'ks2-roundkeys-2x2x4-'+str(len(k))+'-keys.sage'  
    print "fname",fname
    f = open(fname,"w")
    f.write("# round keys according to DeCanniere-Biruykov arrangement of lex-2x2x4 for [")
    f.write(str(len(k)))
    f.write("] keys \n")
    f.write("# the key variables are input and output of sboxes \n")
    f.write("# the input key x is expected to be in the following format \n")
    f.write("# x = " + str(k) + "\n")
    f.write("# the output is the round keys added at the end of each round according to DCB arrangement\n") 
    f.write("def ks2_round_keys(k):\n") 
    f.write("    e = ")
    f.write(s)
    f.write("\n")
    f.write("    return e")
    f.close()

def ks2_test2():

    # input: normal key; output: dcb key
    dcb_k = ks2_normal_to_dcb(ekey,rc)
    print "dcb_k",dcb_k
    # input: dcb key; output: normal key
    nml_k = ks2_dcb_to_normal(ekey2,rc)
    print "nml_k",nml_k
    print "dcb_k==ekey2",dcb_k==ekey2
    print "nml_k==ekey",nml_k==ekey
    assert dcb_k==ekey2
    assert nml_k==ekey
    e=ks2_equations(ekey2,rc)
    print "e",e
    print "e==[0]",e==[0]*len(e)
    assert e==[0]*len(e)

#
#ks2_test2()
#nml_k = ks2_dcb_to_normal(k,rc)

# generate equations for 2,3,...,10 round keys
for nk in range(2,KEYS+1):
    kk = k[0:nk]
    write_ks2_eqs_to_file(kk,rc)
    write_dcb_round_keys_to_file(kk,rc)

