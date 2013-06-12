# 
# Original file: lin-compute-2x2x4.sage
# 

#
# computing equations for the linear part in one round (ShiftRows+MixColumns+AddRoundKey) of lexs 2x2x4 according to the DeCanniere-Biruykov arrangement of the round variables (see next)
#
# used lexs-compute-2x2x4.sage, sbox-compute-2x2x4.sage, sbox-solve-2x2x4.sage, ks-compute-2x2x4.sage
#
# the ciphertext variables are organized as follow:
#
# x[0] = [initial_input[0:16]==sbox-input[0:16]]
# x[1] = [sbox-output[0:16],sbox-input[0:16]]
# x[2] = [sbox-output[0:16],sbox-input[0:16]]
# ...
# x[R] = [sbox-output[0:16],sbox-input[0:16]]
# ...
#
# for one round we have 32+16=48 variables (x[0] = [sbox-input[0:16],sbox-output[0:16]] + x[1] = [sbox-input[0:16]])
# every (sbox-input[0:16],sbox-output[0:16]) pair results in 4*11=44 non-linear equations
# every (sbox-output[0:16],sbox-input[0:16]) pair results in 16 linear equations
# thus for one round we have: 48 variables and 44+16=60 equations
#
#load "sbox-eqs-2x2x4.sage"
load "sbox-eqs-precomputed.sage"

BITS = 4
ROWS = 2
COLS = 2
KEYS=2
ROUNDS = 1
NSTATE = ROWS*COLS*BITS

# ciphertext variables for R rounds: initial ab-input + (ab-output and sb-input for this round)
Nk = NSTATE*KEYS     # only 8 bits sbox output for the last key ie.e the last key does not have 8 bit sbox input
Nc = NSTATE + ROUNDS*(NSTATE + NSTATE)

# total number of variables
N = Nk+Nc

# ring
P = BooleanPolynomialRing(N,'x',order='lex')

#load "leaks-2x2x4.sage"
load "leaked-bits-precomputed.sage"

# field of fractions
R.<z> = P.fraction_field()[]
#rijndael polynomial for GF16
rp = z^4+z+1

# keys
k = []
# first KEYS-1
for r in range(0,KEYS):
    k.append([])
    k[r].append([])
    for i in range((NSTATE*r)+0,(NSTATE*r)+NSTATE/2):
        k[r][0].append(P.gen(i))
    k[r].append([])
    for i in range((NSTATE*r)+NSTATE/2,(NSTATE*r)+NSTATE):
        k[r][1].append(P.gen(i))
# ciphertexts
x=[]
# initial input
for r in range(0,1):
    x.append([])
    x[r].append([])             # x[r][0]=[]
    x[r].append([])             # x[r][1]
    for i in range(Nk,Nk+NSTATE):
        x[r][1].append(P.gen(i))
# the rest inputs
for r in range(1,ROUNDS+1):
    x.append([])
    x[r].append([])
    for i in range(Nk+NSTATE+2*NSTATE*(r-1),Nk+NSTATE+2*NSTATE*(r-1)+NSTATE):
        x[r][0].append(P.gen(i))
    x[r].append([])
    for i in range(Nk+NSTATE+2*NSTATE*(r-1)+NSTATE,Nk+NSTATE+2*NSTATE*(r-1)+2*NSTATE):
        x[r][1].append(P.gen(i))

# ShiftRows4
# left shift constants for each row
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

def mc4(in_state):
    # state is 2x2
    mz2 = matrix ([[z+1,z],
                   [z,z+1]])
    out_state = []
    # list of polynomials in z
    zpolys = []
    # represent each byte of the state as a polynomial in z
    # for 2x2 the order should be: 
    # zpoly[0] = in_state[3]*z^3  + in_state[2]*z^2  + in_state[1]*z^1  + in_state[0]*z^0
    # zpoly[1] = in_state[7]*z^3  + in_state[6]*z^2  + in_state[5]*z^1  + in_state[4]*z^0
    # zpoly[2] = in_state[11]*z^3 + in_state[10]*z^2 + in_state[9]*z^1  + in_state[8]*z^0
    # zpoly[3] = in_state[15]*z^3 + in_state[14]*z^2 + in_state[13]*z^1 + in_state[12]*z^0
    for ibyte in range(0,ROWS*COLS):
        ipoly = 0
        for ibit in range(0,BITS):
            ipoly+= (in_state[ibyte*BITS+ibit])*z^ibit
        zpolys.append(ipoly)
    # put the polynomials of the state in 2x2 matrix sz
    mz = mz2
    sz = matrix([[zpolys[0], zpolys[1]],
                 [zpolys[2], zpolys[3]]])
    # multiply matrices modulo the rijndael polynomial
    mc = (mz*sz) % rp
    # put the coefficients of the polynomials of mc in the output state
    for irow in range(0,ROWS):
        for ibyte in range(0,COLS):
            mcpoly = mc[irow][ibyte].list() # WARNING!
            # get the degree of the polynomial; if it is less than 7, padd the missing coeffs with zeros
            degree = mc[irow][ibyte].degree()     
            assert degree < BITS
            # pad the "missing" coefficients from degree+1 up to 7 with zeros
            for icoeff in range(degree+1,BITS):
                mcpoly.append(0)
            for ibit in range(0,BITS):
                out_state.append(mcpoly[ibit])
    return out_state

# AddRoundKey4
def ark4(in_state,rk):
    out_state = []
    # add round key
    for ibit in range(0,ROWS*COLS*BITS):
        out_bit = in_state[ibit] + rk[ibit]
        out_state.append(out_bit)
    return out_state

# linear part of the round transformation (SR+MC+ARK)
def r4_lin(x, rk):

    #print "r4: sx",sx
    rx = sr4(x)
    #print "r4: rx",rx
    mx = mc4(rx)
    #print "r4: mx",mx
    y = ark4(mx,rk)
    #print "r4: y",y
    return y

# precomputed round constant for the key schedule (computed with rcon4() from lexs-compute-2x2x4.sage)
rc = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

# add two BITS-bit strings in GF(2) ie. xor their elements
def add(a,b):

    c = [None]*BITS
    # constant as polygon
    for i in range(0,BITS):
        c[i] = a[i]+b[i]
    return c

# new key schedule based on the new variable arrangement (only 8 unknown bits) (used in ks-compute-2x2x4.sage)
def ks2(k,rc):
    ek = []                      # extended key
    #kleq = []                    # key linear equations
    #knleq = []                   # key non-linear equations
    keq = []

    # copy the initial key as the first element of the extended key
    ek.append(k[0][0]+k[0][1])
    #print "ek",ek
    for i in range(1,KEYS):
        ek.append([])
        # calculate first element: ki00,ki10
        t0 = add(k[i][0][0:4],rc[i])                  # L(s0) + d + rc_i
        ki00 = add(t0,ek[i-1][0:BITS])      # L(s0) + d + rc_i + k_{i-1},00
        #print "ki00",ki00
        ki10 = add(k[i][0][4:8],ek[i-1][BITS:2*BITS]) # L(s1) + d + k_{i-1},10
        #print "ki10",ki10
                    
        # calculate second element: ki01,ki11
        ki01 = add(ki00,ek[i-1][2*BITS:3*BITS]) # L(s0) + d + rc_i + k_{i-1},00 + k_{i-1},01 
        #print "ki01",ki01
        ki11 = add(ki10,ek[i-1][3*BITS:4*BITS]) # L(s1) + d + k_{i-1},10 + k_{i-1},11                
        #print "ki11",ki11

        ek[i] = ki00 + ki10 + ki01 + ki11
        #print "ek",ek
    return ek

# constructing the linear equations
def lineq(x,ek):
    e = []
    for r in range(1,ROUNDS+1):

        x0 = x[r-1][1] # output from round 0 = input to the sbox of round 1
        x1 = x[r][0] # output from the sbox of round 1
        x2 = x[r][1] # output from round 1 == input to the sbox of round 2
        rk = ek[r]

        sbeq = sb4_eqs(x0,x1)
        print "sbeq",sbeq

        print "\n"
        print "x0[",r-1,"]",x0
        print "x1[",r,"]",x1
        print "x2[",r,"]",x2
        print "rk[",r,"]",rk
        print "\n"

        y = r4_lin(x1,rk)       # output from the round
        print "y",y

        for i in range(0,NSTATE):
            e.append(x2[i]+y[i])
    print "e",e
    return e

# xin - input to the linear round; xout - output from the linear round, key - round key to the round
#def write_eqs_to_file(xin,xout,key):
def write_eqs_to_file():

    e = []
    #x0 = x[0][1] # output from round 0 = input to the sbox of round 1
    x1 = x[1][0] # output from the sbox of round 1
    x2 = x[1][1] # output from round 1 == input to the sbox of round 2
    rk= ek[0]

    print "rk",rk
    #print "x0",x0
    print "x1",x1
    print "x2",x2

    y = r4_lin(x1,rk)       # output from the round
    for i in range(0,NSTATE):
        e.append(x2[i]+y[i])
    
    print "e",e

    s = str(e)
    # replace x[0:15] with k[0:15]
    for i in range(0,NSTATE):
        s1 = 'x'+str(i)+' '
        s2 = 'k['+str(i)+'] '
        s = s.replace(s1,s2)

    # replace x[32:47] with x[16:31]
    for i in range(2*NSTATE,2*NSTATE+NSTATE):
        s1 = 'x'+str(i)+' '
        s2 = 'x['+str(i-2*NSTATE)+'] '
        s = s.replace(s1,s2)

    # replace x[48:63] with x[32:47]
    for i in range(3*NSTATE,3*NSTATE+NSTATE):
        s1 = 'x'+str(i)
        s2 = 'y['+str(i-3*NSTATE)+']'
        s = s.replace(s1,s2)

    print "s",s

    fname = 'lin-eqs-2x2x4.sage'  
    print "fname",fname

    f = open(fname,"w")
    f.write("le = ")
    f.write(s)
    f.close()

    return s

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
#s=write_eqs_to_file()

# cheat sheets for testing
#k0 = [P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0)] # initial key
#k1 = [P(1),P(1),P(1),P(1),P(1),P(1),P(1),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(0)]
# solutions of the key equations for the two keys k0 and k1
#p=[None]*32;p[31]=P(0);p[30]=P(0);p[29]=P(0);p[28]=P(0);p[27]=P(1);p[26]=P(1);p[25]=P(1);p[24]=P(0);p[23]=P(1);p[22]=P(0);p[21]=P(1);p[20]=P(1);p[19]=P(1);p[18]=P(0);p[17]=P(1);p[16]=P(0);p[15]=P(0);p[14]=P(1);p[13]=P(1);p[12]=P(1);p[11]=P(0);p[10]=P(0);p[9]=P(0);p[8]=P(1);p[7]=P(1);p[6]=P(1);p[5]=P(0);p[4]=P(0);p[3]=P(0);p[2]=P(1);p[1]=P(1);p[0]=P(1)
# plaintext and leaks for two rounds
#pt = [P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1)] # plaintext (ie. x,IV,etc)
#l0 = [P(0),P(1),P(0),P(1),P(0),P(0),P(0),P(0),P(0),P(0),P(1),P(1),P(0),P(0),P(1),P(1)]
#l1 = [P(1),P(1),P(0),P(1),P(1),P(1),P(0),P(1),P(0),P(1),P(1),P(0),P(1),P(0),P(0),P(1)]
# cheat sheet - snoxes on the inputs of two rounds
#sage: load "/lien/vvesseli/in/box/work/src/lex/sbox-compute-2x2x4.sage"
#spt = [P(1), P(0), P(1), P(1), P(0), P(1), P(0), P(0), P(1), P(0), P(1), P(0), P(1), P(0), P(1), P(1)]
#sl0 = [P(1), P(1), P(1), P(1), P(0), P(1), P(1), P(0), P(1), P(1), P(0), P(0), P(1), P(1), P(0), P(0)]
#sl1 = [P(0), P(0), P(1), P(1), P(0), P(0), P(1), P(1), P(1), P(1), P(1), P(0), P(1), P(0), P(1), P(1)]

# assign the cheat values
#k[0][0]=p[0:8]
#k[0][1]=p[8:16]
#k[1][0]=p[16:24]
#k[1][1]=p[24:32]
#x[0][1]=l0
#x[1][0]=sl0
#x[1][1]=l1

k = ekey2[:KEYS]

x[0][1]=l[0]
x[1][0]=sbl[0]
x[1][1]=l[1]


#x[0][1]=pt
#x[1][0]=spt
#x[1][1]=l0
#x[2][0]=sl0
#x[2][1]=l1
#x[3][1]=l2

#nml_k = ks2_dcb_to_normal(k,rc)
ek=ks2(k,rc)
print "k",k
print "x",x
print "ek",ek

e=lineq(x,ek)
print "e",e


