# 
# Computes the algebraic equations for one round of lexs-2x2x4 according to De Canniere-Biruyikov; works in the field of fractions of polynomials over GF(2) modulo a small-scale rijndael polynomial
#
# 
# Original file: lexs-compute-2x2x4.sage
#
# used in lexs-compute-2x2x4.sage, sbox-compute-2x2x4.sage, sbox-solve-2x2x4.sage, ks-compute-2x2x4.sage
# used in sbox-eqs.sage, sbox-solve-2x2x4.sage, keysched-eqs.sage
#
# the ciphertext variables are organized as follow:
#
# x[0] = [initial_input[0:16]==sbox-input[0:16]==leak 0]
# x[1] = [sbox-output[0:16],sbox-input[0:16]==leak 1]
# x[2] = [sbox-output[0:16],sbox-input[0:16]==leak 2]
# ...
# x[R] = [sbox-output[0:16],sbox-input[0:16]==leak R]
# ...
#
# for one round we have 32+16=48 variables (x[0] = [sbox-input[0:16],sbox-output[0:16]] + x[1] = [sbox-input[0:16]])
# every (sbox-input[0:16],sbox-output[0:16]) pair results in 4*11=44 non-linear equations
# every (sbox-output[0:16],sbox-input[0:16]) pair results in 16 linear equations
# thus for one round we have: 48 variables and 44+16=60 equations
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
#
t_total_cpu = cputime()
t_total_wall = walltime()

# 
# What we call DeCanniere-Biruykov (DCB) arrangement refers to the following arrangements of the round keys for two rounds:
# 
# The key variables are input and output of sboxes.
# 
# The input key x is expected to be in the following format: 
# 
# x = [[[x0, x1, x2, x3, x4, x5, x6, x7], [x8, x9, x10, x11, x12, x13, x14, x15]], [[x16, x17, x18, x19, x20, x21, x22, x23], [x24, x25, x26, x27, x28, x29, x30, x31]]]
# 
# The output is the round keys added at the end of each round according to DCB arrangement
# def ks2_round_keys(k):
#     e = [[[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]], [k[8], k[9], k[10], k[11], k[12], k[13], k[14], k[15]]], [[k[0] + k[16], k[1] + k[17] + 1, k[2] + k[18], k[3] + k[19], k[4] + k[20], k[5] + k[21], k[6] + k[22], k[7] + k[23]], [k[24], k[25], k[26], k[27], k[28], k[29], k[30], k[31]]]]
#     return e
# # # 

# SR(r,2,2,4)
BITS = 4
ROWS = 2
COLS = 2
KEYS = 2

# 1 aes round is composed of one application of (SubBytes+ShiftRows+MixColumns+AddRoundKey)
AES_ROUNDS = 9

ROUNDS = 5

#
#UNKNOWN = NSTATE - LEAK \in {8,9,10,...} i.e. must be >7
#UNKNOWN=8
#UNKNOWN=10
#UNKNOWN=11
UNKNOWN=12  # <- LEAK = 4 bits i.e. the official leak
#UNKNOWN=13
#UNKNOWN=14
#UNKNOWN=15
#UNKNOWN=16

NSTATE = ROWS*COLS*BITS
LEAK=NSTATE-UNKNOWN

# ciphertext variables for R rounds: initial ab-input + (ab-output and sb-input for this round)
Nk = NSTATE*KEYS     # only 8 bits sbox output for the last key ie.e the last key does not have 8 bit sbox input
Nc = NSTATE + ROUNDS*(NSTATE + NSTATE)

# total number of variables
N = Nk+Nc

# Define the ring of Boolean polynomials
P = PolynomialRing(GF(2), N, 'x',order='degrevlex')

#rijndael polynomial for GF16
R.<z> = P.fraction_field()[]
rp = z^4+z+1

# pre-computed equations using lin-compute-2x2x4.sage for the linear part (SR+MC+ARK) of one round of lexs-2x2x4
#load "lin-eqs-2x2x4.sage"  
load "linear-eqs-precomputed.sage"  
# pre-computed 11 quadratic equations using sbox-compute-2x2x4.sage for the rijndael sbox in GF(2^4)
#load "sbox-eqs-2x2x4.sage" 
load "sbox-eqs-precomputed.sage" 
# pre-computed equations using ks2-compute-2x2x4.sage for the key schedule 2 of lex-2x2x4 for [2] keys
#load "ks2-eqs-2x2x4-2-keys.sage" 
load "keysched-eqs-precomputed.sage" 
# round keys according to DeCanniere-Biruykov arrangement of lex-2x2x4 for [2] keys 
#load "ks2-roundkeys-2x2x4-2-keys.sage" 
load "roundkeys-dcb.sage" 
# Compute the full output from ROUNDS rounds of Lex224
# for a given master key (16 bits) and plaintext (16 bits)
# From these outputs later the leaks are extracted
load "leaked-bits.sage"
# leaks for LEX: ROUNDS 20 ROWS 2 COLS 2 BITS 4 KEYS 2
#load "leaks-2x2x4.sage"  
#load "leaked-bits-precomputed.sage"  

# 
# Values of the key and the plaintext from the example given in Appendix of the paper:
# Velichkov, Rijmen, Preneel, "Algebraic Cryptanalysis of a Small-scale Version of Stream Cipher Lex", IET Journal paper
# 
pt  = [P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1)] # plaintext (ie. x,IV,etc)
key = [P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0)] # initial key

# Uncomment the folowing, to generate random key and plaintext

# Generate random plaintext
#pt  = [P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1))] 

# Generate random key
#key = [P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),
#       P(randint(0,1)),P(randint(0,1)),P(randint(0,1)),P(randint(0,1))] 

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

# transform key from 3d arrey in to a long 1d arrey
def key_to_1d_arrey(k):
    kt = []
    NK = len(k)
    #print "NK==KEYS",NK==KEYS
    for i in range(0,NK):
        for j in range(0,2):
            for l in range(0,NSTATE/2):
                kt.append(k[i][j][l])
    return kt

# lexs-2x2x4 state outline: 
# 
# x[0:3]  x[4:7]
# x[8:11] x[12:15]
# 

# note on the leaks:
# the leak after odd rounds: 1,3,5,... is  x[0:3] and x[8:11]
# the leak after even rounds: 2,4,6,... is x[4:7] and x[12:15]
# 
# case 1
# 
# ODD     | EVEN 
# --------------
# x[0:3]  | x[4:7]
# x[8:11] | x[12:15]
# 
# case 2
# 
# ODD     | EVEN 
# --------------
# x[0:3]  | x[4:7]
# x[12:15]| x[8:11]
# 

# constructing the equations
def eqs2(x,kt):
    e = []
    cheats = []
 
    nsbeq=0                     # counter of the sbox equations
    nleq=0                      # counter of the linear equations

    #print "x",x
    # key equations
    keq = ks2_eqs(kt)
    rks = ks2_round_keys(kt)
    rks1d = key_to_1d_arrey(rks)

    e = e + keq

    #full_leak = lex_leaks(pt,key,16)
    l = lex_leaks(pt,key,16)
    sbl = sbox_leaks(l)

    print "l  ", l
    print "sbl", sbl

    #f = open(fname,"a")
    #s="\n\nkey equations\n\n $"+str(latex(keq))+"$"
    #f.write(s)
    #f.write("\n")
    #f.close()

    print "key equations",len(keq)

    # {---vpv 20090119
    #for i in range(0,NSTATE-UNKNOWN):
    #    cheats.append(x[0][1][i] + l[0][i])
    # vpv 20090119---}

    assert (LEAK < 9)

    # EVEN leak : x[4:7] | x[8:11]
    # leak first 4 bits
    print "even leak round",0
    if LEAK < 4:
        for i in range(4,4+LEAK): # x[4:7]
            cheats.append(x[0][1][i] + l[0][i])
            print i
    else:
        for i in range(4,8): # x[4:7]
            cheats.append(x[0][1][i] + l[0][i])
            print i

    # leak next 4 bits
    if LEAK > 4:
        remaining = LEAK-4
        #for i in range(8,8+remaining): # x[8:11]
        for i in range(12,12+remaining): # x[12:15]
            cheats.append(x[0][1][i] + l[0][i])
            print i

    for r in range(1,ROUNDS+1):

        #x[r-1][1] = l[r-1]
        #x[r][0] = sbl[r-1]
        #x[r][1] = l[r]

        # ODD leak : x[0:3] | x[12:15]
        if r%2 == 1:          
            print "odd leak round",r
            # leak first 4 bits
            if LEAK < 4:
                for i in range(0,LEAK): # x[0:3]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i
            else:
                for i in range(0,4): # x[0:3]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i

            # leak next 4 bits
            if LEAK > 4:
                remaining = LEAK-4
                #for i in range(12,12+remaining): # x[12:15]
                for i in range(8,8+remaining): # x[8:11]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i

        # EVEN leak : x[4:7] | x[8:11]
        if r%2 == 0:          
            print "even leak round",r
            # leak first 4 bits
            if LEAK < 4:
                for i in range(4,4+LEAK): # x[4:7]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i
            else:
                for i in range(4,8): # x[4:7]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i

            # leak next 4 bits
            if LEAK > 4:
                remaining = LEAK-4
                #for i in range(8,8+remaining): # x[8:11]
                for i in range(12,12+remaining): # x[12:15]
                    cheats.append(x[r][0][i] + sbl[r-1][i])
                    cheats.append(x[r][1][i] + l[r][i])
                    print i

        #print "x[r-1][1] + l[r-1]",x[r-1][1] + l[r-1]
        #print "x[r][0] + sbl[r-1]",x[r][0] + sbl[r-1]
        #print "x[r][1] + l[r]",x[r][1] + l[r]
        #print "cheats",cheats

        x0 = x[r-1][1] # output from round 0 = input to the sbox of round 1
        x1 = x[r][0] # output from the sbox of round 1
        x2 = x[r][1] # output from round 1 == input to the sbox of round 2
        rk = key_to_1d_arrey([rks[r%KEYS]])

        # sbox equations 
        sbeq0 = sb4_eqs(x0[0*BITS:0*BITS+BITS],x1[0*BITS:0*BITS+BITS])
        sbeq1 = sb4_eqs(x0[1*BITS:1*BITS+BITS],x1[1*BITS:1*BITS+BITS])
        sbeq2 = sb4_eqs(x0[2*BITS:2*BITS+BITS],x1[2*BITS:2*BITS+BITS])
        sbeq3 = sb4_eqs(x0[3*BITS:3*BITS+BITS],x1[3*BITS:3*BITS+BITS])
        sbeq = sbeq0+sbeq1+sbeq2+sbeq3
        # print "\nround#",r,"sbox equations",sbeq
        leq = lin4_eqs(rk,x1,x2)
        # print "\nround#",r,"linear equations",leq

        e = e + sbeq + leq

        nsbeq=nsbeq+len(sbeq)
        nleq=nleq+len(leq)

        #f = open(fname,"a")
        #s="\n\nround sbox equations\n\n $"+str(latex(sbeq))+"$"
        #f.write(s)
        #s="\n\nround linear equations\n\n $"+str(latex(leq))+"$"
        #f.write(s)
        #s="\n\nleak equations\n\n $"+str(latex(cheats))+"$"
        #f.write(s)
        #f.write("\n")
        #f.close()

    # print "\n all equations",e
    print "sbox equations",nsbeq
    print "linear equations",nleq
    print "leak equations",len(cheats)
    return e + cheats

# generate field equations
def field_equations():
    fe=[]
    for i in range(0,P.ngens()):
        fe.append(P.gen(i)^2-P.gen(i))
    return fe

# extract the values of the bits of the recovered key
# from the variety; d is of type dict
def print_recovered_key(d):

    dk = d.keys()
    dv = d.values()
    key = []

    for j in range(0,16):
        x_str = 'x'+str(j)
        for i in range(0,len(dk)):
            if str(dk[i])==x_str: key.append(dv[i])

    return key


print "\nLEX-MQ: ROUNDS",ROUNDS,"ROWS",ROWS,"COLS",COLS,"BITS",BITS,"KEYS",KEYS,"LEAK",LEAK

print "plaintext", pt
print "key      ", key

fe = field_equations()
print "field equations",len(fe)
#print "fe",fe

kt = key_to_1d_arrey(k)
e = eqs2(x,kt)

e = fe + e

#print "quadratic equations stored in",fname
print "total equations",len(e)
print "total equations excl. field equations ne=",len(e)-len(fe)
print "total unknowns u=",P.ngens()
print "ratio (ne/u)",float(len(e)-len(fe))/float(P.ngens())
#print "\n"

I=ideal(e)
#print "calculating ideal dimention..."
#d=I.dimension()
#print "I dimension",d
#timeit('I.dimension()')
print "calculating groebner basis..."

t_groebner = cputime()
G=I.groebner_basis()
t_groebner = cputime(t_groebner)
print "calculated groebner basis in",t_groebner,"seconds"

I2=ideal(G)
print "calculating ideal dimention..."
d2=I2.dimension()
print "I2 dimension",d2
print "calculating variety..."

t_variety = singular.cputime()
V = I2.variety()
t_variety = singular.cputime(t_variety)
print "calculated variety in",t_variety,"seconds"


#print "\nLEX-MQ: ROUNDS",ROUNDS,"ROWS",ROWS,"COLS",COLS,"BITS",BITS,"KEYS",KEYS,"LEAK",LEAK
legal_leak=4
u = P.ngens() - (ROUNDS+1)*legal_leak
guessed_bits = (ROUNDS+1)*(LEAK-legal_leak)
#print "unknowns",u,"equations",len(e)-len(fe)-guessed_bits,"guessing",guessed_bits,"bits" # not guessing initial plaintext x[0]
ns = len(V)


print "guessing",guessed_bits,"bits" # not guessing initial plaintext x[0]
print "solutions",ns
print "exhaustive search on the unknowns [2^",u,"]"
print "exhaustive search on key+iv [2^",2*NSTATE,"]"
print "exhaustive search on key [2^",NSTATE,"]"
print "current effort",ns,"*2^",guessed_bits,"= [2^", (log_b(ns,2) + guessed_bits),"]"

fname = 'loglex.txt'  
f = open(fname,"w")
s = "\nLEX-MQ: ROUNDS= "+str(ROUNDS)+" ROWS= "+str(ROWS)+" COLS= "+str(COLS)+" BITS= "+str(BITS)+" KEYS= "+str(KEYS)+" LEAK= "+str(LEAK)
f.write(s)
s="\ntotal equations excl. field equations ne= "+str(len(e)-len(fe))
f.write(s)
s="\ntotal unknowns u= "+str(P.ngens())
f.write(s)
s="\nratio (ne/u) "+str(float(len(e)-len(fe))/float(P.ngens()))
f.write(s)
s="\nguessing bits = "+str(guessed_bits)
f.write(s)
s="\nsolutions = "+str(ns)
f.write(s)
s="\ncurrent effort = 2^"+str((log_b(ns,2) + guessed_bits))
f.write(s)
s="\ncalculated Groebner Basis in "+str(t_groebner)+" sec"
f.write(s)
s="\ncalculated the Variety in "+str(t_variety)+" sec"
f.write(s)
t_total_cpu = cputime(t_total_cpu)
s="\ntotal cpu time "+str(t_total_cpu)+" sec"
f.write(s)
t_total_wall = walltime()
s="\ntotal wall time "+str(t_total_wall)+"sec"
f.write(s)
f.write("\n")
f.close()

print "total cpu time ",t_total_cpu,"sec"
print "total wall time ",t_total_wall,"sec"
print "results stored in",fname

if ns == 1:
    rec_key = print_recovered_key(V[0])
    print "    plaintext=", pt
    print "   secret key=", key
    print "recovered key=", rec_key
    assert key == rec_key
