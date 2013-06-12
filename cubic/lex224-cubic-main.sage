# 
# Generate quadratic equations for Lex(2,2,4) and solves them using Groebner bases to recover the secret key of the cipher. 
# Original name: lexs-2x2x4.sage
# 

# 1 aes round is composed of one application of (SubBytes+ShiftRows+MixColumns+AddRoundKey)
AES_ROUNDS = 9

ROUNDS = 5


#UNKNOWN=0
#UNKNOWN=1
#UNKNOWN=2
#UNKNOWN=3
#UNKNOWN=4
#UNKNOWN=5
#UNKNOWN=6
#UNKNOWN=7
# ---
#UNKNOWN=8 #8,9,10,... ; must be >7
#UNKNOWN=8
UNKNOWN=9
#UNKNOWN=10
#UNKNOWN=11
#UNKNOWN=12
# ---
#UNKNOWN=13
#UNKNOWN=14
#UNKNOWN=15
#UNKNOWN=16

# a test version of small lex 2x4x4 which uses ready explicit equations for the round transformation and the key schedule
t_total_cpu = cputime()
t_total_wall = walltime()

BITS = 4
ROWS = 2
COLS = 2
NSTATE = ROWS*COLS*BITS
KEYS=2

Nk = NSTATE*KEYS

# ciphertext variables x (we add 1 for x[0] - the plaintext)
Nc = NSTATE + NSTATE*ROUNDS

# total number of variables
N = Nk + Nc

# create variables for a ring of N generators
def set_ring(N):

    # ring
    P = PolynomialRing(GF(2), N, 'x',order='degrevlex')
    #P = PolynomialRing(GF(2), N, 'x',order='lex')

    # keys
    k = []
    for r in range(0,KEYS):
        k.append([])
        for i in range((NSTATE*r)+0,(NSTATE*r)+NSTATE):
            k[r].append(P.gen(i))

    # ciphertexts (x[0] is the plaintext)
    x = []
    for r in range(0,ROUNDS+1):
        x.append([])
        for i in range((Nk+NSTATE*r)+0,(Nk+NSTATE*r)+NSTATE):
            x[r].append(P.gen(i))    

    return P,k,x

# generate field equations
def field_equations(k,x):
    fe=[]
    for i in range(0,KEYS):
        for j in range(0,NSTATE):
            fe.append(k[i][j]^2-k[i][j])
    for i in range(0,ROUNDS+1):
        for j in range(0,NSTATE):
            fe.append(x[i][j]^2-x[i][j])
    return fe

# generate key equations containing the relations between all round keys
def key_equations(k):
    ke=[]
    for i in range(0,KEYS-1):
        ke = ke + keqs(k[i],k[i+1])
    return ke

def round_equations(k,x):
    re=[]
    # round equations
    for r in range(0,ROUNDS):
        # round eqations x0+k0=x1
        re = re + reqs(x[r],k[r%KEYS],x[r+1])
    return re

# ring
P,k,x = set_ring(N)

#rijndael polynomial for GF16
R.<z> = P.fraction_field()[]
rp = z^4+z+1

# cipher equations for 1 round # def reqs(x,k,y)
#load "eqss-2x2x4.sage"          
load "round-eqs-precomputed.sage"          
# key schedule equations for 1 round # keqs(k1,k2)
#load "ks-2x2x2x4.sage"
load "keysched-eqs-precomputed.sage"
# leaks for 10 rounds
#load "leaks-2x2x4.sage"
#load "leaked-bits-precomputed.sage"
# Compute the leaks for given key and plaintext
load "leaked-bits.sage"

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

# Generate the full leak for the given key and plaintext
l = lex_leaks(pt,key,16)

# initial plaintext x[0]
#pt = [P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1),P(0),P(0),P(1)] # plaintext (ie. x,IV,etc)
# two keys (KEYS==2)
#k0 = [P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(1),P(1),P(1),P(0)] # initial key
#k1 = [P(1),P(1),P(1),P(1),P(1),P(1),P(1),P(0),P(0),P(1),P(1),P(1),P(0),P(0),P(0),P(0)]
# three leaks (ROUNDS==3)
#l0 = [P(0),P(1),P(0),P(1),P(0),P(0),P(0),P(0),P(0),P(0),P(1),P(1),P(0),P(0),P(1),P(1)]
#l1 = [P(1),P(1),P(0),P(1),P(1),P(1),P(0),P(1),P(0),P(1),P(1),P(0),P(1),P(0),P(0),P(1)]
#l2 = [P(1),P(0),P(0),P(1),P(1),P(0),P(1),P(0),P(0),P(1),P(1),P(1),P(1),P(0),P(1),P(0)] 

# main
LEAK=NSTATE-UNKNOWN
print "\nLEX-CU: ROUNDS",ROUNDS,"ROWS",ROWS,"COLS",COLS,"BITS",BITS,"KEYS",KEYS,"LEAK",LEAK
#print "LEAK:",LEAK,"bits (out of",ROWS*COLS*BITS,"\n"   #,") initial leak", INITIAL_LEAK,"\n"

# solutions: for tests
#x[0]=pt
#x[1]=l0
#x[2]=l1
#x[3]=l2
#k[0]=k0
#k[1]=k1

#print "k",k
#print "x",x

# field equations ki^2=ki xi^2=x1 
fe = field_equations(k,x)
print "field equations",len(fe)
#print "fe",fe,"\n"
# key equations k_i=k_{i+1}
ke = key_equations(k)
print "key equations",len(ke)
#print "ke",ke,"\n"
# round equations
re = round_equations(k,x)
print "round equations",len(re)
#print "re",re,"\n"

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

# leak equations
le = []
# add leaks for the initial plaintext - this will be "the bits of the leak after round 0" ie. EVEN leak

assert (LEAK < 12)

print "even leak"

# leak first 4 bits
for j in range(4,8): # x[4:7]
    le.append(x[0][j]+pt[j])
    print j

# leak next 4 bits
if LEAK > 4 and LEAK < 8:
    remaining = LEAK-4
    #for j in range(12,12+remaining): # x[12:15]
    for j in range(8,8+remaining): # x[8:11]
        le.append(x[0][j]+pt[j])
        print j

if LEAK == 8:
    #for j in range(12,16): # x[12:15]
    for j in range(8,12): # x[8:12]
        le.append(x[i+1][j]+l[i][j])
        print j

if LEAK > 8:
    #for j in range(12,16): # x[12:15]
    for j in range(8,12): # x[8:12]
        le.append(x[0][j]+pt[j])
        print j

    remaining = LEAK-8
    #for j in range(8,8+remaining): # x[8:11]
    for j in range(12,12+remaining): # x[12:15]
        le.append(x[0][j]+pt[j])
        print j

# inputs to rounds 0,1,2,...,ROUNDS
for i in range(0,ROUNDS):

    # ODD leak : x[0:3] | x[12:15]
    if (i+1)%2 == 1:          
        print "odd leak"
        # leak first 4 bits
        for j in range(0,4): # x[0:3]
            le.append(x[i+1][j]+l[i][j])
            print j

        # leak next 4 bits
        if LEAK > 4 and LEAK < 8:
            remaining = LEAK-4
            #for j in range(8,8+remaining): # x[8:11]
            for j in range(12,12+remaining): # x[12:15]
                le.append(x[i+1][j]+l[i][j])
                print j

        if LEAK == 8:
            #for j in range(8,12): # x[8:11]
            for j in range(12,16): # x[12:15]
                le.append(x[i+1][j]+l[i][j])
                print j

        if LEAK > 8:
            #for j in range(8,12): # x[8:11]
            for j in range(12,16): # x[12:15]
                le.append(x[0][j]+pt[j])
                print j

            remaining = LEAK-8
            for j in range(8,8+remaining): # x[8:11]
                le.append(x[0][j]+pt[j])
                print j

    # EVEN leak : x[4:7] | x[8:11]
    if (i+1)%2 == 0:          
        print "even leak"
        # leak first 4 bits
        for j in range(4,8): # x[4:7]
            le.append(x[i+1][j]+l[i][j])
            print j

        # leak next 4 bits
        if LEAK > 4 and LEAK < 8:
            remaining = LEAK-4
            #for j in range(12,12+remaining): # x[12:15]
            for j in range(8,8+remaining): # x[8:11]
                le.append(x[i+1][j]+l[i][j])
                print j

        if LEAK == 8:
            #for j in range(12,16): # x[12:15]
            for j in range(8,12): # x[8:12]
                le.append(x[i+1][j]+l[i][j])
                print j

        if LEAK > 8:
            #for j in range(12,16): # x[12:15]
            for j in range(8,12): # x[8:12]
                le.append(x[0][j]+pt[j])
                print j

            remaining = LEAK-8
            #for j in range(8,8+remaining): # x[8:11]
            for j in range(12,12+remaining): # x[12:15]
                le.append(x[0][j]+pt[j])
                print j
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

#print "le",le
print "leak equations",len(le)

# all equations
e = fe+re+ke+le
#print "e",e

print "total equations",len(e)
print "total equations excl. field equations ne=",len(e)-len(fe)
print "total unknowns u=",P.ngens()
print "ratio (ne/u)",float(len(e)-len(fe))/float(P.ngens())
#print "\n"

I=ideal(e)
print "calculating groebner basis..."

t_groebner = cputime()
G=I.groebner_basis()
t_groebner = cputime(t_groebner)
print "calculated groebner basis in",t_groebner,"seconds"

#print "G",G
I2=ideal(G)
print "calculating variety..."

t_variety = singular.cputime()
V = I2.variety()
t_variety = singular.cputime(t_variety)
print "calculated variety in",t_variety,"seconds"


print "\nLEX-CU: ROUNDS",ROUNDS,"ROWS",ROWS,"COLS",COLS,"BITS",BITS,"KEYS",KEYS,"LEAK",LEAK
legal_leak=4
u = P.ngens() - (ROUNDS+1)*legal_leak
guessed_bits = (ROUNDS+1)*(LEAK-legal_leak)
print "unknowns",u,"equations",len(e)-len(fe)-guessed_bits,"guessing",guessed_bits,"bits" # not guessing initial plaintext x[0]
ns = len(V)
print "solutions",ns
print "exhaustive search on the unknowns [2^",u,"]"
print "exhaustive search on key+iv [2^",2*NSTATE,"]"
print "exhaustive search on key [2^",NSTATE,"]"
print "current effort",ns,"*2^",guessed_bits,"= [2^", (log_b(ns,2) + guessed_bits),"]"

fname = 'loglex.txt'  
print "results stored in",fname
f = open(fname,"w")
s = "\nLEX-CU: ROUNDS= "+str(ROUNDS)+" ROWS= "+str(ROWS)+" COLS= "+str(COLS)+" BITS= "+str(BITS)+" KEYS= "+str(KEYS)+" LEAK= "+str(LEAK)
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
