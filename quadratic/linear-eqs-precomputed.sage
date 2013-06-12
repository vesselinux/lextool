# 
# Original file: lin-eqs-2x2x4.sage
# 

# pre-computed equations for the linear part (SR+MC+ARK) of one round of lexs-2x2x4 ; computed with file linear-eqs-compute.sage (old name lin-compute-2x2x4.sage)
# x - input to the linear part of the round (ie. ouput from the SBoxes)
# y - output from the round transformation (ie. input to the next round)
# k - round key
# the final equations are: x+k=y or equivalently x+k-y=0
def lin4_eqs(k,x,y):
    e = [k[0] + x[0] + x[3] + x[15] + y[0], k[1] + x[0] + x[1] + x[3] + x[12] + x[15] + y[1], k[2] + x[1] + x[2] + x[13] + y[2], k[3] + x[2] + x[3] + x[14] + y[3], k[4] + x[4] + x[7] + x[11] + y[4], k[5] + x[4] + x[5] + x[7] + x[8] + x[11] + y[5], k[6] + x[5] + x[6] + x[9] + y[6], k[7] + x[6] + x[7] + x[10] + y[7], k[8] + x[3] + x[12] + x[15] + y[8], k[9] + x[0] + x[3] + x[12] + x[13] + x[15] + y[9], k[10] + x[1] + x[13] + x[14] + y[10], k[11] + x[2] + x[14] + x[15] + y[11], k[12] + x[7] + x[8] + x[11] + y[12], k[13] + x[4] + x[7] + x[8] + x[9] + x[11] + y[13], k[14] + x[5] + x[9] + x[10] + y[14], k[15] + x[6] + x[10] + x[11] + y[15]]
    return e
