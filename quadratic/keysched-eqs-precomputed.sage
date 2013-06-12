# 
# Original file: ks2-eqs-2x2x4-2-keys.sage
# 
# explicit equations for the key schedule 2 of lex-2x2x4 for [2] keys 
# the key variables are input and output of sboxes 
# the input key x is expected to be in the following format 
# x = [[[x0, x1, x2, x3, x4, x5, x6, x7], [x8, x9, x10, x11, x12, x13, x14, x15]], [[x16, x17, x18, x19, x20, x21, x22, x23], [x24, x25, x26, x27, x28, x29, x30, x31]]]
def ks2_eqs(k):
    e = [k[12]*k[16] + k[12]*k[17] + k[12]*k[19] + k[12] + k[13]*k[16] + k[13]*k[17] + k[13] + k[14]*k[16] + k[14]*k[19] + k[15]*k[18] + k[15]*k[19] + k[15], k[12]*k[16] + k[12]*k[17] + k[12]*k[18] + k[13]*k[16] + k[13]*k[17] + k[13]*k[19] + k[13] + k[14]*k[16] + k[14]*k[17] + k[14] + k[15]*k[16] + k[15]*k[19], k[12]*k[17] + k[12]*k[18] + k[12]*k[19] + k[13]*k[16] + k[13]*k[17] + k[13]*k[18] + k[14]*k[16] + k[14]*k[17] + k[14]*k[19] + k[14] + k[15]*k[16] + k[15]*k[17] + k[15], k[12]*k[16] + k[12]*k[18] + k[12]*k[19] + k[13]*k[16] + k[13]*k[17] + k[13]*k[18] + k[14]*k[16] + k[14]*k[17] + k[14] + k[15]*k[18] + k[15]*k[19] + k[15], k[12]*k[16] + k[12]*k[17] + k[12]*k[19] + k[12] + k[13]*k[16] + k[13]*k[19] + k[13] + k[14]*k[19] + k[15]*k[16] + k[15]*k[18] + k[15], k[12]*k[16] + k[12]*k[17] + k[12]*k[18] + k[13]*k[16] + k[13]*k[17] + k[13] + k[14]*k[18] + k[14]*k[19] + k[15]*k[17] + k[15]*k[19] + k[15], k[12]*k[17] + k[12]*k[18] + k[12]*k[19] + k[13]*k[16] + k[13]*k[17] + k[13]*k[19] + k[13] + k[14]*k[16] + k[14]*k[19] + k[15]*k[19] + k[15], k[12]*k[17] + k[12]*k[19] + k[12] + k[13]*k[17] + k[13]*k[18] + k[13]*k[19] + k[14]*k[16] + k[14]*k[18] + k[14] + k[15]*k[16] + k[15]*k[17] + k[15]*k[18] + k[16] + k[18] + k[19] + 1, k[12]*k[16] + k[12]*k[17] + k[12]*k[18] + k[13]*k[18] + k[13] + k[14]*k[16] + k[14]*k[17] + k[14]*k[19] + k[14] + k[15]*k[17] + k[15] + k[16] + k[17] + k[19] + 1, k[12]*k[16] + k[12]*k[18] + k[12] + k[13]*k[16] + k[13]*k[17] + k[13]*k[18] + k[14]*k[18] + k[14] + k[15]*k[16] + k[15]*k[17] + k[15]*k[19] + k[15] + k[16] + k[17] + k[18], k[12]*k[17] + k[12]*k[18] + k[12]*k[19] + k[13]*k[16] + k[13]*k[18] + k[13] + k[14]*k[16] + k[14]*k[17] + k[14]*k[18] + k[15]*k[18] + k[15] + k[17] + k[18] + k[19], k[8]*k[20] + k[8]*k[21] + k[8]*k[23] + k[8] + k[9]*k[20] + k[9]*k[21] + k[9] + k[10]*k[20] + k[10]*k[23] + k[11]*k[22] + k[11]*k[23] + k[11], k[8]*k[20] + k[8]*k[21] + k[8]*k[22] + k[9]*k[20] + k[9]*k[21] + k[9]*k[23] + k[9] + k[10]*k[20] + k[10]*k[21] + k[10] + k[11]*k[20] + k[11]*k[23], k[8]*k[21] + k[8]*k[22] + k[8]*k[23] + k[9]*k[20] + k[9]*k[21] + k[9]*k[22] + k[10]*k[20] + k[10]*k[21] + k[10]*k[23] + k[10] + k[11]*k[20] + k[11]*k[21] + k[11], k[8]*k[20] + k[8]*k[22] + k[8]*k[23] + k[9]*k[20] + k[9]*k[21] + k[9]*k[22] + k[10]*k[20] + k[10]*k[21] + k[10] + k[11]*k[22] + k[11]*k[23] + k[11], k[8]*k[20] + k[8]*k[21] + k[8]*k[23] + k[8] + k[9]*k[20] + k[9]*k[23] + k[9] + k[10]*k[23] + k[11]*k[20] + k[11]*k[22] + k[11], k[8]*k[20] + k[8]*k[21] + k[8]*k[22] + k[9]*k[20] + k[9]*k[21] + k[9] + k[10]*k[22] + k[10]*k[23] + k[11]*k[21] + k[11]*k[23] + k[11], k[8]*k[21] + k[8]*k[22] + k[8]*k[23] + k[9]*k[20] + k[9]*k[21] + k[9]*k[23] + k[9] + k[10]*k[20] + k[10]*k[23] + k[11]*k[23] + k[11], k[8]*k[21] + k[8]*k[23] + k[8] + k[9]*k[21] + k[9]*k[22] + k[9]*k[23] + k[10]*k[20] + k[10]*k[22] + k[10] + k[11]*k[20] + k[11]*k[21] + k[11]*k[22] + k[20] + k[22] + k[23] + 1, k[8]*k[20] + k[8]*k[21] + k[8]*k[22] + k[9]*k[22] + k[9] + k[10]*k[20] + k[10]*k[21] + k[10]*k[23] + k[10] + k[11]*k[21] + k[11] + k[20] + k[21] + k[23] + 1, k[8]*k[20] + k[8]*k[22] + k[8] + k[9]*k[20] + k[9]*k[21] + k[9]*k[22] + k[10]*k[22] + k[10] + k[11]*k[20] + k[11]*k[21] + k[11]*k[23] + k[11] + k[20] + k[21] + k[22], k[8]*k[21] + k[8]*k[22] + k[8]*k[23] + k[9]*k[20] + k[9]*k[22] + k[9] + k[10]*k[20] + k[10]*k[21] + k[10]*k[22] + k[11]*k[22] + k[11] + k[21] + k[22] + k[23], k[0] + k[8] + k[16] + k[24], k[1] + k[9] + k[17] + k[25] + 1, k[2] + k[10] + k[18] + k[26], k[3] + k[11] + k[19] + k[27], k[4] + k[12] + k[20] + k[28], k[5] + k[13] + k[21] + k[29], k[6] + k[14] + k[22] + k[30], k[7] + k[15] + k[23] + k[31]]
    return e