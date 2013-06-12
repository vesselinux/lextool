# 
# Original file: ks2-roundkeys-2x2x4-2-keys.sage
# 
# round keys according to DeCanniere-Biruykov arrangement of lex-2x2x4 for [2] keys 
# the key variables are input and output of sboxes 
# the input key x is expected to be in the following format 
# x = [[[x0, x1, x2, x3, x4, x5, x6, x7], [x8, x9, x10, x11, x12, x13, x14, x15]], [[x16, x17, x18, x19, x20, x21, x22, x23], [x24, x25, x26, x27, x28, x29, x30, x31]]]
# the output is the round keys added at the end of each round according to DCB arrangement
def ks2_round_keys(k):
    e = [[[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]], [k[8], k[9], k[10], k[11], k[12], k[13], k[14], k[15]]], [[k[0] + k[16], k[1] + k[17] + 1, k[2] + k[18], k[3] + k[19], k[4] + k[20], k[5] + k[21], k[6] + k[22], k[7] + k[23]], [k[24], k[25], k[26], k[27], k[28], k[29], k[30], k[31]]]]
    return e
