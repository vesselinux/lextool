Lex Toolkit Usage Instructions

NOTE: Tested with Sage Version 4.1.1, Release Date: 2009-08-14

The Lex Toolkit is a collection of Python programs for the computer algebra system Sage. The programs generate Boolean algebraic equations for a small-scale version of stream cipher LEX. The small-scale version of Lex, called Lex(2,2,4) is based on the small-scale version of AES-128 - the block cipher SR(10,2,2,4) proposed by Murphy an Robshaw. Lex(2,2,4) uses 4-bit words and has a state of 4 words (16 bits). At every round Lex(2,2,4) leaks 4 bits of the state i.e. one-fourth of the whole state as the full-scale version of Lex. The Lex Toolkit was used for the algebraic cryptanalysis of small-scale LEX. For more detials refer to: 

[1] V. Velichkov, V. Rijmen, and B. Preneel, "Algebraic Cryptanalysis of a Small-Scale Version of Stream Cipher LEX," IET Information Security 4(49), 16 pages, 2010. 

The toolkit is composed of three programs: lex224.sage, lex224-quadratic-main.sage and lex224-cubic-main.sage. Each of them can be started by executing the corresponding command in Sage from within the directory of the toolkit:

sage: load lex224.sage
sage: load quadratic/lex224-quadratic-main.sage
sage: load cubic/lex224-cubic-main.sage

Brief description of the three programs follows.

- The program 'lex224.sage' implements several small-scale variants of Lex (including Lex(2,2,4)), based on small-scale variants of AES proposed in Murphy,Robshaw paper "Small Scale Variants of the AES".

- The program 'lex224-quadratic-main.sage' generates quadratic equations for Lex(2,2,4) and solves them using Groebner bases to recover the secret key of the cipher. The default values of the key, the plaintext, the number of rounds, and the size of the leak are set to the example presented in the Appendix of the paper [1]. Those values can be easily modified within the code. 

- The program 'lex224-cubic-main.sage' generate quadratic equations for Lex(2,2,4) and solves them using Groebner bases to recover the secret key of the cipher.

Below is provided a description of the files present in the toolkit directory tree.

./

lex224.sage                       : Several small-scale variants of Lex (including Lex(2,2,4)), based on small-scale variants of AES proposed in Murphy,Robshaw paper "Small Scale Variants of the AES".

./quadratic

lex224-quadratic-main.sage        : Main file 1: Generate quadratic equations for Lex(2,2,4). Uses the following files:
|
L linear-eqs.sage                 : Generate equations for the linear part (ShiftRows + MixColumns + AddRoundKkey) of one round of Lex(2,2,4)
L linear-eqs-precomputed.sage     : Pre-computed equations using linear-eqs.sage
L sbox-eqs.sage                   : Generate 11 quadratic equations for the SR(10,2,2,4) sbox in GF(2^4)
L sbox-eqs-precomputed.sage       : Pre-computed 11 quadratic equations using sbox-eqs.sage
L keysched-eqs.sage               : Generate equations for the key schedule of Lex(2,2,4) for 2 keys
L keysched-eqs-precomputed.sage   : Pre-computed equations using keysched-eqs.sage
L roundkeys-dcb.sage              : Ordering of the round key bits according to De Canniere, Biruykov paper "Block ciphers and systems of quadratic equations"
L leaked-bits-precomputed.sage    : Precomputed leaks for Lex(2,2,4) for a given key. Computed using leaked-bits.sage.
L leaked-bits.sage                : Compute the full output from several rounds of Lex(2,2,4) for a given master key (16 bits) and plaintext (16 bits)

./cubic

lex224-cubic-main.sage            : Main file 2: Generate quadratic equations for Lex(2,2,4). Uses the following files:
|
L round-eqs-precomputed.sage      : Precomputed cipher equations for 1 round, computed using lex224.sage
L keysched-eqs.sage               : Generate key schedule equations for 1 round
L keysched-eqs-precomputed.sage   : Precomputed key schedule equations using keysched-eqs.sage
L leaked-bits-precomputed.sage    : Generate leaks for 10 rounds. Computed using leaked-bits.sage.
L leaked-bits.sage                : Compute the full output from several rounds of Lex224 for a given master key (16 bits) and plaintext (16 bits)


Although all files provided with the tool were renamed, in the comments within the code the old names are sometimes used. Because of that below we provide a table with the new names and the corresponding old names.

OLD NAME                          : NEW NAME

lexs.sage                         : lex224.sage

lexs-compute-2x2x4.sage           : lex224-quadratic-main.sage
|
L lin-compute-2x2x4.sage          : linear-eqs.sage
L lin-eqs-2x2x4.sage              : linear-eqs-precomputed.sage
L sbox-compute-2x2x4.sage         : sbox-eqs.sage
L sbox-eqs-2x2x4.sage             : sbox-eqs-precomputed.sage
L ks2-compute-2x2x4.sage          : keysched-eqs.sage
L ks2-eqs-2x2x4-2-keys.sage       : keysched-eqs-precomputed.sage
L ks2-roundkeys-2x2x4-2-keys.sage : roundkeys-dcb.sage
L leaks-2x2x4.sage                : leaked-bits-precomputed.sage

lexs-2x2x4.sage                   : lex224-cubic-main.sage
|
L eqss-2x2x4.sage                 : round-eqs-precomputed.sage
L ks-compute-2x2x4.sage           : keysched-eqs.sage
L ks-2x2x2x4.sage                 : keysched-eqs-precomputed.sage
L leaks-2x2x4.sage                : leaked-bits-precomputed.sage
