#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define REQUIRED_SAMPLES    (1<<21)
#define INFORMATION_SPAN    21
#define CODE_LENGTH         7

//#define OSX_64

uint64_t mask = (~0ull)<<(63-INFORMATION_SPAN);
uint64_t code_mask = 0x7full;
uint64_t lookup[1<<7];

/* print bits in decending magnitude */
void print_bits(uint64_t dword) {
    int i;
    printf("msb -> ");
    for(i = 63; i >= 0; --i) {
        printf("%d", (int)((dword >> i) & 1));
    }
    printf("\n");
}


/* generate a vector from a hammingball of size p */
inline uint64_t in_ball_rand(int p) {
    int i, x;
    uint64_t res = 0x0ull;
    for(i = 0; i < p; ++i) {
        x = 63-INFORMATION_SPAN+(rand() % INFORMATION_SPAN);
        res ^= (0x1ull << x);
    }
    return res;
}

/* random function for 32/64-bit systems */
inline uint64_t rand64() {
    
    #ifdef OSX_64
        rand();
    #else
        uint64_t num;
        num = rand();
        num = (num << 32) | rand();
        return num;
    #endif

}

inline uint64_t parity(uint64_t v) {
    /* calculate parity */
    v ^= v >> 1;
    v ^= v >> 2;
    v = (v & 0x1111111111111111ul) * 0x1111111111111111ul;
    return (v >> 60) & 1;
}



/* make a table of the code */
void make_table() {
    
    /* we use a [7,4,3] hamming code 
    
    0001111
    0110011  = H
    1010101
    
    */
    int x[8] = {0,0,0,0,0,0,0,0};
    uint64_t code[3] = {0xfull,
                        0x33ull,
                        0x55ull};
    
    uint64_t i, result = 0, row;
    for(i = 0; i < (1<<7); ++i) {
        /* i is the received codeword */
        result = 0;
        for (row = 0; row < 3; ++row) {
            result ^= parity(code[row]&i) << row;
        }
        if (result != 0)
            lookup[i] = i;// ^ (1 << (result-1));
    }
    
}

/* function for testing a hypothesis (counting incorrect) */
int test_hypothesis(uint64_t* ptr, uint64_t secret) {
    int i = 0, incorrect = 0;
    uint64_t v;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        incorrect += (ptr[i] & 1) != parity(ptr[i]&secret);
    }
    return incorrect;
}

/* this is our oracle function */
void generate_queries(uint64_t* ptr, uint64_t secret) {
    
    int i = 0;
    uint64_t v, rnd_vect;
    for(i = 0; i < 2; ++i) {
        rnd_vect = rand64() & mask;
        ptr[i] = parity(rnd_vect&secret); // ^ noise
        v = lookup[(rnd_vect >> (63-INFORMATION_SPAN)) & code_mask] << (63-INFORMATION_SPAN);
        v ^= lookup[(rnd_vect >> (63-INFORMATION_SPAN+CODE_LENGTH)) & code_mask] << (63-INFORMATION_SPAN+CODE_LENGTH);
        v ^= lookup[(rnd_vect >> (63-INFORMATION_SPAN+2*CODE_LENGTH)) & code_mask] << (63-INFORMATION_SPAN+2*CODE_LENGTH);
        ptr[i] ^= v;
        
        print_bits(rnd_vect);
        print_bits((rnd_vect >> (63-INFORMATION_SPAN+7)));
        print_bits((rnd_vect >> (63-INFORMATION_SPAN+7))&code_mask);
        print_bits((rnd_vect >> (63-INFORMATION_SPAN+7))&code_mask);
        print_bits(((rnd_vect >> (63-INFORMATION_SPAN+7)) & code_mask) << (63-INFORMATION_SPAN+7));
        print_bits(v^rnd_vect);
        printf("\n");
    }
}

/* main, sweet main */
int main ( int arc, char **argv ) {
    
    /* set the random seed */
    srand(2009712);
    
    /*  define secret and query list */
    uint64_t secret = in_ball_rand(1) & mask;
    uint64_t* queries = malloc(sizeof(uint64_t)*REQUIRED_SAMPLES);
    
    /* make lookup table */
    make_table();
    
    /* generate the list of queries */
    generate_queries(queries, secret);
    
    printf("Observed[#incorrect | correct hyp] = ");
    printf("%d",test_hypothesis(queries, secret));
    printf("\n");
    printf("Observed[#incorrect | incorrect hyp] = ");
    printf("%d",test_hypothesis(queries, secret+1));
    
    return 0; // Indicates that everything went well.
}