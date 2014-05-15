#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>


#define INFORMATION_SPAN    21
#define REQUIRED_SAMPLES    (1<<INFORMATION_SPAN)
#define CODE_LENGTH         7
//#define INLINE

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
#ifdef INLINE
inline 
#endif
uint64_t in_ball_rand(int p) {
    int i, x;
    uint64_t res = 0x0ull;
    for(i = 0; i < p; ++i) {
        x = 63-INFORMATION_SPAN+(rand() % INFORMATION_SPAN);
        res ^= (0x1ull << x);
    }
    return res;
}

/* random function for 32/64-bit systems */
#ifdef INLINE
inline 
#endif
uint64_t rand64() {
    #ifdef OSX_64
        rand();
    #else
        uint64_t num;
        num = rand();
        num = (num << 32) | rand();
        return num;
    #endif
}

#ifdef INLINE
inline 
#endif
uint64_t parity(uint64_t v) {
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
        lookup[i] = i ^ ((result != 0) << (result-1));
    }
}

void transform_to_distr(uint64_t* from_queries, int32_t* to_distr) {
    int i;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        to_distr[from_queries[i] >> (63-INFORMATION_SPAN)] += (0x1ull & from_queries[i]) ? -1 : 1;
    }
}

/* Thanks to Paul */
void FWHT(int32_t *distr) {
    uint64_t i, j, k;
    for (i = 0; i < INFORMATION_SPAN; ++i) {  
        for (j = 0; j < (0x1ull << INFORMATION_SPAN); j += (0x1ull << (i+1))) {
            for (k = 0; k < (0x1ull << i); ++k) {
                int32_t a = distr[(uint32_t)(j+k)];
                int32_t b = distr[(uint32_t)(j+k+(0x1ull << i))];
                distr[(uint32_t)(j+k)] = a + b;
                distr[(uint32_t)(j+k+(0x1ull << i))] = a - b;
            }
        }
    }
}

uint64_t get_max_solution(int32_t *distr) {
    int32_t max_val = -100000, entry;
    int i;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        if (distr[i] > max_val) {
            max_val = distr[i];
            entry = i;
        }
    }
    return entry;
}

/* function for testing a hypothesis (counting incorrect) */
int test_hypothesis(uint64_t* ptr, uint64_t secret) {
    int i = 0, incorrect = 0;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        incorrect += (ptr[i]&0x1ull) != parity(ptr[i]&secret);
    }
    return incorrect;
}

/* this is our oracle function */
void generate_queries(uint64_t* ptr, uint64_t secret) {
    int i;
    uint64_t v, rnd_vect;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        rnd_vect = rand64() & mask;
        ptr[i] = parity(rnd_vect&secret); // ^ noise
        v = lookup[(rnd_vect >> (63-INFORMATION_SPAN)) & code_mask] << (63-INFORMATION_SPAN);
        v ^= lookup[(rnd_vect >> (63-INFORMATION_SPAN+CODE_LENGTH)) & code_mask] << (63-INFORMATION_SPAN+CODE_LENGTH);
        v ^= lookup[(rnd_vect >> (63-INFORMATION_SPAN+2*CODE_LENGTH)) & code_mask] << (63-INFORMATION_SPAN+2*CODE_LENGTH);
        ptr[i] ^= v;
    }
}

/* main, sweet main */
int main ( int arc, char **argv ) { 
    int pN = 7, discrepancies, i;
    
    /* set the random seed */
    srand(2009382);
    
    /*  define secret and query list */
    uint64_t secret = in_ball_rand(pN) & mask;
    uint64_t* queries = malloc(sizeof(uint64_t)*REQUIRED_SAMPLES);
    int32_t* distr = malloc(sizeof(int32_t)*REQUIRED_SAMPLES);
    
    /* make lookup table */
    make_table();
    
    /* generate the list of queries */
    generate_queries(queries, secret);
    transform_to_distr(queries, distr);
    
    printf("p = %f, bias € = %f\n\n", 1.0*3/INFORMATION_SPAN, pow(2.0*3/INFORMATION_SPAN-1,pN));
    
    discrepancies = test_hypothesis(queries, secret);
    printf("Observed[ #incorrect | correct hyp ] = ");
    printf("%d",discrepancies);
    printf(", bias € = %f\n", pow(2.0*discrepancies/REQUIRED_SAMPLES-1,1));
    
    //for(i=0;i<50;++i) {
     //   printf("%d, ",((queries[i]&0x1ull) != parity(queries[i]&secret)));
     //   print_bits(queries[i]);
    //}
    print_bits(secret);
    secret = (rand64()& mask);
/*    printf("   ");
    print_bits(secret);*/
    discrepancies = test_hypothesis(queries, secret);
    printf("Observed[ #incorrect | incorrect hyp ] = ");
    printf("%d of total %d",discrepancies, REQUIRED_SAMPLES);
    printf(", bias € = %f\n", 2.0*discrepancies/REQUIRED_SAMPLES-1);
    
    for(i = 0; i < (1<<7); ++i) {
        /* i is the received codeword */
    //    print_bits(lookup[i]);
    }
    
    FWHT(distr);
    
    secret = get_max_solution(distr) << (63-INFORMATION_SPAN);
    print_bits(secret);
    discrepancies = test_hypothesis(queries, secret);
    printf("Observed[ #incorrect | incorrect hyp ] = ");
    printf("%d of total %d",discrepancies, REQUIRED_SAMPLES);
    printf(", bias € = %f\n", 2.0*discrepancies/REQUIRED_SAMPLES-1);
    
    //print_bits(queries[0]);
    //print_bits(queries[1]);
    //print_bits(queries[2]);
    
    return 0; // Indicates that everything went well.
}