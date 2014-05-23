/*
LPN solver using covering codes
This code simulates the subspace hypothesis testing

report bugs to carl@eit.lth.se

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <time.h>

#define SECRET_WEIGHT       5
#define BKW_ITERATIONS      4
#define QUERY_LENGTH        60
#define REPETITIONS         3

#define REQUIRED_SAMPLES    (1 << 24)
#define INFORMATION_SPAN    (QUERY_LENGTH/REPETITIONS)
#define ERROR_RATIO         0.5*(1-(pow(0.8, (1 << BKW_ITERATIONS))))

#define MASK                ((0x1ull << QUERY_LENGTH)-1)
#define INFO_MASK           ((0x1ull << INFORMATION_SPAN)-1)
#define REP_MASK            ((0x1ull << REPETITIONS)-1)
#define REP_MASK_TABLE      ((0x1ull << (REPETITIONS*TABLE_FACTOR))-1)

#define TABLE_FACTOR        5
#define TRESHOLD            3
#define CONSTANT_FLAG       ((0x1ull<<62)-1)

uint64_t repetition_table[(1 << (REPETITIONS*TABLE_FACTOR))];
int total_discarded = 0;

/* print bits in decending magnitude */
void 
print_bits32(uint32_t word) {
    int i;
    printf("msb -> ");
    for(i = 31; i >= 0; --i) {
        printf("%d", (int)((word >> i) & 1));
    }
    printf("\n");
}

uint64_t 
decode_repetition(uint64_t codeword, int information_bits) {
	uint64_t y = 0x0ull; 
	int i, j, total = 0;
	for(i = 0; i < information_bits; i++) {
		j = codeword & REP_MASK;
    	j = (j & 0x55) + (j >> 1 & 0x55);
    	j = (j & 0x33) + (j >> 2 & 0x33);
    	j = (j & 0x0f) + (j >> 4 & 0x0f);
		if (j > 1) {
			y ^= (j > 1) << (i);
		}
		total += (j == 2) || ( j== 2);
		codeword >>= REPETITIONS;
	}
	/* for simplicity, we check the weight here
	 	even though we may discard some more samples */
	if (total > TRESHOLD) {
		y = CONSTANT_FLAG;
	}
	return y;
}


void 
make_table() {
	uint64_t i;
	for(i = 0; i < (1 << (REPETITIONS*TABLE_FACTOR)); ++i) {
		repetition_table[i] = decode_repetition(i, (REPETITIONS*TABLE_FACTOR));
	}
}

uint64_t 
decode_repetition_table(uint64_t codeword, int information_chunks) {
	uint64_t y = 0x0ull, j = 0x0ull; 
	int i;
	for(i = 0; i < information_chunks; i++) {
		j = repetition_table[codeword & REP_MASK_TABLE] << (i*TABLE_FACTOR);
		if (repetition_table[codeword & REP_MASK_TABLE] == CONSTANT_FLAG) {
			return CONSTANT_FLAG;
		}
		y ^= j;
		codeword >>= (REPETITIONS*TABLE_FACTOR);
	}
	return y;
}

/* generate a vector from a hammingball of size p */
uint64_t in_ball_rand(int p) {
    int i, x;
    uint64_t res = 0x0ull;
    for(i = 0; i < p; ++i) {
        x = rand() % QUERY_LENGTH;
        res ^= (0x1ull << x);
    }
    return res;
}

/* random function for 32/64-bit systems */
uint64_t rand64() {
    #ifdef OSX_64
        return rand();
    #else
        uint64_t num;
        num = rand();
        num = (num << 32) | rand();
        return num;
    #endif
}

/* parity function */
#ifdef INLINE
inline 
#endif
uint64_t parity(uint64_t v) {
    /* calculate parity */
    v ^= v >> 1;
    v ^= v >> 2;
    v = (v & 0x1111111111111111ul) * 0x1111111111111111ul;
    return (v >> 60) & 0x1ull;
}

/* make table into frequency distribution table */
void 
transform_to_distr(uint64_t* from_queries, int32_t* to_distr) {
    int i;
    for(i = 0; i < (1 << INFORMATION_SPAN); ++i) {
        to_distr[i] = 0;
    }
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        to_distr[from_queries[i] & INFO_MASK] += (0x1ull & (from_queries[i] >> 63)) ? -1 : 1;
    }
}

/* fast walsh-hadamard transform */
void 
FWHT(int32_t *distr) {
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

uint64_t 
transform_solution(uint64_t secret) {
	int i;
	uint64_t answer = 0x0ull;
	for (i = 0; i < INFORMATION_SPAN; ++i) { 
		answer ^= (parity(secret & REP_MASK) << i);
		secret = secret >> REPETITIONS;
	}
	return answer;
}

/* find-max */
uint64_t 
get_max_solution(int32_t *distr) {
    int32_t max_val = -100000, entry;
    int i;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
        if (abs(distr[i]) > max_val) {
            max_val = abs(distr[i]);
            entry = i;
        }
    }
    return entry;
}

/* function for testing a hypothesis (counting incorrect) */
int 
test_hypothesis(uint64_t* ptr, uint64_t secret) {
    int i = 0, incorrect = 0;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
		//print_bits((ptr[i] >> 63));
        incorrect += (ptr[i] >> 63) != parity(ptr[i]&secret);
    }
    return incorrect;
}

/* this is our LPN-oracle function */
void 
generate_queries(uint64_t* ptr, uint64_t secret) {
    int i;
    uint64_t rnd_vect, noise;
    for(i = 0; i < REQUIRED_SAMPLES; ++i) {
		noise = (rand()/(double)RAND_MAX) < ERROR_RATIO;
        rnd_vect = rand64() & MASK;
        ptr[i] = decode_repetition_table(rnd_vect, INFORMATION_SPAN/TABLE_FACTOR);
		if (ptr[i] == CONSTANT_FLAG) {
			i--;
			total_discarded++;
		} else {
		ptr[i] ^= (parity(rnd_vect&secret) ^ noise) << 63;
		}
    }
}

int 
compare(const int32_t *a, const int32_t *b) {
    return  (abs(*a) < abs(*b))-(abs(*a) > abs(*b));
}

/* main, sweet main */
int 
main ( int arc, char **argv ) {
	int count = 0; 
	int32_t x;
    /* set the random seed */
	srand(0xCAFE);
    /*  define secret, query list and distribution */
    uint64_t secret = in_ball_rand(SECRET_WEIGHT);
    uint64_t* queries = malloc(sizeof(uint64_t)*REQUIRED_SAMPLES);
    int32_t* distr = malloc(sizeof(int32_t)*REQUIRED_SAMPLES);
	/* make a decoding table */
	make_table();
    /* generate the list of queries */
    generate_queries(queries, secret);
    transform_to_distr(queries, distr);
	secret = transform_solution(secret);
	/* perform hypothesis testing using fast walsh-hadamard transform */
	FWHT(distr);
	/* determine bias of the secret and the best bias solution*/
	printf("\nSample bias:\nTheoretical    \t%f \n",
		(1-2.0*ERROR_RATIO)*pow(1.0-2.0*TRESHOLD/REPETITIONS/TABLE_FACTOR, SECRET_WEIGHT));
	printf("Correct secret \t%f \n",1-2.0*test_hypothesis(queries, secret)/REQUIRED_SAMPLES);
    printf("Found secret   \t%f\n\n", 1-2.0*test_hypothesis(queries, get_max_solution(distr))/REQUIRED_SAMPLES);
	/* find first occurrence of the bias */
	x = distr[secret];
	qsort(distr, REQUIRED_SAMPLES, sizeof(int32_t), compare);
	while(abs(distr[count]) > abs(x)) count++;
	printf("The correct solution is at position %d\n", count);
	printf("Total discarded samples: 2^%f\n", log(total_discarded)/log(2));
	/* clean up */
    free(queries);
    free(distr);
    return 0; 
}