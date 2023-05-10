/*
===========================================================================================================================
* AUTHOR:     Christopher Riccardi, PhD student at Computational Biology Dept., University of Firenze
* PROJECT:    Fast ungapped mapping of sequencing reads using Squid; Feb 2022
* CONTACT:    CHRISTOPHER.RICCARDI@UNIFI.IT
===========================================================================================================================
*/

#pragma GCC diagnostic ignored "-Wunused-variable"
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <pthread.h>

#ifdef __APPLE__
#include <limits.h>
#endif // !__APPLE__
int VERBOSE = 1;
int DIFF = 0;
int MASK_LOWER = 0;
int AVAIL_THREADS = 1;
unsigned short MISMATCH = 15;
unsigned short K = 11;
int STEP = 17;
int FASTQ_OUT = 1;
int BED_OUT = 1;
int NO_DISJOIN = 1;
int IGNORE_N = 0;
int EVALS = 0;
#define REQUIRE_ARG_ERR(message, type) {                            \
    fprintf(stderr,                                                 \
        "[ERROR] %s option requires an argument of type %s.\n",     \
        message, type);                                             \
}
#define MEMORY_FAILURE(pointer) {                                   \
        if (pointer==NULL) {                                        \
            fprintf (stderr, "%s:%d: "                              \
                     " memory allocation failed. (Error: %p)\n",    \
                     __FILE__, __LINE__,                            \
                     strerror (errno));                             \
            exit (EXIT_FAILURE);                                    \
        }                                                           \
    }
#define IO_FAILURE(bytes) {                                         \
        if (bytes==0) {                                             \
            fprintf (stderr, "\n%s:%d: "                            \
                     " I/O failure on file. (Error: %p)\n",         \
                     __FILE__, __LINE__,                            \
                     strerror (errno));                             \
            exit (EXIT_FAILURE);                                    \
        }                                                           \
    }
#define VERBOSE_MSG(type, message) {                                \
    if(VERBOSE){                                                    \
        if(type=='W'){                                              \
            fprintf(stderr, "[Warning] %s\n", message);             \
        }                                                           \
        else if(type=='E'){                                         \
            fprintf(stderr, "[Error] %s\n", message);               \
            exit(EXIT_FAILURE);                                     \
        }                                                           \
        else if(type=='L'){                                         \
            fprintf(stderr, "%s", message);                         \
        }                                                           \
    }                                                               \
}
#define INDEPENDENT_MSG(type, message) {                            \
      if(type=='W'){                                                \
          fprintf(stderr, "[Warning] %s\n", message);               \
      }                                                             \
      else if(type=='E'){                                           \
          fprintf(stderr, "[Error] %s\n", message);                 \
          exit(EXIT_FAILURE);                                       \
      }                                                             \
      else if(type=='L'){                                           \
          fprintf(stderr, "%s", message);                           \
    }                                                               \
}
#define VERBOSE_STR(message, str) {                                 \
    if(VERBOSE) fprintf(stderr, "%s %s\n", message, str);           \
}
#define READ_SIZE 512
#define BUF_SIZE 8192   // buffer for reading a line in fgets/gzgets
#define PATH 256        // max size for a file path

static const char LIB_TYPE[9][4] = { "ISF",  "ISR",  "IU",
				  "OSF",  "OSR",  "OU",
				  "SF",   "SR",   "U" };

/* Structures */
typedef struct {
	char header[READ_SIZE];
	char sequence[READ_SIZE];
	char placeholder[READ_SIZE];
	char quality[READ_SIZE];
} READ;

typedef struct {
	char* chrom;            // The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)
	uint32_t chromStart;    // The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
	uint32_t chromEnd;      // The ending position of the feature in the chromosome or scaffold.
	char* name;             // Defines the name of the BED line
} BED;

typedef struct {
	char* chrom1;       // The name of the chromosome on which the first end of the feature exists
	uint32_t start1;    // The zero-based starting position of the first end of the feature on chrom1
	uint32_t end1;      // The one-based ending position of the first end of the feature on chrom1
	char* chrom2;       // The name of the chromosome on which the second end of the feature exists
	uint32_t start2;    // The zero-based starting position of the second end of the feature on chrom2
	uint32_t end2;      // The one-based ending position of the second end of the feature on chrom2
	char* name;         // Defines the name of the BEDPE feature
	int score;          // The UCSC definition requires that a BED score range from 0 to 1000, inclusive
	char strand1;       // Defines the strand for the first end of the feature. Either ‘+’ or ‘-‘
	char strand2;       // Defines the strand for the second end of the feature. Either ‘+’ or ‘-‘
} BEDPE;

typedef struct {
	char* db;
	char* input_R1;
	char* input_R2;
	char output_R1[PATH];
	char output_R2[PATH];
	char output_BED[PATH];
	char basename[PATH];
	unsigned short quiet;
	short lib;
} PARAMS;

typedef struct {
	char header[1024];
	size_t size;
	double GC_content;
	char* sequence;
} FASTA;

typedef struct {
	uint32_t id;
	uint32_t size;
	uint32_t** pos;
} HASH;

typedef struct {
	HASH* hash;
	uint32_t size;
} HASH_TBL;

typedef struct {
	FASTA* fasta;
	uint32_t size;
} INDEX;

typedef struct {
	PARAMS* params;
	INDEX* index;
	HASH_TBL* hash_tbl;
	long bytes_R1_start;
	long bytes_R2_start;
	size_t lines;
	char temp_out_R1[PATH];
	char temp_out_R2[PATH];
	char temp_out_BED[PATH];
	uint32_t hashfunc;
} THREAD;


/* General Functions */
void print_usage() {
	fprintf(stderr,
		"Thank you for using Squid\n"
		"Usage: squid -i <str> -R1 <str> [-R2 <str>] -o <str> -l <str> [Options]\n\n"
		"Please type \"squid -h\" to see a detailed help menu\n");
	exit(EXIT_FAILURE);
}

void print_help() {
	fprintf(stderr, "\n"
		"Thank you for using Squid\n\n"
		"Usage: squid -i <str> -R1 <str> [-R2 <str>] -l <str> -o <str> [Options]\n\n"
		"Mandatory arguments:\n"
		"   -i         <str>         input database in FASTA format (can be gzipp'd)\n"
		"   -R1        <str>         read in forward direction (R1) (can be gzipp'd)\n"
		"   -R2        <str>         read in reverse direction (R2) (can be gzipp'd)\n"
		"   -o         <str>         basename, Squid will add \"_R1.fastq\", \"_R2.fastq\" and/or \".bed\"\n"
		"\n   \"-l\" argument is also mandatory, with one of the following format strings:\n"
		"   -l SF                    Stranded Forward. R1 comes from the forward strand, R2 from the reverse strand\n"
		"   -l SR                    Stranded Reverse. R1 comes from the reverse strand, R2 from the forward strand\n"
		"   -l U                     Unstranded. R1 or R2 can derive from both strands.\n"
		"   -l ISF                   Inward Stranded Forward. R1 and R2 behave as in SF. R1 must map upstream to R2\n"
		"   -l ISR                   Inward Stranded Reverse. R1 and R2 behave as in SR. R1 must map downstream to R2\n"
		"   -l IU                    Inward Unstranded. R1 and R2 behave as in U. With this option Squid tries ISF and ISR\n"
		"   -l OSF                   Outward Stranded Forward. R1 and R2 behave as in SF. R1 must map downstream to R2\n"
		"   -l OSR                   Outard Stranded Reverse. R1 and R2 behave as in SR. R1 must map upstream to R2\n"
		"   -l OU                    Outard Unstranded. R1 and R2 behave as in U. With this option Squid tries OSF and OSR\n"
		"\nSquid also provides a number of additional arguments for a more flexible mapping.\n\n"
		"Boolean arguments:\n"
		"   --diff                   when FASTQ(s) output is enabled, return reads that do not map to database.\n"
		"                            By default this is switched off, meaning that only mapping reads will be written.\n"
		"   --disjoin                when database is a multi-FASTA, allow R1 and R2 to map to different sequences.\n"
		"                            Default is to coerce R1 and R2 to map to the the same seqid. When on,\n"
		"                            disjoined read pairs will switch the score field in the BEDPE to 1 instead of 0.\n"
		"   --ignore_N               do not treat Ns as mismatches, simply ignore them (default: OFF)\n"		
		"   --mask-lower             do not capitalize lowercase letters in database (default is to make them uppercase)\n"

		"   --no-bed                 do not produce BED/BEDPE output file (default is to write it)\n"
		"   --no-fastq               do not produce FASTQ output file(s) (default is to write them)\n"
		"   --quiet                  do not print log to stderr (default is to be verbose)\n\n"
		"Scanning and performance arguments:\n"
		"   -e         <int>         evaluate <int> number of alternative positionings of R1 and R2, looking for a better match.\n"
		"                            Default is to break as soon as a suitable match is found (-e 0). This option is meaningful\n"
		"                            when BED/BEDPE output is enabled. Greater values of <int> affect performance but could\n"
		"                            report a more accurate mapping when higly similar sequences are present in the database.\n"
		"   -k         <int>         kmer size: 9, 11, 13 or 15 (default: 11)\n"
		"   -m         <int>         max %% of mismatches allowed during ungapped extension\n"
		"                            Default is to force 85%% sequence identity, hence -m 15.\n"
		"   -s         <int>         step size while sliding over the sequencing reads\n"
		"                            for a perfect match of length k.\n"
		"                            Lower s increases sensitivity but decreases speed.\n"
		"                            Min=1 (sliding window of 1), default: 17.\n"
		"   -t         <int>         number of threads (default: 1)\n\n"
	);
	exit(EXIT_FAILURE);
}

int is_int(char* string) {
	if (strlen(string) == 0) return 0;
	if (strlen(string) == 1 && !isdigit(*string)) return 0;
	if (strlen(string) > 1 && string[0] == '-') {
		++string;
		while (*string != 0)
		{
			if (!isdigit(*string++))  return 0;
		}
	}
	else if (strlen(string) > 1 && string[0] != '-') {
		while (*string != 0) {
			if (!isdigit(*string++)) return 0;
		}
	}
	return 1;
}


/* Pseudo-Hash Functions */
uint32_t calc_hash_9(const char* s) {
	uint32_t arr[9] = { 0 };
	switch (s[8])
	{
	case 'A':
	{
		arr[7] = 4;
	}break;

	case 'C':
	{
		arr[8] = 1;
	}break;
	case 'G':
	{
		arr[8] = 2;
	}break;
	case 'T':
	{
		arr[8] = 3;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[7])
	{
	case 'A':
	{
		arr[7] = arr[8];
	}break;
	case 'C':
	{
		arr[7] = arr[8] + 4;
	}break;
	case 'G':
	{
		arr[7] = arr[8] + 8;
	}break;
	case 'T':
	{
		arr[7] = arr[8] + 12;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[6])
	{
	case 'A':
	{
		arr[6] = arr[7];
	}break;
	case 'C':
	{
		arr[6] = arr[7] + 16;
	}break;
	case 'G':
	{
		arr[6] = arr[7] + 32;
	}break;
	case 'T':
	{
		arr[6] = arr[7] + 48;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[5])
	{
	case 'A':
	{
		arr[5] = arr[6];
	}break;
	case 'C':
	{
		arr[5] = arr[6] + 64;
	}break;
	case 'G':
	{
		arr[5] = arr[6] + 128;
	}break;
	case 'T':
	{
		arr[5] = arr[6] + 192;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[4])
	{
	case 'A':
	{
		arr[4] = arr[5];
	}break;
	case 'C':
	{
		arr[4] = arr[5] + 256;
	}break;
	case 'G':
	{
		arr[4] = arr[5] + 512;
	}break;
	case 'T':
	{
		arr[4] = arr[5] + 768;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[3])
	{
	case 'A':
	{
		arr[3] = arr[4];
	}break;
	case 'C':
	{
		arr[3] = arr[4] + 1024;
	}break;
	case 'G':
	{
		arr[3] = arr[4] + 2048;
	}break;
	case 'T':
	{
		arr[3] = arr[4] + 3072;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[2])
	{
	case 'A':
	{
		arr[2] = arr[3];
	}break;
	case 'C':
	{
		arr[2] = arr[3] + 4096;
	}break;
	case 'G':
	{
		arr[2] = arr[3] + 8192;
	}break;
	case 'T':
	{
		arr[2] = arr[3] + 12288;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[1])
	{
	case 'A':
	{
		arr[1] = arr[2];
	}break;
	case 'C':
	{
		arr[1] = arr[2] + 16384;
	}break;
	case 'G':
	{
		arr[1] = arr[2] + 32768;
	}break;
	case 'T':
	{
		arr[1] = arr[2] + 49152;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[0])
	{
	case 'A':
	{
		arr[0] = arr[1];
	}break;
	case 'C':
	{
		arr[0] = arr[1] + 65536;
	}break;
	case 'G':
	{
		arr[0] = arr[1] + 131072;
	}break;
	case 'T':
	{
		arr[0] = arr[1] + 196608;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	return arr[0];
}

uint32_t calc_hash_11(const char* s) {
	uint32_t arr[11] = { 0 };
	switch (s[10])
	{
	case 'A':
	{
		arr[9] = 4;
	}break;

	case 'C':
	{
		arr[10] = 1;
	}break;
	case 'G':
	{
		arr[10] = 2;
	}break;
	case 'T':
	{
		arr[10] = 3;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[9])
	{
	case 'A':
	{
		arr[9] = arr[10];
	}break;
	case 'C':
	{
		arr[9] = arr[10] + 4;
	}break;
	case 'G':
	{
		arr[9] = arr[10] + 8;
	}break;
	case 'T':
	{
		arr[9] = arr[10] + 12;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[8])
	{
	case 'A':
	{
		arr[8] = arr[9];
	}break;
	case 'C':
	{
		arr[8] = arr[9] + 16;
	}break;
	case 'G':
	{
		arr[8] = arr[9] + 32;
	}break;
	case 'T':
	{
		arr[8] = arr[9] + 48;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[7])
	{
	case 'A':
	{
		arr[7] = arr[8];
	}break;
	case 'C':
	{
		arr[7] = arr[8] + 64;
	}break;
	case 'G':
	{
		arr[7] = arr[8] + 128;
	}break;
	case 'T':
	{
		arr[7] = arr[8] + 192;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[6])
	{
	case 'A':
	{
		arr[6] = arr[7];
	}break;
	case 'C':
	{
		arr[6] = arr[7] + 256;
	}break;
	case 'G':
	{
		arr[6] = arr[7] + 512;
	}break;
	case 'T':
	{
		arr[6] = arr[7] + 768;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[5])
	{
	case 'A':
	{
		arr[5] = arr[6];
	}break;
	case 'C':
	{
		arr[5] = arr[6] + 1024;
	}break;
	case 'G':
	{
		arr[5] = arr[6] + 2048;
	}break;
	case 'T':
	{
		arr[5] = arr[6] + 3072;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[4])
	{
	case 'A':
	{
		arr[4] = arr[5];
	}break;
	case 'C':
	{
		arr[4] = arr[5] + 4096;
	}break;
	case 'G':
	{
		arr[4] = arr[5] + 8192;
	}break;
	case 'T':
	{
		arr[4] = arr[5] + 12288;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[3])
	{
	case 'A':
	{
		arr[3] = arr[4];
	}break;
	case 'C':
	{
		arr[3] = arr[4] + 16384;
	}break;
	case 'G':
	{
		arr[3] = arr[4] + 32768;
	}break;
	case 'T':
	{
		arr[3] = arr[4] + 49152;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[2])
	{
	case 'A':
	{
		arr[2] = arr[3];
	}break;
	case 'C':
	{
		arr[2] = arr[3] + 65536;
	}break;
	case 'G':
	{
		arr[2] = arr[3] + 131072;
	}break;
	case 'T':
	{
		arr[2] = arr[3] + 196608;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[1])
	{
	case 'A':
	{
		arr[1] = arr[2];
	}break;
	case 'C':
	{
		arr[1] = arr[2] + 262144;
	}break;
	case 'G':
	{
		arr[1] = arr[2] + 524288;
	}break;
	case 'T':
	{
		arr[1] = arr[2] + 786432;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[0])
	{
	case 'A':
	{
		arr[0] = arr[1];
	}break;
	case 'C':
	{
		arr[0] = arr[1] + 1048576;
	}break;
	case 'G':
	{
		arr[0] = arr[1] + 2097152;
	}break;
	case 'T':
	{
		arr[0] = arr[1] + 3145728;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	return arr[0];
}

uint32_t calc_hash_13(const char* s) {
	uint32_t arr[13] = { 0 };
	switch (s[12])
	{
	case 'A':
	{
		arr[11] = 4;
	}break;

	case 'C':
	{
		arr[12] = 1;
	}break;
	case 'G':
	{
		arr[12] = 2;
	}break;
	case 'T':
	{
		arr[12] = 3;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[11])
	{
	case 'A':
	{
		arr[11] = arr[12];
	}break;
	case 'C':
	{
		arr[11] = arr[12] + 4;
	}break;
	case 'G':
	{
		arr[11] = arr[12] + 8;
	}break;
	case 'T':
	{
		arr[11] = arr[12] + 12;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[10])
	{
	case 'A':
	{
		arr[10] = arr[11];
	}break;
	case 'C':
	{
		arr[10] = arr[11] + 16;
	}break;
	case 'G':
	{
		arr[10] = arr[11] + 32;
	}break;
	case 'T':
	{
		arr[10] = arr[11] + 48;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[9])
	{
	case 'A':
	{
		arr[9] = arr[10];
	}break;
	case 'C':
	{
		arr[9] = arr[10] + 64;
	}break;
	case 'G':
	{
		arr[9] = arr[10] + 128;
	}break;
	case 'T':
	{
		arr[9] = arr[10] + 192;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[8])
	{
	case 'A':
	{
		arr[8] = arr[9];
	}break;
	case 'C':
	{
		arr[8] = arr[9] + 256;
	}break;
	case 'G':
	{
		arr[8] = arr[9] + 512;
	}break;
	case 'T':
	{
		arr[8] = arr[9] + 768;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[7])
	{
	case 'A':
	{
		arr[7] = arr[8];
	}break;
	case 'C':
	{
		arr[7] = arr[8] + 1024;
	}break;
	case 'G':
	{
		arr[7] = arr[8] + 2048;
	}break;
	case 'T':
	{
		arr[7] = arr[8] + 3072;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[6])
	{
	case 'A':
	{
		arr[6] = arr[7];
	}break;
	case 'C':
	{
		arr[6] = arr[7] + 4096;
	}break;
	case 'G':
	{
		arr[6] = arr[7] + 8192;
	}break;
	case 'T':
	{
		arr[6] = arr[7] + 12288;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[5])
	{
	case 'A':
	{
		arr[5] = arr[6];
	}break;
	case 'C':
	{
		arr[5] = arr[6] + 16384;
	}break;
	case 'G':
	{
		arr[5] = arr[6] + 32768;
	}break;
	case 'T':
	{
		arr[5] = arr[6] + 49152;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[4])
	{
	case 'A':
	{
		arr[4] = arr[5];
	}break;
	case 'C':
	{
		arr[4] = arr[5] + 65536;
	}break;
	case 'G':
	{
		arr[4] = arr[5] + 131072;
	}break;
	case 'T':
	{
		arr[4] = arr[5] + 196608;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[3])
	{
	case 'A':
	{
		arr[3] = arr[4];
	}break;
	case 'C':
	{
		arr[3] = arr[4] + 262144;
	}break;
	case 'G':
	{
		arr[3] = arr[4] + 524288;
	}break;
	case 'T':
	{
		arr[3] = arr[4] + 786432;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[2])
	{
	case 'A':
	{
		arr[2] = arr[3];
	}break;
	case 'C':
	{
		arr[2] = arr[3] + 1048576;
	}break;
	case 'G':
	{
		arr[2] = arr[3] + 2097152;
	}break;
	case 'T':
	{
		arr[2] = arr[3] + 3145728;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[1])
	{
	case 'A':
	{
		arr[1] = arr[2];
	}break;
	case 'C':
	{
		arr[1] = arr[2] + 4194304;
	}break;
	case 'G':
	{
		arr[1] = arr[2] + 8388608;
	}break;
	case 'T':
	{
		arr[1] = arr[2] + 12582912;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[0])
	{
	case 'A':
	{
		arr[0] = arr[1];
	}break;
	case 'C':
	{
		arr[0] = arr[1] + 16777216;
	}break;
	case 'G':
	{
		arr[0] = arr[1] + 33554432;
	}break;
	case 'T':
	{
		arr[0] = arr[1] + 50331648;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	return arr[0];
}

uint32_t calc_hash_15(const char* s) {
	uint32_t arr[15] = { 0 };
	switch (s[14])
	{
	case 'A':
	{
		arr[13] = 4;
	}break;

	case 'C':
	{
		arr[14] = 1;
	}break;
	case 'G':
	{
		arr[14] = 2;
	}break;
	case 'T':
	{
		arr[14] = 3;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[13])
	{
	case 'A':
	{
		arr[13] = arr[14];
	}break;
	case 'C':
	{
		arr[13] = arr[14] + 4;
	}break;
	case 'G':
	{
		arr[13] = arr[14] + 8;
	}break;
	case 'T':
	{
		arr[13] = arr[14] + 12;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[12])
	{
	case 'A':
	{
		arr[12] = arr[13];
	}break;
	case 'C':
	{
		arr[12] = arr[13] + 16;
	}break;
	case 'G':
	{
		arr[12] = arr[13] + 32;
	}break;
	case 'T':
	{
		arr[12] = arr[13] + 48;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[11])
	{
	case 'A':
	{
		arr[11] = arr[12];
	}break;
	case 'C':
	{
		arr[11] = arr[12] + 64;
	}break;
	case 'G':
	{
		arr[11] = arr[12] + 128;
	}break;
	case 'T':
	{
		arr[11] = arr[12] + 192;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[10])
	{
	case 'A':
	{
		arr[10] = arr[11];
	}break;
	case 'C':
	{
		arr[10] = arr[11] + 256;
	}break;
	case 'G':
	{
		arr[10] = arr[11] + 512;
	}break;
	case 'T':
	{
		arr[10] = arr[11] + 768;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[9])
	{
	case 'A':
	{
		arr[9] = arr[10];
	}break;
	case 'C':
	{
		arr[9] = arr[10] + 1024;
	}break;
	case 'G':
	{
		arr[9] = arr[10] + 2048;
	}break;
	case 'T':
	{
		arr[9] = arr[10] + 3072;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[8])
	{
	case 'A':
	{
		arr[8] = arr[9];
	}break;
	case 'C':
	{
		arr[8] = arr[9] + 4096;
	}break;
	case 'G':
	{
		arr[8] = arr[9] + 8192;
	}break;
	case 'T':
	{
		arr[8] = arr[9] + 12288;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[7])
	{
	case 'A':
	{
		arr[7] = arr[8];
	}break;
	case 'C':
	{
		arr[7] = arr[8] + 16384;
	}break;
	case 'G':
	{
		arr[7] = arr[8] + 32768;
	}break;
	case 'T':
	{
		arr[7] = arr[8] + 49152;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[6])
	{
	case 'A':
	{
		arr[6] = arr[7];
	}break;
	case 'C':
	{
		arr[6] = arr[7] + 65536;
	}break;
	case 'G':
	{
		arr[6] = arr[7] + 131072;
	}break;
	case 'T':
	{
		arr[6] = arr[7] + 196608;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[5])
	{
	case 'A':
	{
		arr[5] = arr[6];
	}break;
	case 'C':
	{
		arr[5] = arr[6] + 262144;
	}break;
	case 'G':
	{
		arr[5] = arr[6] + 524288;
	}break;
	case 'T':
	{
		arr[5] = arr[6] + 786432;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[4])
	{
	case 'A':
	{
		arr[4] = arr[5];
	}break;
	case 'C':
	{
		arr[4] = arr[5] + 1048576;
	}break;
	case 'G':
	{
		arr[4] = arr[5] + 2097152;
	}break;
	case 'T':
	{
		arr[4] = arr[5] + 3145728;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[3])
	{
	case 'A':
	{
		arr[3] = arr[4];
	}break;
	case 'C':
	{
		arr[3] = arr[4] + 4194304;
	}break;
	case 'G':
	{
		arr[3] = arr[4] + 8388608;
	}break;
	case 'T':
	{
		arr[3] = arr[4] + 12582912;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[2])
	{
	case 'A':
	{
		arr[2] = arr[3];
	}break;
	case 'C':
	{
		arr[2] = arr[3] + 16777216;
	}break;
	case 'G':
	{
		arr[2] = arr[3] + 33554432;
	}break;
	case 'T':
	{
		arr[2] = arr[3] + 50331648;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[1])
	{
	case 'A':
	{
		arr[1] = arr[2];
	}break;
	case 'C':
	{
		arr[1] = arr[2] + 67108864;
	}break;
	case 'G':
	{
		arr[1] = arr[2] + 134217728;
	}break;
	case 'T':
	{
		arr[1] = arr[2] + 201326592;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	switch (s[0])
	{
	case 'A':
	{
		arr[0] = arr[1];
	}break;
	case 'C':
	{
		arr[0] = arr[1] + 268435456;
	}break;
	case 'G':
	{
		arr[0] = arr[1] + 536870912;
	}break;
	case 'T':
	{
		arr[0] = arr[1] + 805306368;
	}break;
	default:
	{
		return UINT_MAX;
	}
	}
	return arr[0];
}

uint32_t(*hashfunc)(const char*) = calc_hash_11;


/* Sorting and Manipulation Functions */
void merge(uint32_t* hashes, uint32_t* index, uint32_t* pos, int p, int m, int q) {
	int i, j, k;
	int n1 = m - p + 1;
	int n2 = q - m;

	// create temp arrays
	// unsigned int L[n1], R[n2];
	uint32_t* L = malloc(n1 * sizeof(uint32_t));
	MEMORY_FAILURE(L);
	uint32_t* R = malloc(n2 * sizeof(uint32_t));
	MEMORY_FAILURE(R);

	uint32_t* L_index = malloc(n1 * sizeof(uint32_t));
	MEMORY_FAILURE(L_index);
	uint32_t* R_index = malloc(n2 * sizeof(uint32_t));
	MEMORY_FAILURE(R_index);

	uint32_t* L_pos = malloc(n1 * sizeof(uint32_t));
	MEMORY_FAILURE(L_pos);
	uint32_t* R_pos = malloc(n2 * sizeof(uint32_t));
	MEMORY_FAILURE(R_pos);

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++) {
		L[i] = hashes[p + i];
		L_index[i] = index[p + i];
		L_pos[i] = pos[p + i];
	}
	for (j = 0; j < n2; j++) {
		R[j] = hashes[m + 1 + j];
		R_index[j] = index[m + 1 + j];
		R_pos[j] = pos[m + 1 + j];
	}

	// Merge the temp arrays back into hashes[p..q]
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = p; // Initial index of merged subarray
	while (i < n1 && j < n2) {
		if (L[i] <= R[j]) {
			hashes[k] = L[i];
			index[k] = L_index[i];
			pos[k] = L_pos[i];
			i++;
		}
		else {
			hashes[k] = R[j];
			index[k] = R_index[j];
			pos[k] = R_pos[j];
			j++;
		}
		k++;
	}
	while (i < n1) {
		hashes[k] = L[i];
		index[k] = L_index[i];
		pos[k] = L_pos[i];
		i++;
		k++;
	}
	while (j < n2) {
		hashes[k] = R[j];
		index[k] = R_index[j];
		pos[k] = R_pos[j];
		j++;
		k++;
	}
	free(L);
	free(R);
	free(L_index);
	free(R_index);
	free(L_pos);
	free(R_pos);
}

void mergeSort(uint32_t* hashes, uint32_t* index, uint32_t* pos, int p, int q) {
	// Merge sort takes three unsigned int arrays, possibly containing:
	// (1) pseudo-hashes;
	// (2) seqid indices;
	// (3) kmer start positions.
	if (p < q) {
		int m = p + (q - p) / 2;
		mergeSort(hashes, index, pos, p, m);
		mergeSort(hashes, index, pos, m + 1, q);
		merge(hashes, index, pos, p, m, q);
	}
}

int compare(const void* a, const void* b) {
	// Function for comparing values in structs
	HASH* hash_a = (HASH*)a;
	HASH* hash_b = (HASH*)b;
	if (hash_a->id < hash_b->id) return -1;
	if (hash_a->id > hash_b->id) return 1;
	return 0;
}

void revcmp(char* dest, char* src, size_t size) {
	// A simple reverse complement
	char* d = &dest[0];
	char* p = &src[size - 1];
	do {
		switch (*p)
		{
		case 'A':
		{
			*d = 'T';
		}break;
		case 'C':
		{
			*d = 'G';
		}break;
		case 'G':
		{
			*d = 'C';
		}break;
		case 'T':
		{
			*d = 'A';
		}break;
		default:
		{
			*d = *p;
		}
		}
		--p;
		++d;
		--size;
	} while (size);
	*d = 0x00;
}


/* Write Functions */
void WriteREAD(READ* read, FILE* fptr) {
	// Write a READ struct
	const char nln = '\n';
	fwrite(read->header, sizeof(char), strlen(read->header), fptr);
	fwrite(read->sequence, sizeof(char), strlen(read->sequence), fptr);
	fwrite(&nln, sizeof(char), 1, fptr);
	fwrite(read->placeholder, sizeof(char), strlen(read->placeholder), fptr);
	fwrite(read->quality, sizeof(char), strlen(read->quality), fptr);
}

void WriteBED(BED* bed, FILE* fptr) {
	// Write a BED struct
	fprintf(fptr, "%s\t%u\t%u\t%s\n",
		bed->chrom,
		bed->chromStart,
		bed->chromEnd,
		bed->name);
}

void WriteBEDPE(BEDPE* bedpe, FILE* fptr) {
	// Write a BEDPE struct
	const char nln = '\n';
	const char tab = '\t';

	fprintf(fptr, "%s\t%u\t%u\t%s\t%u\t%u\t%s\t%d\t%c\t%c\n",
		bedpe->chrom1,
		bedpe->start1,
		bedpe->end1,
		bedpe->chrom2,
		bedpe->start2,
		bedpe->end2,
		bedpe->name,
		bedpe->score,
		bedpe->strand1,
		bedpe->strand2);
}


/* Search Functions */
int UngappedSearch(const char* s1, const char* s2, const size_t s, unsigned short mismatches) {
	// Ungapped search
	mismatches = (mismatches * s) / 100;
	size_t score = 0;
	size_t i = s;
	do {
		--i;
        if(s1[i]=='N' && IGNORE_N) continue;
		if (s1[i] != s2[i]) ++score;
		if (score > mismatches) return 0;
	} while (i > 0);
	return 1;
}

int SeedSearch(HASH_TBL* hash_tbl, const char* read, size_t s, int* i) {
	uint32_t B;
	uint32_t index;
	const int lim = s - K + 1;  //read length - kmer length - 1
	do {
		if (*i > lim) return -1;
		B = hashfunc(&read[*i]);
		if (B == UINT_MAX) {  // could not create hash
			*i += STEP;
			continue;
		}
		uint32_t* p = bsearch(&B, hash_tbl->hash, hash_tbl->size, sizeof(*hash_tbl->hash), compare);
		if (p == NULL) { // hash not found
			*i += STEP;
			continue;
		}
		index = ((uintptr_t)p - (uintptr_t)hash_tbl->hash) / sizeof(*hash_tbl->hash);
		return (int)index;
	} while (1);
	return -1;
}

int UngappedSearch2(const char* s1, const char* s2, const size_t s, unsigned short mismatches) {
	// Ungapped search
	mismatches = (mismatches * s) / 100;
	size_t score = 0;
	size_t i = s;
	do {
		--i;
        if(s1[i]=='N' && IGNORE_N) continue;
		if (s1[i] != s2[i]) ++score;
		if (score > mismatches) return 0;
	} while (i > 0);
	return (1 + score);
}

int no_disjoin_inward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	int at1 = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1) return 0;
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		for (j = 0; j < R1_count; ++j) {
			int at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c) {
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;
					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					if (same_seq) {
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] <= R2_pos[i][1] + s2) {  // when it's NO_DISJOIN R1 must come before R2
									bedpe->chrom1 = index->fasta[i1].header;
									bedpe->start1 = R1_pos[j][1] - at1;
									bedpe->end1 = bedpe->start1 + s1;
									bedpe->chrom2 = index->fasta[i2].header;
									bedpe->start2 = R2_pos[i][1] - at2;
									bedpe->end2 = bedpe->start2 + s2;
									bedpe->score = 0;
									return 1;
								}
							}
						}
					}
					at2 += STEP;
				}
			}
		}
		at1 += STEP;
	}
	return 0;
}

int disjoin_inward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	int at1 = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1) return 0;
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		for (j = 0; j < R1_count; ++j) {
			int at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c) {
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;

					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					// give same seq priority
					int pos = i;
					if (same_seq) {
						// seq indices are sorted in the hash table!
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) break; // very important to break here
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] <= R2_pos[i][1] + s2) { // when it's NO_DISJOIN R1 must come before R2
									bedpe->chrom1 = index->fasta[i1].header;
									bedpe->start1 = R1_pos[j][1] - at1;
									bedpe->end1 = bedpe->start1 + s1;
									bedpe->chrom2 = index->fasta[i2].header;
									bedpe->start2 = R2_pos[i][1] - at2;
									bedpe->end2 = bedpe->start2 + s2;
									bedpe->score = 0;
									return 1;
								}
							}
						}
						// program moves here when a kmer was found on same seq, but seed extension did not produce a match
						// now check that SEQi+1 > SEQi, and not R1 R2 positions (since they're relative) 
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] < i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								bedpe->chrom1 = index->fasta[i1].header;
								bedpe->start1 = R1_pos[j][1] - at1;
								bedpe->end1 = bedpe->start1 + s1;
								bedpe->chrom2 = index->fasta[i2].header;
								bedpe->start2 = R2_pos[i][1] - at2;
								bedpe->end2 = bedpe->start2 + s2;
								if (i1 == i2) bedpe->score = 0;
								else bedpe->score = 1;
								return 1;
							}
						}
					}
					else {
						// not found on the same SEQ
						for (i = 0; i < R2_count; ++i) {
							i2 = R2_pos[i][0];
							if (i2 < i1) continue; // must be downstream!
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								bedpe->chrom1 = index->fasta[i1].header;
								bedpe->start1 = R1_pos[j][1] - at1;
								bedpe->end1 = bedpe->start1 + s1;
								bedpe->chrom2 = index->fasta[i2].header;
								bedpe->start2 = R2_pos[i][1] - at2;
								bedpe->end2 = bedpe->start2 + s2;
								if (i1 == i2) bedpe->score = 0;
								else bedpe->score = 1;
								return 1;
							}
						}
					}
					at2 += STEP;
				}
			}
		}
		at1 += STEP;
	}
	return 0;
}

int no_disjoin_outward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	int at1 = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1) return 0;
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		for (j = 0; j < R1_count; ++j) {
			int at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c) {
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;
					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					if (same_seq) {
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] >= R2_pos[i][1] + s2) {  // when it's NO_DISJOIN R1 must come after R2
									bedpe->chrom1 = index->fasta[i1].header;
									bedpe->start1 = R1_pos[j][1] - at1;
									bedpe->end1 = bedpe->start1 + s1;
									bedpe->chrom2 = index->fasta[i2].header;
									bedpe->start2 = R2_pos[i][1] - at2;
									bedpe->end2 = bedpe->start2 + s2;
									bedpe->score = 0;
									return 1;
								}
							}
						}
					}
					at2 += STEP;
				}
			}
		}
		at1 += STEP;
	}
	return 0;
}

int disjoin_outward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	int at1 = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1) return 0;
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		for (j = 0; j < R1_count; ++j) {
			int at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c) {
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;

					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					// give same seq priority
					int pos = i;
					if (same_seq) {
						// seq indices are sorted in the hash table!
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) break; // very important to break here
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] >= R2_pos[i][1] + s2) { // when it's NO_DISJOIN R1 must come after R2
									bedpe->chrom1 = index->fasta[i1].header;
									bedpe->start1 = R1_pos[j][1] - at1;
									bedpe->end1 = bedpe->start1 + s1;
									bedpe->chrom2 = index->fasta[i2].header;
									bedpe->start2 = R2_pos[i][1] - at2;
									bedpe->end2 = bedpe->start2 + s2;
									bedpe->score = 0;
									return 1;
								}
							}
						}
						// program moves here when a kmer was found on same seq, but seed extension did not produce a match
						// now check that SEQi+1 < SEQi, and not R1 R2 positions (since they're relative) 
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] > i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								bedpe->chrom1 = index->fasta[i1].header;
								bedpe->start1 = R1_pos[j][1] - at1;
								bedpe->end1 = bedpe->start1 + s1;
								bedpe->chrom2 = index->fasta[i2].header;
								bedpe->start2 = R2_pos[i][1] - at2;
								bedpe->end2 = bedpe->start2 + s2;
								if (i1 == i2) bedpe->score = 0;
								else bedpe->score = 1;
								return 1;
							}
						}
					}
					else {
						// not found on the same SEQ
						for (i = 0; i < R2_count; ++i) {
							i2 = R2_pos[i][0];
							if (i2 > i1) continue; // must be upstream!
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								bedpe->chrom1 = index->fasta[i1].header;
								bedpe->start1 = R1_pos[j][1] - at1;
								bedpe->end1 = bedpe->start1 + s1;
								bedpe->chrom2 = index->fasta[i2].header;
								bedpe->start2 = R2_pos[i][1] - at2;
								bedpe->end2 = bedpe->start2 + s2;
								if (i1 == i2) bedpe->score = 0;
								else bedpe->score = 1;
								return 1;
							}
						}
					}
					at2 += STEP;
				}
			}
		}
		at1 += STEP;
	}
	return 0;
}

int eval_inward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	BEDPE* beds = malloc(EVALS * sizeof * beds);
	MEMORY_FAILURE(beds);
	int* scores = malloc(EVALS * sizeof(int));
	MEMORY_FAILURE(scores);
	int e = 0;
	for(e=0;e<EVALS;++e)
	{
		scores[e] = BUF_SIZE;
	}
	e = 0;
	int at1 = 0;
	int exit = 0;
	int ever_found = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1){
			if(ever_found)
				continue;
		}
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		int at2 = 0;
		int found = 0;
		for (j = 0; j < R1_count; ++j) {
			if(found){
				++e;
				found = 0;
			}
			if(e==EVALS) {
				exit = 1; break;
			}
			at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch2(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c){
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						if(ever_found)
							continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;

					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					if (same_seq) {
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch2(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] <= R2_pos[i][1] + s2) {  // when it's NO_DISJOIN R1 must come before R2
									beds[e].chrom1 = index->fasta[i1].header;
									beds[e].start1 = R1_pos[j][1] - at1;
									beds[e].end1 = beds[e].start1 + s1;
									beds[e].chrom2 = index->fasta[i2].header;
									beds[e].start2 = R2_pos[i][1] - at2;
									beds[e].end2 = beds[e].start2 + s2;
									beds[e].score = 0;
									found = 1;
									ever_found = 1;
									scores[e] = c+d;
									break;
								}
							}
						}
					}
					at2 += STEP;
					if(found) break;
				}
			}
		}
		at1 += STEP;
		if(exit) break;
	}
	if(ever_found){
		int val = scores[0];
		int win = 0;
		for(e=0;e<EVALS;++e){
			if(scores[e]<val){
				val = scores[e];
				win = e;
			}
		}
		bedpe->chrom1 = beds[win].chrom1;
		bedpe->start1 = beds[win].start1;
		bedpe->end1 = beds[win].end1;
		bedpe->chrom2 = beds[win].chrom2;
		bedpe->start2 = beds[win].start2;
		bedpe->end2 = beds[win].end2;
		bedpe->score = beds[win].score;
		free(scores);
		free(beds);
		return 1;	
	}
	free(scores);
	free(beds);
	return 0;
}

int eval_outward_s(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) {
	BEDPE* beds = malloc(EVALS * sizeof * beds);
	MEMORY_FAILURE(beds);
	int* scores = malloc(EVALS * sizeof(int));
	MEMORY_FAILURE(scores);
	int e = 0;
	for(e=0;e<EVALS;++e)
	{
		scores[e] = BUF_SIZE;
	}
	e = 0;
	int at1 = 0;
	int exit = 0;
	int ever_found = 0;
	while (at1 < s1 - K + 1) {
		// bsearch R1. If found bsearch R2. If found enter while loop
		int a, b;
		a = SeedSearch(hash_tbl, read_R1, s1, &at1);
		if (a == -1){
			if(ever_found)
				continue;
		}
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;
		R1_pos = hash_tbl->hash[a].pos;
		int c, d;
		int j, i, i1, i2;
		int at2 = 0;
		int found = 0;
		for (j = 0; j < R1_count; ++j) {
			if(found){
				++e;
				found = 0;
			}
			if(e==EVALS) {
				exit = 1; break;
			}
			at2 = 0;
			i1 = R1_pos[j][0];
			if ((R1_pos[j][1] - at1) < 0) continue;
			if ((R1_pos[j][1] - at1) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch2(&index->fasta[i1].sequence[R1_pos[j][1] - at1], read_R1, s1, MISMATCH);
			if (c){
				// First in pair was found
				while (at2 < s2 - K + 1) {
					// find occurrence of kmer in sequencing read
					b = SeedSearch(hash_tbl, read_R2, s2, &at2);
					if (b == -1){
						if(ever_found)
							continue;
					}
					uint32_t R2_count;
					uint32_t** R2_pos;
					R2_count = hash_tbl->hash[b].size;  // collisions
					R2_pos = hash_tbl->hash[b].pos;

					// check if k on R2 belongs to the same sequence as R1
					int same_seq = 0;
					for (i = 0; i < R2_count; ++i) {
						if (R2_pos[i][0] == i1) {
							i2 = R2_pos[i][0];
							same_seq = 1;
							break;
						}
					}
					if (same_seq) {
						for (i=i; i < R2_count; ++i) {
							if (R2_pos[i][0] != i1) continue;
							i2 = R2_pos[i][0];
							if ((R2_pos[i][1] - at2) < 0) { continue; }
							if ((R2_pos[i][1] - at2) + s2 > index->fasta[i2].size) { continue; }
							if (index->fasta[i2].size < s2) { continue; }
							d = UngappedSearch2(&index->fasta[i2].sequence[R2_pos[i][1] - at2], read_R2, s2, MISMATCH);
							if (d) {
								// mate matched
								if (R1_pos[j][1] >= R2_pos[i][1] + s2) {  // when it's NO_DISJOIN R2 must come before R1
									beds[e].chrom1 = index->fasta[i1].header;
									beds[e].start1 = R1_pos[j][1] - at1;
									beds[e].end1 = beds[e].start1 + s1;
									beds[e].chrom2 = index->fasta[i2].header;
									beds[e].start2 = R2_pos[i][1] - at2;
									beds[e].end2 = beds[e].start2 + s2;
									beds[e].score = 0;
									found = 1;
									ever_found = 1;
									scores[e] = c+d;
									break;
								}
							}
						}
					}
					at2 += STEP;
					if(found) break;
				}
			}
		}
		at1 += STEP;
		if(exit) break;
	}
	if(ever_found){
		int val = scores[0];
		int win = 0;
		for(e=0;e<EVALS;++e){
			if(scores[e]<val){
				val = scores[e];
				win = e;
			}
		}
		bedpe->chrom1 = beds[win].chrom1;
		bedpe->start1 = beds[win].start1;
		bedpe->end1 = beds[win].end1;
		bedpe->chrom2 = beds[win].chrom2;
		bedpe->start2 = beds[win].start2;
		bedpe->end2 = beds[win].end2;
		bedpe->score = beds[win].score;
		free(scores);
		free(beds);
		return 1;	
	}
	free(scores);
	free(beds);
	return 0;
}

int singlet(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, size_t s1, BED* bed) {
	int at = 0;
	while (at < s1 - K + 1) {
		int a, c;
		a = SeedSearch(hash_tbl, read_R1, s1, &at); // a is the index of the hash (kmer) with a perfect match
		if (a == -1) return 0;
		uint32_t R1_count;
		uint32_t** R1_pos;
		R1_count = hash_tbl->hash[a].size;  //hash_tbl->hash[a].size is the number of collisions
		R1_pos = hash_tbl->hash[a].pos; //hash_tbl->hash[a].pos is the array with all positions with that kmer in db
		int i, i1;
		for (i = 0; i < R1_count; ++i) {
			i1 = R1_pos[i][0];  // i1 is the replicon index
			if ((R1_pos[i][1] - at) < 0) continue; // R1_pos[i][1] is the position 
			if ((R1_pos[i][1] - at) + s1 > index->fasta[i1].size) continue;
			if (index->fasta[i1].size < s1) continue;
			c = UngappedSearch(&index->fasta[i1].sequence[R1_pos[i][1] - at], read_R1, s1, MISMATCH);
			if (c) {
				bed->chrom = index->fasta[i1].header;
				bed->chromStart = R1_pos[i][1] - at;
				bed->chromEnd = R1_pos[i][1] + at + s1;
				return 1;
			}
		}
		at += STEP;
	}
	return 0;
}

int(*inward_searchfunc)(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) = no_disjoin_inward_s;
int(*outward_searchfunc)(INDEX* index, HASH_TBL* hash_tbl, char* read_R1, char* read_R2, size_t s1, size_t s2, BEDPE* bedpe) = no_disjoin_outward_s;

/* Unstranded Libraries */
void* LibIUSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int isf, isr;
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R2.sequence, s2);
			BEDPE bedpe;

			isf = inward_searchfunc(thread->index, thread->hash_tbl, read_R1.sequence, reverse_str, s1, s2, &bedpe);
			if (isf) {
				bedpe.strand1 = '+';
				bedpe.strand2 = '-';
				matched = 1;
			}
			else {
				revcmp(reverse_str, read_R1.sequence, s1);
				isr = inward_searchfunc(thread->index, thread->hash_tbl, read_R2.sequence, reverse_str, s2, s1, &bedpe);
				if (isr) {
					bedpe.strand1 = '-';
					bedpe.strand2 = '+';
					matched = 1;
				}
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibOUSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int osf, osr;
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R2.sequence, s2);
			BEDPE bedpe;
			osf = outward_searchfunc(thread->index, thread->hash_tbl, read_R1.sequence, reverse_str, s1, s2, &bedpe);
			if (osf) {
				bedpe.strand1 = '-';
				bedpe.strand2 = '+';
				matched = 1;
			}
			else {
				revcmp(reverse_str, read_R1.sequence, s1);
				osr = outward_searchfunc(thread->index, thread->hash_tbl, read_R2.sequence, reverse_str, s2, s1, &bedpe);
				{
					if (osr)
					{
						bedpe.strand1 = '+';
						bedpe.strand2 = '-';
						matched = 1;
					}
				}
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibUSearchR1(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	FILE* fout_R1;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0)
		{
			strcpy(read_R1.header, buf1);
		}
		if (line == 1)
		{
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
		}
		if (line == 2)
		{
			strcpy(read_R1.placeholder, buf1);
		}
		if (line == 3)
		{
			strcpy(read_R1.quality, buf1);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s1 = strlen(read_R1.sequence);
			sgl = singlet(thread->index, thread->hash_tbl, read_R1.sequence, s1, &bed);
			if (sgl) {
				matched = 1;
			}
			else {
				revcmp(reverse_str, read_R1.sequence, s1);
				sgl = singlet(thread->index, thread->hash_tbl, reverse_str, s1, &bed);
				if (sgl) matched = 1;

			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bed.name = read_R1.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	fclose(fout_R1);
	fclose(BEDptr);
	return 0;
}

void* LibUSearchR2(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R2;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R2;
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s2 = strlen(read_R2.sequence);
			sgl = singlet(thread->index, thread->hash_tbl, read_R2.sequence, s2, &bed);
			if (sgl) {
				matched = 1;
			}
			else {
				revcmp(reverse_str, read_R2.sequence, s2);
				sgl = singlet(thread->index, thread->hash_tbl, reverse_str, s2, &bed);
				if (sgl) matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R2.header[strcspn(read_R2.header, "\n")] = 0;
					bed.name = read_R2.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R2);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}


/* Stranded Libraries */
void* LibISFSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int isf;
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R2.sequence, s2);
			BEDPE bedpe;

			isf = inward_searchfunc(thread->index, thread->hash_tbl, read_R1.sequence, reverse_str, s1, s2, &bedpe);
			if (isf) {
				bedpe.strand1 = '+';
				bedpe.strand2 = '-';
				matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibISRSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int isr;
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R1.sequence, s1);
			BEDPE bedpe;
			isr = inward_searchfunc(thread->index, thread->hash_tbl, read_R2.sequence, reverse_str, s2, s1, &bedpe); {
				if (isr) {
					bedpe.strand1 = '-';
					bedpe.strand2 = '+';
					matched = 1;
				}
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibOSFSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int osf;   // lib IU tries both
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R2.sequence, s2);
			BEDPE bedpe;
			osf = outward_searchfunc(thread->index, thread->hash_tbl, read_R1.sequence, reverse_str, s1, s2, &bedpe);
			if (osf) {
				bedpe.strand1 = '-';
				bedpe.strand2 = '+';
				matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibOSRSearch(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	gzFile fastq_R2;
	FILE* fout_R1;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	READ read_R2;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL && gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			int osr;
			size_t s1 = strlen(read_R1.sequence);
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R1.sequence, s1);
			BEDPE bedpe;
			osr = outward_searchfunc(thread->index, thread->hash_tbl, read_R2.sequence, reverse_str, s1, s2, &bedpe);
			if (osr) {
				bedpe.strand1 = '+';
				bedpe.strand2 = '-';
				matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bedpe.name = read_R1.header;
					WriteBEDPE(&bedpe, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	gzclose(fastq_R2);
	fclose(fout_R1);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}


/* Single Ends Libraries */
void* LibSFSearchR1(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	FILE* fout_R1;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	READ read_R1;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s1 = strlen(read_R1.sequence);

			sgl = singlet(thread->index, thread->hash_tbl, read_R1.sequence, s1, &bed);
			if (sgl) {
				matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bed.name = read_R1.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	fclose(fout_R1);
	fclose(BEDptr);
	return 0;
}

void* LibSFSearchR2(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R2;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf2[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R2;
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s2 = strlen(read_R2.sequence);
			revcmp(reverse_str, read_R2.sequence, s2);
			sgl = singlet(thread->index, thread->hash_tbl, reverse_str, s2, &bed); {
				if (sgl) matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF)
					{
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R2.header[strcspn(read_R2.header, "\n")] = 0;
					bed.name = read_R2.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R2);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

void* LibSRSearchR1(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R1;
	FILE* fout_R1;
	FILE* BEDptr;
	fastq_R1 = gzopen(thread->params->input_R1, "rb");
	MEMORY_FAILURE(fastq_R1);
	fout_R1 = fopen(thread->temp_out_R1, "w");
	MEMORY_FAILURE(fout_R1);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf1[READ_SIZE];
	char reverse_str[READ_SIZE];
	READ read_R1;
	gzseek(fastq_R1, thread->bytes_R1_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R1, buf1, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R1.header, buf1);
		}
		if (line == 1) {
			buf1[strcspn(buf1, "\r\n")] = 0; strcpy(read_R1.sequence, buf1);
		}
		if (line == 2) {
			strcpy(read_R1.placeholder, buf1);
		}
		if (line == 3) {
			strcpy(read_R1.quality, buf1);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s1 = strlen(read_R1.sequence);
			revcmp(reverse_str, read_R1.sequence, s1);
			sgl = singlet(thread->index, thread->hash_tbl, reverse_str, s1, &bed); {
				if (sgl) matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R1, fout_R1);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R1.header[strcspn(read_R1.header, " \n")] = 0;
					bed.name = read_R1.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R1);
	fclose(fout_R1);
	fclose(BEDptr);
	return 0;
}

void* LibSRSearchR2(void* arg) {
	THREAD* thread = arg;
	gzFile fastq_R2;
	FILE* fout_R2;
	FILE* BEDptr;
	fastq_R2 = gzopen(thread->params->input_R2, "rb");
	MEMORY_FAILURE(fastq_R2);
	fout_R2 = fopen(thread->temp_out_R2, "w");
	MEMORY_FAILURE(fout_R2);
	BEDptr = fopen(thread->temp_out_BED, "w");
	MEMORY_FAILURE(BEDptr);
	char buf2[READ_SIZE];
	READ read_R2;
	gzseek(fastq_R2, thread->bytes_R2_start, SEEK_SET);
	int total = 0;
	int line = 0;
	while (gzgets(fastq_R2, buf2, READ_SIZE) != NULL) {
		int matched = 0;
		if (line == 0) {
			strcpy(read_R2.header, buf2);
		}
		if (line == 1) {
			buf2[strcspn(buf2, "\r\n")] = 0; strcpy(read_R2.sequence, buf2);
		}
		if (line == 2) {
			strcpy(read_R2.placeholder, buf2);
		}
		if (line == 3) {
			strcpy(read_R2.quality, buf2);
		}
		++line;
		if (line == 4) {
			BED bed;
			int sgl;
			size_t s2 = strlen(read_R2.sequence);
			sgl = singlet(thread->index, thread->hash_tbl, read_R2.sequence, s2, &bed);
			if (sgl) {
				matched = 1;
			}
			if (!matched) {
				if (FASTQ_OUT) {
					if (DIFF) {
						WriteREAD(&read_R2, fout_R2);
					}
				}
			}
			else {
				if (FASTQ_OUT) {
					if (!DIFF) {
						WriteREAD(&read_R2, fout_R2);
					}
				}
				if (BED_OUT && !DIFF) {
					read_R2.header[strcspn(read_R2.header, "\n")] = 0;
					bed.name = read_R2.header;
					WriteBED(&bed, BEDptr);
				}
			}
			matched = 0;
			total += line;
			if (total >= thread->lines) break;
			line = 0;
		}
	}
	gzclose(fastq_R2);
	fclose(fout_R2);
	fclose(BEDptr);
	return 0;
}

/* Multithreading Functions */
int MultiThreadManager(PARAMS* params, INDEX* index, HASH_TBL* hash_tbl, void(*func)) {
	VERBOSE_MSG('L', "Buffering data, please wait ");
	long* bytes_R1_end = calloc(AVAIL_THREADS, sizeof(long));
	long* bytes_R2_end = calloc(AVAIL_THREADS, sizeof(long));
	long* bytes_R1_start = calloc(AVAIL_THREADS, sizeof(long));
	long* bytes_R2_start = calloc(AVAIL_THREADS, sizeof(long));
	size_t* GZTELL_LINES = calloc(AVAIL_THREADS + 1, sizeof(size_t));
	unsigned int i = 0;
	long filesize = 0;
	char buf[BUF_SIZE];
	if (params->input_R1 != NULL && params->input_R2 != NULL) {
		gzFile f;
		f = gzopen(params->input_R1, "rb");
		while (gzgets(f, buf, 8192) != NULL) {}
		filesize = gztell(f);
		gzrewind(f);
		long approx = filesize / AVAIL_THREADS;
		long tot = approx;
		size_t this_line = 0;
		long current = 0;
		while (gzgets(f, buf, 8192) != NULL) {
			++this_line;
			current = gztell(f);
			if (current >= tot) {
				if (this_line % 4 == 0) {
					tot += approx;
					bytes_R1_end[i] = current;
					GZTELL_LINES[i] = this_line;
					this_line = 0;
					++i;
				}
			}
		}
		bytes_R1_start[0] = 0;
		for (i = 0; i < AVAIL_THREADS - 1; ++i) {
			bytes_R1_start[i + 1] = bytes_R1_end[i];
		}
		gzclose(f);
		this_line = 0;
		i = 0;
		f = gzopen(params->input_R2, "rb");
		while (gzgets(f, buf, 8192) != NULL) {
			++this_line;
			if (this_line == GZTELL_LINES[i])
			{
				bytes_R2_end[i] = gztell(f);
				this_line = 0;
				++i;
			}
		}
		bytes_R2_start[0] = 0;
		for (i = 0; i < AVAIL_THREADS - 1; ++i) {
			bytes_R2_start[i + 1] = bytes_R2_end[i];
		}
		gzclose(f);
	}
	else {
		char* file = params->input_R2;
		if (params->input_R2 == NULL) file = params->input_R1;
		gzFile f;
		f = gzopen(file, "rb");
		while (gzgets(f, buf, 8192) != NULL) {}
		filesize = gztell(f);
		gzrewind(f);
		long approx = filesize / AVAIL_THREADS;
		long tot = approx;
		size_t this_line = 0;
		long current = 0;
		while (gzgets(f, buf, 8192) != NULL) {
			++this_line;
			current = gztell(f);
			if (current >= tot) {
				if (this_line % 4 == 0) {
					tot += approx;
					bytes_R1_end[i] = current;
					bytes_R2_end[i] = current;
					GZTELL_LINES[i] = this_line;
					this_line = 0;
					++i;
				}
			}
		}
		for (i = 0; i < AVAIL_THREADS - 1; ++i) {
			bytes_R1_start[i + 1] = bytes_R1_end[i];
			bytes_R2_start[i + 1] = bytes_R2_end[i];
		}
		gzclose(f);
	}
	VERBOSE_MSG('L', "[OK]\n");
	VERBOSE_MSG('L', "Starting Execution\n");
	THREAD* threads;
	threads = malloc(AVAIL_THREADS * sizeof * threads);
	MEMORY_FAILURE(threads);
	pthread_t* tids = malloc(AVAIL_THREADS * sizeof * tids);
	MEMORY_FAILURE(tids);
	int t;
	for (t = 0; t < AVAIL_THREADS; ++t) {
		// get thread id to string
		int length = 11;
		char* thr = malloc(length);
		MEMORY_FAILURE(thr);
		snprintf(thr, length, "%d", t);

		threads[t].params = params;
		threads[t].index = index;
		threads[t].hash_tbl = hash_tbl;

		strcpy(threads[t].temp_out_R1, params->basename);
		strcat(threads[t].temp_out_R1, "_R1.thread");
		strcat(threads[t].temp_out_R1, thr);

		strcpy(threads[t].temp_out_R2, params->basename);
		strcat(threads[t].temp_out_R2, "_R2.thread");
		strcat(threads[t].temp_out_R2, thr);

		strcpy(threads[t].temp_out_BED, params->basename);
		strcat(threads[t].temp_out_BED, "_BED.thread");
		strcat(threads[t].temp_out_BED, thr);

		threads[t].bytes_R1_start = bytes_R1_start[t];
		threads[t].bytes_R2_start = bytes_R2_start[t];
		threads[t].lines = GZTELL_LINES[t];

		pthread_create(&tids[t], NULL, func, &threads[t]);
		free(thr);
	}
	for (t = 0; t < AVAIL_THREADS; ++t) {
		if (VERBOSE) fprintf(stderr, "Joining thread %d\n", t);
		pthread_join(tids[t], NULL);
	}
	VERBOSE_MSG('L', "Finalizing results ");
	//Merge files
	FILE* R1;
	R1 = fopen(params->output_R1, "w");
	MEMORY_FAILURE(R1);

	FILE* R2;
	R2 = fopen(params->output_R2, "w");
	MEMORY_FAILURE(R2);

	FILE* BEDptr;
	BEDptr = fopen(params->output_BED, "w");
	MEMORY_FAILURE(BEDptr);

	for (t = 0; t < AVAIL_THREADS; ++t) {
		FILE* fp1 = fopen(threads[t].temp_out_R1, "r");
		if (fp1 == NULL) continue;
		while (fgets(buf, BUF_SIZE, fp1) != NULL) fwrite(buf, sizeof(char), strlen(buf), R1);
		fclose(fp1);
		remove(threads[t].temp_out_R1);
	}

	for (t = 0; t < AVAIL_THREADS; ++t) {
		FILE* fp2 = fopen(threads[t].temp_out_R2, "r");
		if (fp2 == NULL) continue;
		while (fgets(buf, BUF_SIZE, fp2) != NULL) fwrite(buf, sizeof(char), strlen(buf), R2);
		fclose(fp2);
		remove(threads[t].temp_out_R2);
	}

	for (t = 0; t < AVAIL_THREADS; ++t) {
		FILE* fp3 = fopen(threads[t].temp_out_BED, "r");
		if (fp3 == NULL) continue;
		while (fgets(buf, BUF_SIZE, fp3) != NULL) fwrite(buf, sizeof(char), strlen(buf), BEDptr);
		fclose(fp3);
		remove(threads[t].temp_out_BED);
	}

	fclose(R1);
	fclose(R2);
	fclose(BEDptr);

	R1 = fopen(params->output_R1, "r");
	MEMORY_FAILURE(R1);
	fseek(R1, 0L, SEEK_END);
	if (ftell(R1) == 0) {
		fclose(R1);
		int r = remove(params->output_R1);
	}
	else fclose(R1);

	R2 = fopen(params->output_R2, "r");
	MEMORY_FAILURE(R2);
	fseek(R2, 0L, SEEK_END);
	if (ftell(R2) == 0) {
		fclose(R2);
		int r = remove(params->output_R2);
	}
	else fclose(R2);

	BEDptr = fopen(params->output_BED, "r");
	MEMORY_FAILURE(BEDptr);
	fseek(BEDptr, 0L, SEEK_END);
	if (ftell(BEDptr) == 0) {
		fclose(BEDptr);
		int r = remove(params->output_BED);
	}
	else fclose(BEDptr);

	free(threads);
	free(tids);

	free(bytes_R1_end);
	free(bytes_R2_end);
	free(bytes_R1_start);
	free(bytes_R2_start);
	free(GZTELL_LINES);
	VERBOSE_MSG('L', "[OK]\n");
	return 0;
}


/* Other */
void PrintSizes() {
	if (VERBOSE) {
		fprintf(stderr, "Your System Uses:\n");
		fprintf(stderr, "  %zu bytes for type uint_fast64_t\n", sizeof(uint_fast64_t));
		fprintf(stderr, "  %zu bytes for type size_t\n", sizeof(size_t));
		fprintf(stderr, "  %zu bytes for type uint32_t\n", sizeof(uint32_t));
		fprintf(stderr, "  %zu bytes for type int\n", sizeof(int));
		fprintf(stderr, "  %zu bytes for type char\n\n", sizeof(char));
	}
}

void ReadDatabase(PARAMS* params, INDEX* index, HASH_TBL* hash_tbl) {
	gzFile fptr;
	char buffer[BUF_SIZE];
	unsigned int headers_count = 0;

	fptr = gzopen(params->db, "rb");
	MEMORY_FAILURE(fptr);

	// find out how much memory needs to be allocated
	while (gzgets(fptr, buffer, BUF_SIZE) != 0) {
		if (buffer[0] == '>') {
			++headers_count;
		}
	}
	gzclose(fptr);
	if (!headers_count) {
		fprintf(stderr, "[ERROR] Could not parse FASTA database\n");
		exit(1);
	}

	// allocate memory for each FASTA database (array of INDEX structs)
	index->fasta = malloc(headers_count * sizeof * index->fasta);
	MEMORY_FAILURE(index->fasta);

	// fill database struct
	index->size = headers_count;

	// declare counters
	size_t i = 0;
	size_t j = 0;
	size_t n = 0;
	int begin_of_file = 1;
	size_t nucl_count = 0;
	double GC_count = 0.0;

	// reopen file handle and copy information
	fptr = gzopen(params->db, "rb");
	MEMORY_FAILURE(fptr);

	while (gzgets(fptr, buffer, BUF_SIZE) != 0) {
		if (buffer[0] == '>') {
			//strcpy(index->fasta[i].header, strtok(strtok(strtok(buffer, ">"), " "), "\n"));
			buffer[strcspn(buffer, " \r\n")] = 0;
			strcpy(index->fasta[i].header, &buffer[1]);
			if (!begin_of_file) {
				index->fasta[i - 1].size = nucl_count;
				index->fasta[i - 1].GC_content = GC_count;
				nucl_count = 0;
				GC_count = 0;
			}
			++i;
			begin_of_file = 0;
		}
		else {
			buffer[strcspn(buffer, "\r\n")] = 0;
			nucl_count += strlen(buffer);
			for (n = 0; n < strlen(buffer); ++n) {
				//++nucl_positions[i - 1];
				if (buffer[n] == 'C') {
					++GC_count;
				}
				if (buffer[n] == 'G') {
					++GC_count;
				}
			}
		}
	}
	gzclose(fptr);
	// store last out of loop
	index->fasta[i - 1].size = nucl_count;
	index->fasta[i - 1].GC_content = GC_count;
	nucl_count = 0;
	GC_count = 0;
	// make sure headers are unique
	for (j = 0; j < index->size; ++j) {
		for (i = 0; i < index->size; ++i) {
			if (j != i && strcmp(index->fasta[j].header, index->fasta[i].header) == 0) {
				fprintf(stderr, "[ERROR] database not contain identical headers\n");
				exit(1);
			}
		}
	}

	// allocate memory for storing whole database
	for (i = 0; i < index->size; ++i) {
		if (index->fasta[i].size < K) {
			fprintf(stderr, "[ERROR] Sequence shorter than kmer in database (%s)\n", index->fasta[i].header);
			exit(1);
		}
		index->fasta[i].sequence = malloc((1 + index->fasta[i].size) * sizeof(char));
		MEMORY_FAILURE(index->fasta[i].sequence);
		memset(index->fasta[i].sequence, 0x00, (1 + index->fasta[i].size) * sizeof(char));
	}

	//reopen file and copy sequence
	fptr = gzopen(params->db, "rb");
	MEMORY_FAILURE(fptr);

	i = 0; j = 0;
	while (gzgets(fptr, buffer, BUF_SIZE) != 0) {
		if (buffer[0] == '>') {
			++i;
			j = 0;
		}
		else {
			buffer[strcspn(buffer, "\r\n")] = 0;
			if (MASK_LOWER) {
				for (n = 0; n < strlen(buffer); ++n) {
					index->fasta[i - 1].sequence[j] = buffer[n];  // store sequence
					++j;
				}
			}
			else {
				for (n = 0; n < strlen(buffer); ++n) {
					index->fasta[i - 1].sequence[j] = toupper(buffer[n]);  // store sequence
					++j;
				}
			}
		}
	}

	// log information
	if (VERBOSE) {
		fprintf(stderr, "Database composition:\n");
		for (i = 0; i < index->size; ++i) {
			fprintf(stderr, "  [*] %s;bp:%ld;GC:%3.2f\n",
				index->fasta[i].header, index->fasta[i].size,
				index->fasta[i].GC_content / index->fasta[i].size * 100);
		}
	}
	gzclose(fptr);

	VERBOSE_MSG('L', "Hashing database sequence ");

	uint32_t* hashes = NULL;
	uint32_t* indices = NULL;
	uint32_t* pos = NULL;

	// pre-calculate needed memory
	unsigned int arr_size = 0;
	for (i = 0; i < index->size; ++i) {
		switch (K)
		{
		case 9:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 9; ++n) {
				++arr_size;
			}
		} break;
		case 11:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 11; ++n) {
				++arr_size;
			}
		} break;
		case 13:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 13; ++n) {
				++arr_size;
			}
		} break;
		case 15:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 15; ++n) {
				++arr_size;
			}
		} break;
		}
	}

	hashes = calloc(arr_size, sizeof(uint32_t));
	MEMORY_FAILURE(hashes);

	indices = calloc(arr_size, sizeof(uint32_t));
	MEMORY_FAILURE(indices);

	pos = calloc(arr_size, sizeof(uint32_t));
	MEMORY_FAILURE(pos);

	size_t c = 0; //global counter
	for (i = 0; i < index->size; ++i) {
		size_t n = 0;
		uint32_t h;
		switch (K)
		{
		case 9:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 9; ++n) {
				h = calc_hash_9(&(index->fasta[i].sequence[n]));
				if (UINT_MAX - h == 0) continue;

				hashes[c] = h;
				indices[c] = i;
				pos[c] = n;
				++c;
			}
		} break;
		case 11:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 11; ++n) {
				h = calc_hash_11(&(index->fasta[i].sequence[n]));
				if (UINT_MAX - h == 0) continue;
				hashes[c] = h;
				indices[c] = i;
				pos[c] = n;
				++c;
			}
		} break;
		case 13:
		{
			for (n = 0; n < 1 + index->fasta[i].size - 13; ++n) {
				h = calc_hash_13(&(index->fasta[i].sequence[n]));
				if (UINT_MAX - h == 0) continue;
				hashes[c] = h;
				indices[c] = i;
				pos[c] = n;
				++c;
			}
		} break;
		case 15:
		{

			for (n = 0; n < 1 + index->fasta[i].size - 15; ++n) {
				h = calc_hash_15(&(index->fasta[i].sequence[n]));
				if (UINT_MAX - h == 0) continue;
				hashes[c] = h;
				indices[c] = i;
				pos[c] = n;
				++c;
			}

		} break;
		}

	}

	mergeSort(hashes, indices, pos, 0, arr_size - 1); //fundamental for binary search
	//for(i=0;i<arr_size;++i){printf("%u\n", hashes[i]);} // for visualizing pseudo hashes
	size_t lines = 0; // for shrunk copy
	unsigned int last = UINT_MAX;
	c = 0;
	for (c = 0; c < arr_size; ++c) {
		if (hashes[c] != last) {
			last = hashes[c];
			++lines;
		}
	}
	VERBOSE_MSG('L', "[OK]\n");
	if (VERBOSE) fprintf(stderr, "%zu pseudo-hashes were produced\n", lines);
	if (VERBOSE) fprintf(stderr, "%d sequences were loaded in memory\n", headers_count);
	last = UINT_MAX;

	// now allocate memory for hash_tbl->hash*

	hash_tbl->hash = malloc(lines * sizeof * hash_tbl->hash);
	MEMORY_FAILURE(hash_tbl->hash);

	// iterate again to get sizes for mallocs
	uint32_t* mem_hash;
	mem_hash = malloc(lines * sizeof(uint32_t));
	MEMORY_FAILURE(mem_hash);
	size_t l = lines;
	while (l) {
		--l;
		mem_hash[l] = 1;  //there is at least one element
	}
	/*
	Visual Representation of what is happening here
	hashes stores sorted pseudo-hashes (unsigned integers)

	[0]
	[0]
	[0]
	[3]
	[4]
	[5]
	[5]
	[.]

	Collisions are grouped:

	hash  coll
	[0]   [3] (seqid and kmer start position are also stored)
	[3]   [1] ( . . . )
	[4]   [1]
	[5]   [2]
	[.]

	*/
	j = 0;
	int start = 1;
	last = hashes[0];
	mem_hash[j] = 0;
	for (c = 0; c < arr_size; ++c) {
		if (hashes[c] == last) mem_hash[j]++;
		else {
			hash_tbl->hash[j].id = 0;
			hash_tbl->hash[j].size = mem_hash[j]; //initialize memory
			hash_tbl->hash[j].pos = malloc(mem_hash[j] * sizeof(uint32_t*));
			MEMORY_FAILURE(hash_tbl->hash[j].pos);
			uint32_t i = 0;
			for (i = 0; i < mem_hash[j]; ++i) {
				hash_tbl->hash[j].pos[i] = malloc(2 * sizeof(uint32_t));
				MEMORY_FAILURE(hash_tbl->hash[j].pos[i]);
				hash_tbl->hash[j].pos[i][0] = 0;
				hash_tbl->hash[j].pos[i][1] = 0;
			}
			++j;
		}
		last = hashes[c];
	}
	if (hashes[c - 2] == hashes[c - 1]) mem_hash[j]++;
	hash_tbl->hash[j].id = hashes[c - 1];
	hash_tbl->hash[j].size = mem_hash[j];
	hash_tbl->hash[j].pos = malloc(mem_hash[j] * sizeof(*hash_tbl->hash[j].pos));
	MEMORY_FAILURE(hash_tbl->hash[j].pos);
	i = 0;
	for (i = 0; i < mem_hash[j]; ++i) {
		hash_tbl->hash[j].pos[i] = malloc(2 * sizeof(uint32_t));
		MEMORY_FAILURE(hash_tbl->hash[j].pos[i]);
		hash_tbl->hash[j].pos[i][0] = 0;
		hash_tbl->hash[j].pos[i][1] = 0;
	}

	last = UINT_MAX;
	j = 0;
	i = 0;
	start = 1;
	for (c = 0; c < arr_size; ++c) {
		if (hashes[c] != last) {
			if (!start) ++j;
			hash_tbl->hash[j].id = hashes[c];
			last = hashes[c];
			i = 0;
			hash_tbl->hash[j].pos[i][0] = (uint32_t)indices[c];
			hash_tbl->hash[j].pos[i][1] = (uint32_t)pos[c];
			start = 0;
		}
		else {
			++i;
			//hash_tbl->hash[j].id = hashes[c];
			hash_tbl->hash[j].pos[i][0] = (uint32_t)indices[c];
			hash_tbl->hash[j].pos[i][1] = (uint32_t)pos[c];
		}
	}

	hash_tbl->size = lines;

	// if you want to take a look(p) at the pseudohashes. Sequences indices are sorted
	// for(c=0;c<hash_tbl->size;++c){
	//   for (i=0;i<hash_tbl->hash[c].size;++i){
	//    printf("%d %d %d\n", 
	//   hash_tbl->hash[c].id, 
	//   hash_tbl->hash[c].pos[i][0], 
	//   hash_tbl->hash[c].pos[i][1]);      
	//   }
	//   printf("\n");
	// }

	free(mem_hash);
	free(hashes);
	free(indices);
	free(pos);
}

void SquidSplit(PARAMS* params, INDEX* index, HASH_TBL* hash_tbl) {
	memset(params->output_R1, 0x00, PATH);
	strcpy(params->output_R1, params->basename);
	strcat(params->output_R1, "_R1.fastq");
	memset(params->output_R2, 0x00, PATH);
	strcpy(params->output_R2, params->basename);
	strcat(params->output_R2, "_R2.fastq");
	memset(params->output_BED, 0x00, PATH);
	strcpy(params->output_BED, params->basename);
	strcat(params->output_BED, ".bed");
	switch (params->lib)
	{
	case -1:
	{
		exit(EXIT_FAILURE);
	}break;
	case 0:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibISFSearch);
	}break;
	case 1:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibISRSearch);
	}break;
	case 2:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibIUSearch);
	}break;
	case 3:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibOSFSearch);
	}break;
	case 4:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibOSRSearch);
	}break;
	case 5:
	{
		if (params->input_R1 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		MultiThreadManager(params, index, hash_tbl, LibOUSearch);
	}break;
	case 6:
	{
		if (params->input_R1 == NULL && params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL)
		{
			/* R1 was specified */
			MultiThreadManager(params, index, hash_tbl, LibSFSearchR1);
		}
		if (params->input_R1 == NULL)
		{
			/* R2 was specified */
			MultiThreadManager(params, index, hash_tbl, LibSFSearchR2);
		}
	}break;
	case 7:
	{
		if (params->input_R1 == NULL && params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL)
		{
			/* R1 was specified */
			MultiThreadManager(params, index, hash_tbl, LibSRSearchR1);
		}
		if (params->input_R1 == NULL)
		{
			/* R2 was specified */
			MultiThreadManager(params, index, hash_tbl, LibSRSearchR2);
		}
	}break;
	case 8:
	{
		if (params->input_R1 == NULL && params->input_R2 == NULL) INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		if (params->input_R2 == NULL)
		{
			/* R1 was specified */
			MultiThreadManager(params, index, hash_tbl, LibUSearchR1);
		}
		if (params->input_R1 == NULL)
		{
			/* R2 was specified */
			MultiThreadManager(params, index, hash_tbl, LibUSearchR2);
		}
	}break;
	default:
		INDEPENDENT_MSG('E', "The library you have selected does not match input files. Check your command line carefully.");
		break;
	}
}

void PrintParams(PARAMS* params) {
	if (VERBOSE) {
		fprintf(stderr, "\nFollowing are your session's parameters:\n");
		fprintf(stderr, "  Database: %s\n", params->db);
		if (params->input_R1) {
			fprintf(stderr, "  R1 file: %s\n", params->input_R1);
		}
		if (params->input_R2) {
			fprintf(stderr, "  R2 file: %s\n", params->input_R2);
		}
		fprintf(stderr, "  Output basename: %s\n", params->basename);
		fprintf(stderr, "  Lib: %s\n", LIB_TYPE[params->lib]);
		if(DIFF){
			fprintf(stderr, "  --diff option: ON\n");
		}
		else{
			fprintf(stderr, "  --diff option: OFF\n");
		}
		if(!NO_DISJOIN){
			fprintf(stderr, "  --disjoin option: ON\n");
		}
		else{
			fprintf(stderr, "  --disjoin option: OFF\n");
		}
		if(IGNORE_N){
			fprintf(stderr, "  --ignore_N option: ON\n");
		}
		else{
			fprintf(stderr, "  --ignore_N option: OFF\n");
		}
		if (MASK_LOWER) {
			fprintf(stderr, "  --mask-lower filter: ON\n");
		}
		else {
			fprintf(stderr, "  --mask-lower filter: OFF\n");
		}
		if (BED_OUT && !DIFF) {
			fprintf(stderr, "  BED output: Enabled\n");
		}
		else {
			fprintf(stderr, "  BED output: Disabled\n");
		}
		if (FASTQ_OUT) {
			fprintf(stderr, "  FASTQ output: Enabled\n");
		}
		else {
			fprintf(stderr, "  FASTQ output: Disabled\n");
		}
		if(EVALS){
			fprintf(stderr, "  Num Evals: %d\n", EVALS);
		}
		else{
			fprintf(stderr, "  -e option set to 0\n");
		}
		fprintf(stderr, "  Kmer size: %d\n", K);
		fprintf(stderr, "  Mismatches: %3.2f%%\n", (double)MISMATCH);
		fprintf(stderr, "  Step size: %d\n", STEP);
		if (AVAIL_THREADS == 1) {
			fprintf(stderr, "  Working on single thread\n");
		}
		else {
			fprintf(stderr, "  Using %d threads\n", AVAIL_THREADS);
		}
		fprintf(stderr, "\n");
	}
}

int CommandLine(PARAMS* params, unsigned int argc, char** argv) {
	if (argc == 1) return 1;
	if (argc == 2 && strcmp(argv[1], "-h") == 0) print_help();
	if (argc == 2 && strcmp(argv[1], "--help") == 0) print_help();

	/* Memory Initialization */
	params->db = NULL;      // required
	params->input_R1 = NULL;
	params->input_R2 = NULL;
	params->output_R1[0] = '\0';
	params->output_R2[0] = '\0';
	params->output_BED[0] = '\0';
	params->basename[0] = '\0';  // NO USER INPUT
	params->quiet = 0;      // default is to be verbose
	params->lib = -1;
	int i;
	for (i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-i") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-i", "string");
				return 1;
			}
			params->db = argv[1 + i];
			++i;
		}
		else if (strcmp(argv[i], "-R1") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-R1", "string");
				return 1;
			}
			params->input_R1 = argv[1 + i];
			++i;
		}
		else if (strcmp(argv[i], "-R2") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-R2", "string");
				return 1;
			}
			params->input_R2 = argv[1 + i];
			++i;
		}
		else if (strcmp(argv[i], "-o") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-o", "string");
				return 1;
			}
			strcpy(params->basename, argv[1 + i]);
			++i;
		}
		else if (strcmp(argv[i], "--mask-lower") == 0) MASK_LOWER = 1;
		else if (strcmp(argv[i], "--no-bed") == 0) BED_OUT = 0;
		else if (strcmp(argv[i], "--no-fastq") == 0) FASTQ_OUT = 0;
		else if (strcmp(argv[i], "--disjoin") == 0) NO_DISJOIN = 0;
		else if (strcmp(argv[i], "--quiet") == 0) VERBOSE = 0;
		else if (strcmp(argv[i], "--diff") == 0) DIFF = 1;
		else if (strcmp(argv[i], "--ignore_N") == 0) IGNORE_N = 1;
		else if (strcmp(argv[i], "-l") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-l", "string");
				return 1;
			}
			if (strcmp(argv[1 + i], "ISF") == 0) {
				params->lib = 0;
			}
			else if (strcmp(argv[1 + i], "ISR") == 0) {
				params->lib = 1;
			}
			else if (strcmp(argv[1 + i], "IU") == 0) {
				params->lib = 2;
			}
			else if (strcmp(argv[1 + i], "OSF") == 0) {
				params->lib = 3;
			}
			else if (strcmp(argv[1 + i], "OSR") == 0) {
				params->lib = 4;
			}
			else if (strcmp(argv[1 + i], "OU") == 0) {
				params->lib = 5;
			}
			else if (strcmp(argv[1 + i], "SF") == 0) {
				params->lib = 6;
			}
			else if (strcmp(argv[1 + i], "SR") == 0) {
				params->lib = 7;
			}
			else if (strcmp(argv[1 + i], "U") == 0) {
				params->lib = 8;
			}
			else {
				REQUIRE_ARG_ERR("-l", "string");
				return 1;
			}
			++i;
		}
		else if (strcmp(argv[i], "-m") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-m", "integer in range [0, 100)");
				return 1;
			}
			if (!is_int(argv[i + 1])) {
				REQUIRE_ARG_ERR("-m", "integer in range [0, 100)");
				return 1;
			}
			char* str;
			long str_to_i;
			str_to_i = strtol(argv[i + 1], &str, 10);
			if (str_to_i > 99 || str_to_i < 0) return 1;
			MISMATCH = (int)str_to_i;
			++i;
		}
		else if (strcmp(argv[i], "-e") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-e", "unsigned integer");
				return 1;
			}
			if (!is_int(argv[i + 1])) {
				REQUIRE_ARG_ERR("-e", "unsigned integer");
				return 1;
			}
			char* str;
			long str_to_i;
			str_to_i = strtol(argv[i + 1], &str, 10);
			if (str_to_i < 0) return 1;
			EVALS = (int)str_to_i;
			++i;
		}
		else if (strcmp(argv[i], "-s") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-s", "integer in range [1, L)");
				return 1;
			}
			if (!is_int(argv[i + 1])) {
				REQUIRE_ARG_ERR("-s", "integer in range [1, L)");
				return 1;
			}
			char* str;
			long str_to_i;
			str_to_i = strtol(argv[i + 1], &str, 10);
			if (str_to_i < 1) return 1;
			STEP = (int)str_to_i;
			++i;
		}
		else if (strcmp(argv[i], "-t") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-t", "integer");
				return 1;
			}
			if (!is_int(argv[i + 1])) {
				REQUIRE_ARG_ERR("-t", "integer");
				return 1;
			}
			char* str;
			long str_to_i;
			str_to_i = strtol(argv[i + 1], &str, 10);
			if (str_to_i < 0) return 1;
			AVAIL_THREADS = (int)str_to_i;
			++i;
		}
		else if (strcmp(argv[i], "-k") == 0) {
			if ((1 + i) >= argc) {
				REQUIRE_ARG_ERR("-k", "integer (9, 11, 13 or 15)");
				return 1;
			}
			if (!is_int(argv[i + 1])) {
				REQUIRE_ARG_ERR("-k", "integer (9, 11, 13, 15)");
				return 1;
			}
			char* str;
			long str_to_i;
			str_to_i = strtol(argv[i + 1], &str, 10);
			if (str_to_i != 9 &&
				str_to_i != 11 &&
				str_to_i != 13 &&
				str_to_i != 15) return 1;
			K = (int)str_to_i;
			++i;
		}
	}

	if (params->db == NULL) return 1;
	if (params->lib == -1) return 1;
	if (params->basename[0] == '\0') return 1;
	if (!FASTQ_OUT && !BED_OUT) {
		INDEPENDENT_MSG('W', "No output will be produced because \"--no-fastq\" and \"--no-bed\" flags were set to true");
	}
	if(EVALS && !NO_DISJOIN){
		INDEPENDENT_MSG('W', "\"--disjoin\" and \"-e\" flags are mutually exclusive. \"--disjoin\" option forced off");
		NO_DISJOIN = 1;
		// if more than 0 evals are requested, two different functions are loaded instead of default
		inward_searchfunc = eval_inward_s;
		outward_searchfunc = eval_outward_s;
	}
	if(EVALS && NO_DISJOIN){
		// if more than 0 evals are requested, two different functions are loaded instead of default
		inward_searchfunc = eval_inward_s;
		outward_searchfunc = eval_outward_s;				
	}
	if(!EVALS && NO_DISJOIN){
		// if more than 0 evals are requested, two different functions are loaded instead of default
		inward_searchfunc = no_disjoin_inward_s;
		outward_searchfunc = no_disjoin_outward_s;				
	}
	if(!EVALS && !NO_DISJOIN){
		// if more than 0 evals are requested, two different functions are loaded instead of default
		inward_searchfunc = disjoin_inward_s;
		outward_searchfunc = no_disjoin_outward_s;				
	}
	if(DIFF && BED_OUT){
		INDEPENDENT_MSG('W', "No BED output will be produced because \"--diff\" and BED output are both enabled");
	}
	switch (K)
	{
	case 9:
	{
		hashfunc = calc_hash_9;
	}break;
	case 11:
	{
		hashfunc = calc_hash_11;
	}break;
	case 13:
	{
		hashfunc = calc_hash_13;
	}break;
	case 15:
	{
		hashfunc = calc_hash_15;
	}break;
	default:
	{
		hashfunc = calc_hash_15;
	}break;
	}
	return 0;
}

int main(int argc, char** argv) {
	PARAMS params;
	int status = CommandLine(&params, argc, argv);
	if (status == 1) print_usage();
	INDEX index;
	HASH_TBL hash_tbl;
	PrintSizes();
	ReadDatabase(&params, &index, &hash_tbl);
	PrintParams(&params);
	SquidSplit(&params, &index, &hash_tbl);

	int i, j = 0;
	for (i = 0; i < index.size; ++i)
	{
		free(index.fasta[i].sequence);
	}
	free(index.fasta);
	for (i = 0; i < hash_tbl.size; ++i)
	{
		for (j = 0; j < hash_tbl.hash[i].size; ++j)
		{
			free(hash_tbl.hash[i].pos[j]);
		}
		free(hash_tbl.hash[i].pos);
	}
	free(hash_tbl.hash);
	exit(EXIT_SUCCESS);
	// valgrind --leak-check=yes --track-origins=yes
}
