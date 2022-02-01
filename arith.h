/*
 * R&D in tree-structure friendly arithmetic coding techniques.
 *
 * Copyright (C) 2021, 2022 Jamie Lokier
 *
 * This file is licensed under either of "MIT license" or "Apache License,
 * Version 2.0", at your option.  Links to each license respectively:
 * - <http://opensource.org/licenses/MIT>
 * - <http://www.apache.org/licenses/LICENSE-2.0>.
 *
 * This file is provided without any warranty.  Use at your own risk.  It is
 * intended that excerpts be used and changed in other programs, subject to the
 * terms of the one or both of the above licenses.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <assert.h>
#include <time.h>

/*
 * The number of bits of the fractional part of arithmetic code ranges.
 * A larger values results in slower operations but better compression ratio.
 */
#define FP_BITS 23

/*
 * Set this to 1 to use a static symbol frequency model made from a histogram
 * of all the input bytes, instead of an adaptive model that updates as
 * compression and decompression progress.
 *
 * In tests this slightly improved the compression ratio, maybe due to the
 * existence of zero-frequency symbols.  But improvement cannot be assumed in
 * general.  A static frequency model is sub-optimal compared with a reverse
 * adaptive model that starts with an accurate frequency model then decreases
 * frequencies as symbols are encountered and won't occur again.
 */
#ifndef USE_STATIC_MODEL
#define USE_STATIC_MODEL 0
#endif

/*
 * Set this to 1 to use a non-linear, monotonic map instead of linear scaling.
 * This reduces compression ratio by about 1%, in exchange for avoiding
 * multiplications and divisions.  On 64-bit x86 the speed improvement is
 * small.  The speed improvement may be significant on older and low-powered
 * hardware.
 */
#ifndef USE_NON_LINEAR
#define USE_NON_LINEAR 1
#endif

/*
 * Set this to 1 to use `recip_table[]` to replace division with a lookup table
 * and carefully rounded fixed-point multiplication.  See `recip_table_init`
 * for more about the algorithm.
 *
 * This is usually faster up to some `FP_BITS`, because most architectures have
 * slower division than multiplication, whether either function is in hardware
 * or not.  On some which have hardware multiplication only, division is much
 * slower.  The table lookup limits the useful size.
 */
#ifndef USE_RECIP_TABLE
#define USE_RECIP_TABLE 1
#endif

/*
 * Set this to 1 to use `to_symbol[]` table lookup for the next symbol when
 * decompressing, instead of binary searching the frequency table.  This is
 * usually faster up to some `FQ_BITS`.  Above that, the size of the table
 * slows it down.
 */
#ifndef USE_TO_SYMBOL_MAP
#define USE_TO_SYMBOL_MAP 1
#endif

/*
 * Set this to 1 to process the compressed stream a 32-bit word at a time,
 * instead of an 8-bit byte at a time.
 *
 * TODO: This speeds up compression but currently slows decompression a little.
 * Figure out why, as this should speed decompression as well.
 */
#ifndef USE_WORD_AT_A_TIME
#define USE_WORD_AT_A_TIME 1
#endif

/* Turn on detailed tracing of each compression/decompression step. */
#ifndef USE_TRACE
#define USE_TRACE 0
#endif

/* Turn on slow integrity checks to catch arithmetic issues. */
#ifndef USE_SLOW_CHECKS
#define USE_SLOW_CHECKS 0
#endif

/* XOR the uncompressed input and output bytes with this value for tests. */
#ifndef USE_BYTE_XOR
#define USE_BYTE_XOR 0
#endif

typedef uint8_t  byte;
typedef uint64_t count_t;
typedef uint32_t xy_type;
typedef uint32_t lookup_type;

/*
 * Word-at-a-time mode needs 64-bit x and y. Byte-at-a-time mode works with
 * 32-bit x and y if `FP_BITS` <= 23.  (Decompression works with `FP_BITS` <=
 * 24.)  Bit-at-a-time mode works up to `FP_BITS` <= 30, but we don't support
 * bit-at-a-time any more.
 */
#if USE_WORD_AT_A_TIME || FP_BITS > 23
typedef uint64_t xy_extended_type;
#else
typedef xy_type xy_extended_type;
#endif

/*
 * When `FP_BITS` is > 30 the algorithm values no longer fit in 32-bit types.
 */
#if FP_BITS > 30
#error This arithmetic coder breaks if FP_BITS > 30
#endif

/*
 * When `FP_BITS` is >= 21 the `recip_table[]` fixed-point calculation no
 * longer fits in 64-bit integers, but the division equivalent works fine.
 */
#if FP_BITS >= 21
#undef  USE_RECIP_TABLE
#define USE_RECIP_TABLE 0
#endif

/*
 * When `FP_BITS` is large `recip_table[]` uses a lot of memory and decoding
 * speed is reduced.
 *
 * The cut-off when it becomes slower than division depends (a lot) on the
 * architecture, but in my tests (late 2013 Haswell i7 @ 2.3GHz) the transition
 * happened at >= 16.
 *
 * Splitting the table into accurate reciprocal for upper bits and non-linear
 * or slightly different linear map for lower bits might change this.
 */
#if FP_BITS >= 16
#undef  USE_RECIP_TABLE
#define USE_RECIP_TABLE 0
#endif

/*
 * Allow `FQ_BITS` to be smaller than `FP_BITS`.  Increasing both of these
 * improves compression but the effect is more significant for `FP_BITS` at
 * larger values, and `FQ_BITS` determines the size of `to_symbol[]`, so it's
 * worth limiting `FQ_BITS`.
 */
#define FQ_BITS_MAX 15
#define FQ_BITS_MIN 8
#define FQ_BITS (FP_BITS <= FQ_BITS_MAX ? FP_BITS : FQ_BITS_MAX)
#if FQ_BITS < FQ_BITS_MIN
#error This arithmetic coder requires FQ_BITS >= FQ_BITS_MIN
#endif

/*
 * When `FQ_BITS` is large `to_symbol[]` uses a lot of memory and decoding
 * speed is reduced.
 *
 * The cut-off when it becomes slower than binary search depends on the
 * architecture, how often `to_symbol[]` is updated and the compressed data,
 * but in my tests (late 2013 Haswell i7 @ 2.3GHz) the transition happened at
 * >= 22.  A two-level table would change this.
 *
 * When switching between multiple models and/or using symbols with more than 8
 * bits, this threshold is likely to be smaller due to more memory being used.
 */
#if FQ_BITS >= 22
#undef  USE_TO_SYMBOL_MAP
#define USE_TO_SYMBOL_MAP 0
#endif

#define FP_ONE   ((uint32_t)1 << FP_BITS)
#define FQ_ONE   ((uint32_t)1 << FQ_BITS)

#if USE_RECIP_TABLE && !USE_NON_LINEAR
/*
 * At `FP_BITS <= 14` the entries fit into 32 bits.  At `FP_BITS == 15` all but
 * one entry fits and `recip_table[0]` has value 2^32 exactly.
 */
#if FP_BITS <= 14
static uint32_t recip_table[FP_ONE];
#else
static uint64_t recip_table[FP_ONE];
#endif
#else
/* Allow the code in unused branches to still be compiled. */
static byte recip_table[1];
#endif

#if USE_WORD_AT_A_TIME

#define DIGIT_BITS 32
#define DIGIT_HIGH ((digit_type)0xffffffffUL)
typedef uint32_t digit_type;
#define compress_put_digit compress_put_word32
#define decompress_get_digit decompress_get_word32

#else

#define DIGIT_BITS 8
#define DIGIT_HIGH ((byte)0xff)
typedef byte digit_type;
#define compress_put_digit compress_put_byte
#define decompress_get_digit decompress_get_byte

#endif

/*
 * API structures.
 */

typedef struct StaticModel {
	count_t count[256];
} StaticModel;

typedef struct Model {
	count_t count[256];
	int countdown;
	uint32_t accum_normalized[256+1];
#if USE_TO_SYMBOL_MAP
	byte to_symbol[FQ_ONE];
#endif
	count_t total;
	bool for_decompression;
} Model;

typedef struct Compress {
	Model model;
	FILE *output;
	xy_extended_type x;
	xy_type y;
	int cbit_count;
	count_t carry_out;
	digit_type carry_digit;
	count_t compressed_bytes;
} Compress;

typedef struct Decompress {
	Model model;
	FILE *input;
	byte bits_in;
	int n_bits_in;
	xy_extended_type x, y;
} Decompress;

/* This macro because GCC requires different builtins depending on the type. */
#define count_leading_zeros(n) (sizeof(n) > sizeof(unsigned int) \
				? __builtin_clzl(n) : __builtin_clz(n))
#define count_trailing_zeros(n) (sizeof(n) > sizeof(unsigned int) \
				 ? __builtin_ctzl(n) : __builtin_ctz(n))
/* Subtle: `count_leading_zeros(0)` is undefined. */
#define count_non_zero_value_width(n) (8 * sizeof(n) - count_leading_zeros(n))

static void
#ifdef __GNUC__
__attribute__((noinline, noreturn, cold))
#endif
exit_perror(const char *str)
{
	perror(str);
	exit(EXIT_FAILURE);
}

static void
#ifdef __GNUC__
__attribute__((noinline, noreturn, cold, format(printf, 1, 2)))
#endif
exit_printf(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	exit(EXIT_FAILURE);
}

int64_t get_time(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (1000000000 * ts.tv_sec) + ts.tv_nsec;
}

void * xmalloc(size_t size)
{
	void *ptr = malloc(size);
	if (!ptr && size)
		exit_perror("malloc");
	return ptr;
}
