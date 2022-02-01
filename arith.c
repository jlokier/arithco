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

#include "arith.h"

static void model_normalize(Model *model)
{
	if (USE_TRACE)
		fprintf(stderr, "Note: Updating symbol frequency model\n");

	count_t total = model->total - model->countdown;
	if (total == 0) {
		for (int symbol = 0; symbol < 256; symbol++)
			total += model->count[symbol];
	}
	model->countdown = (total < 128 ? 1 : total < 1024 ? 4 : total < 8192 ? 16 : 65536);
	model->total = total + model->countdown;

	/*
	 * Ensure `accum_scaled` and `accum_normalized` can't overflow 64-bit
	 * arithmetic with very large symbol counts and `FQ_BITS`.
	 */
	const int num_symbols = 256;
	int count_shift = 0;
	if (FQ_ONE > num_symbols) {
		const uint64_t divide_by =
			((uint64_t)(FQ_ONE - num_symbols) << 1) + 1;
		const uint64_t count_overflow_threshold =
			(((uint64_t)1 << 63) / divide_by) << 2;
		while ((total >> count_shift) >= count_overflow_threshold)
			count_shift++;
	}

	/*
	 * All symbols with a non-zero probabilities, no matter how small, must
	 * have normalised lookup range at least 1 otherwise they cannot be
	 * encoded.
	 *
	 * In forward-adaptive models, where probabilities are learned as the
	 * data stream progresses, all symbols may occur.  This is accounted
	 * for by starting all symbols at count 1, and they always have
	 * non-zero probabilities.  In backward-adaptive models where counts
	 * are encoded at the start and reduce to zero as the data stream
	 * progress, it's ok for symbols with zero count to have zero lookup
	 * range and be unable to be encoded.  Allowing that results in more
	 * efficient compression.  In static models the same applies, lookup
	 * range for symbols that don't occur can be zero and it's more
	 * efficient compression to allow that.
	 *
	 * The total of `accum_normalized[]` deltas must be less than or equal
	 * to `FQ_ONE`, ideally equal for best compression.  To make it equal
	 * we distribute the rounding errors by calculating `accum_normalized`
	 * from the cumulative count, rather than calculating each interval
	 * separately from its count.  Adding `accum_non_zero` to each
	 * `accum_normalized` ensures each lookup interval with non-zero
	 * probability has size at least 1.
	 */
	uint32_t scale_multiplier = FQ_ONE - num_symbols;
	if (USE_STATIC_MODEL && 1) {
		for (int symbol = 0; symbol < num_symbols; symbol++) {
			if (model->count[symbol] == 0)
				scale_multiplier++;
		}
	}

	bool for_decompression = model->for_decompression;
	uint32_t accum_non_zero = 0;
	uint64_t accum_scaled = 0;
	uint32_t accum_normalized = 0, j = 0;
	model->accum_normalized[0] = 0;
	for (int symbol = 0; symbol < num_symbols; symbol++) {
		count_t count = model->count[symbol];
		if (!USE_STATIC_MODEL || count != 0) {
			count >>= count_shift;
			accum_scaled += (uint64_t)count * scale_multiplier;
			accum_normalized = (accum_scaled + (total >> 1)) / total;
			accum_non_zero += 1;
			accum_normalized += accum_non_zero;
#if USE_TO_SYMBOL_MAP
			if (for_decompression) {
				for (; j < accum_normalized; j++)
					model->to_symbol[j] = symbol;
			}
#endif
		}
		model->accum_normalized[symbol+1] = accum_normalized;
	}

	assert(accum_normalized == FQ_ONE);
}

static void model_init(Model *model, StaticModel *static_model,
		       bool for_decompression)
{
	for (int symbol = 0; symbol < 256; symbol++) {
		model->count[symbol] =
			static_model ? static_model->count[symbol] : 1;
	}
	model->total = 0;
	model->countdown = 0;
	model->for_decompression = for_decompression;
	model_normalize(model);
}

static inline void model_update(Model *model, byte b)
{
	model->count[b]++;
	if (--model->countdown == 0)
		model_normalize(model);
}

static void static_model_init(StaticModel *static_model)
{
	for (int symbol = 0; symbol < 256; symbol++)
		static_model->count[symbol] = 0;
}

static inline void static_model_update(StaticModel *static_model, byte b)
{
	static_model->count[b]++;
}

/*
 * The `map` functions map between x/y number space and lookup/N number space.
 * All of x, y, lookup and N are integers.  When mapping x/y -> lookup/N, the
 * fractions are values in the semi-open interval [0,1), and when mapping
 * lookup/N -> x/y, they are values in the closed interval [0,1].
 *
 * Therefore depending on the map direction, x has range [0,y-1] or [0,y], and
 * lookup has range [0,N-1] or [0,N].
 *
 * In these maps, y has range [`FP_ONE`,`2*FP_ONE-1`], and N is `FP_ONE`.
 *
 * In standard arithmetic coding or range coding, the maps are linear scaling,
 * with some rounding errors in any practical implementation.  Linear mode is
 * selected by default (when `USE_NON_LINEAR` is not set).  In linear mode,
 * lookup/N to x/y uses a multiply.  x/y to lookup/N uses a divide unless
 * `USE_RECIP_TABLE` is set, in which case it uses a carefully designed lookup
 * table and multiply.
 *
 * Meaning of x/y and lookup/N
 * ---------------------------
 *
 * At each new symbol, x/y track bits of the compressed bitstream.  y is the
 * size of the current arithmetic code interval, left-shifted by enough bits to
 * bring y into the range stated above.  x is bits of the compressed bitstream,
 * left-shifted the same amount as y, minus the low-bound of the current
 * arithmetic code interval.
 *
 * The interval shrinks as symbols are processed, and its bounds after
 * left-shifting are rounded so that x and y are always integers.  y is a
 * function of only previous symbols.  In the view of arithmetic coding as a
 * variable-length bit encoding that supports fractions-of-a-bit, y represents
 * the fractional part of the position in the bitstream and x represents a
 * window over the bitstream.
 *
 * When decompressing, lookup has that name because it is used to look up the
 * next symbol.  It is x/y mapped to lookup/N in the fixed interval [0,N-1],
 * then lookup selects the next symbol using the symbol frequency model.  In
 * reverse lookup/N is mapped to x/y in both compression and decompression,
 * with lookup is the smallest value that would select the current symbol in
 * the frequency model.
 *
 * Constant or per-model N and model interleaving
 * ----------------------------------------------
 *
 * N is constant from one symbol to the next, unlike y, to make lookup more
 * efficient than using x and y directly in the frequency model.  It is
 * possible to switch between different models in the stream, in which case
 * each model may use a different N.  This is useful when interleaving symbols
 * from multiple dictionaries in a kind of symbol syntax.
 *
 * Rare flags
 * ----------
 *
 * Rare flags are mentioned here because they are an important case of
 * per-model N and model interleaving.  Almost zero bits are required on
 * average to store "rare flags", boolean flags with probability very close to
 * 0 or 1.  Interleaving symbols with rare flags allows flags to be used in
 * many places with almost no space overhead, adding a lot of syntax
 * flexibility.
 *
 * Non-linear option
 * -----------------
 *
 * It turns out any non-linear, monotonic map with sufficiently separated
 * values works, with some loss of compression efficiency.  If `USE_NON_LINEAR`
 * is set, the maps use an alternative map which avoids multiplication and
 * division entirely.
 *
 * This reduces compression ratio by about 1%, in exchange for avoiding
 * multiplications and divisions.  On 64-bit x86 the speed improvement is
 * small.  The speed improvement may be significant on older and low-powered
 * hardware, especially when multiplication is done in software.
 *
 * After developing this we found an excellent paper by Barney Stratford (2006)
 * which contains the same technique:
 *
 *   [Stratford2006] "A Faster Arithmetic Coder", Barney Stratford, Oxford
 *   University Computing Laboratory, 2006.
 *   http://www.cs.ox.ac.uk/research/pdt/ap/papers/fast.pdf
 *
 * That well-written paper is worth reading for showing "Rule of Floors" and
 * way of presenting clear and concise proofs of the linear arithmetic rounding
 * used in this code, as well as the non-linear method enabled by
 * `USE_NON_LINEAR`.  (It doesn't go into the `USE_RECIP_TABLE` method, though.
 * For that, see comment at `recip_table_init`.)
 *
 * However the ~1% compression loss claim omits to mention: The non-linear map
 * does not treat all symbols uniformly.  It is biased towards symbols with
 * lower lookup values, so the compression ratio depends on the order of the
 * symbol dictionary.
 *
 * This code inverts the kink shape compared with [Stratford2006], because the
 * compression ratio is slightly better when lower-valued symbols are more
 * common, as in for example source code (though context modelling or LZ
 * compression should really be used with that).  This code generally
 * encourages placing more common symbols at the lower end of the symbol range.
 *
 * Other non-monotonic maps are possible and some decompress faster, but
 * compression is more complicated, and for some maps takes exponential time.
 * "Asymmetric Numeral Systems"
 * <https://en.wikipedia.org/wiki/Asymmetric_numeral_systems> is equivalent to
 * a type of non-monotonic map here, but compression in the forward direction
 * is too expensive, so ANS compression is done in reverse.
 *
 */
static inline lookup_type map_xy_to_lookup(xy_type x, xy_type y)
{
	lookup_type lookup;
	if (0 && USE_NON_LINEAR) {
		/*
		 * Non-linear, monotonic map, centred kink, hoping
		 * out to be worse compression as well.
		 */
		xy_type excess = y - FP_ONE;
		xy_type t1 = FP_ONE >> 1, t2a = t1 + (excess << 1);
		lookup = (x <= t1 ? x : x >= t2a ? x - excess : t1 + ((x - t1) >> 1));
	} else if (USE_NON_LINEAR) {
		/* Non-linear, monotonic map . */
		xy_type excess = y - FP_ONE, t1 = (FP_ONE << 1) - y;
		if (1) {
			/* Inverted kink, slightly better compression, slower. */
			lookup = (x <= t1 ? x : t1 + ((x - t1) >> 1));
		} else {
			/* Kink same as in [Stratford2006]. */
			lookup = (x >= (excess << 1)) ? x - excess : x >> 1;
		}
	} else if (USE_RECIP_TABLE) {
		/* See `recip_table_init` for more about this. */
		xy_type pn = (xy_type)(y - x - 1);
		uint64_t pm = recip_table[y - FP_ONE];
		lookup = ((uint64_t)pm * pn) >> (2 * FP_BITS + 2);
		lookup ^= (FP_ONE - 1);
		if (USE_SLOW_CHECKS) {
			uint32_t div = (((uint64_t)x << FP_BITS) + (FP_ONE - 1)) / y;
			if (lookup != div) {
				fprintf(stderr, "Note: x %lu y %lu FP_BITS %d -> lookup %lu != div %lu\n",
					(unsigned long)x, (unsigned long)y, FP_BITS,
					(unsigned long)lookup, (unsigned long)div);
				abort();
			}
		}
	} else if (FP_BITS <= 15) {
		/* Division method, 32-bit optimisation. */
		lookup = (((uint32_t)x << FP_BITS) + (FP_ONE - 1)) / y;
	} else {
		/* Division method, 64-bit. */
		lookup = (((uint64_t)x << FP_BITS) + (FP_ONE - 1)) / y;
	}
	return lookup;
}

static inline xy_type map_lookup_to_xy(xy_type y, lookup_type lookup)
{
	xy_type x;
	if (0 && USE_NON_LINEAR) {
		/* Non-linear, monotonic map. */
		xy_type excess = y - FP_ONE, t1 = FP_ONE >> 1, t2 = t1 + excess;
		x = lookup + (lookup <= t1 ? 0 : lookup >= t2 ? excess : lookup - t1);
	} else if (USE_NON_LINEAR) {
		/* Non-linear, monotonic map. */
		xy_type excess = y - FP_ONE, t1 = (FP_ONE << 1) - y;
		if (1) {
			/* Inverted kink, slightly better compression, slower. */
			x = lookup + (lookup <= t1 ? 0 : lookup - t1);
		} else {
			/* Kink same as in [STRATFORD2006]. */
			x = lookup + (lookup < excess ? lookup : excess);
		}
	} else if (FP_BITS <= 15) {
		/* 32-bit optimisation. */
		x = ((uint32_t)lookup * y) >> FP_BITS;
	} else {
		/* 64-bit multiplication. */
		x = ((uint64_t)lookup * y) >> FP_BITS;
	}
	return x;
}

/*
 * Initialise `recip_table[]`.
 *
 * The division algorithm uses the unsigned division method in Figure 4.1 of:
 *
 *   [PLDI1994] "Division by Invariant Integers using Multiplication" by
 *   Torbjörn Granlund and Peter L. Montgomery, PLDI 1994.
 *
 * These values and the reciprocal-multiplication step are equivalent to a
 * pseudo-rounding-up version of fixed-point division x/y for all integers x in
 * the range [0..y-1], all integers y in the range [`FP_ONE`..`2*FP_ONE-1`],
 * and treating the denominator y as a fixed-point representation of a number
 * with `FP_BITS` fractional bits.
 *
 * The fixed-point aspect is handled by really dividing `(x << FP_BITS) / y`.
 * The shift comes from the fixed-point precision of the denominator only.  The
 *
 * The pseudo-rounding-up aspect adds an offset to the numerator as if shifting
 * in all ones, so `((x << FP_BITS) + ((1 << FP_BITS) - 1)) / y`.  It's not
 * exactly the same as rounding-up the fractional result of division, which
 * would add `y - 1` instead.  This offset is to compensate for the
 * rounding-down quantization in the `map_lookup_to_xy` multiplication,
 * which maps to new `x`, right-shifting the product by `FP_BITS`.
 *
 * The Granlund-Montgomery unsigned division method doesn't include our
 * fixed-point or pseudo-rounding-up.  If the adjusted numerator is used
 * directly with that method, it works but requires more bits in some
 * intermediate values and more steps.  By negating the input `x`, adding
 * appropriate offsets, selecting a constant shift `l`, folding the `n` offset
 * in the method's final step `(t + ((n - t) >> 1)) >> (l - 1)` into `pm2`, and
 * controlling the range of values so they can't overflow, the required
 * intermediate bits and number of steps are reduced.
 *
 * This algorithm works up to `FP_BITS` <= 21 using some 64-bit unsigned
 * intermediates.  At higher `FP_BITS`, use division instead.  The division
 * works up to `FP_BITS` <= 31.
 */
static void recip_table_init(void)
{
#if USE_RECIP_TABLE && !USE_NON_LINEAR
	for (uint32_t y = FP_ONE; y <= 2 * FP_ONE - 1; y++) {
		int N = 2 * FP_BITS + 1; /* Bit-range of shifted numerator. */
		int pl = FP_BITS + 1;    /* Reciprocal calculation shift `l`. */
		uint32_t pd = y;         /* Invariant denominator. */
		uint64_t pm = (((((uint64_t)1 << pl) - pd) << N) / pd) + 1;
		uint64_t pm2 = ((uint64_t)1 << N) + (uint64_t)pm;
		recip_table[y - FP_ONE] = pm2;
	}
#endif

	/*
	 * If `FP_BITS` is > 2^14, sample a 2^14-size random subset of `y`,
	 * `l` and `x` to keep testing time down.
	 */
	long random_mask = FP_BITS <= 14 ? 0 : (1L << (FP_BITS - 14)) - 1;
	(void)random_mask; /* Prevent unused variable warning. */

#if USE_SLOW_CHECKS && USE_RECIP_TABLE && !USE_NON_LINEAR
	/*
	 * Verify the reciprocal-multiplication method, which is exactly equal
	 * to pseudo-rounding-up division of fixed-point `x/y` for all integers
	 * `x` in the range [`0`..`y-1`] and all integers `y` in the range
	 * [`FP_ONE`..`2*FP_ONE-1`], treating `y` as a fixed-point
	 * representation of a number with `FP_BITS` fractional bits.
	 */
	for (uint32_t y = FP_ONE; y <= 2 * FP_ONE - 1; y++) {
		if (random_mask != 0 && (random_mask & random()) != 0)
			continue;
		for (uint32_t x = 0; x < y; x++) {
			if (random_mask != 0 && (random_mask & random()) != 0)
				continue;
			uint32_t div = (((uint64_t)x << FP_BITS) + (FP_ONE - 1)) / y;
			uint64_t pm = recip_table[y - FP_ONE];
			uint32_t pn = (uint32_t)(y - x - 1);
			uint32_t pt = ((uint64_t)pm * pn) >> (2 * FP_BITS + 2);
			uint32_t sh = (FP_ONE - 1) ^ (uint32_t)pt;
			if (sh != div) {
				fprintf(stderr, "Note: Recip mismatch, y 0x%08lx x 0x%08lx -> pm 0x%08llx pn 0x%08lx pm*pn 0x%016llx sh 0x%08lx div 0x%08lx\n",
					(unsigned long)y, (unsigned long)x,
					(unsigned long long)pm, (unsigned long)pn,
					(unsigned long long)(pm*pn),
					(unsigned long)sh, (unsigned long)div);
			}
		}
	}
#endif

#if 0
	/*
	 * Verify an old reciprocal-multiplication method, which successfully
	 * reverses the mapping `l` -> `x` to `x` -> `sh` for all `x` that
	 * arise from the forward mapping from all integers `l` in the range
	 * [`1`...`FP_ONE`] and all integers `y` in the range
	 * [`FP_ONE`..`2*FP_ONE-1`].
	 *
	 * This proved to be not good enough because the decompressor needs the
	 * reverse map `x` -> `sh` to return pseudo-rounding-up division for
	 * all integers `x` in the range 0 to `y-1`.  That is, the in-between
	 * values of `x` too, not just those from the mapping `l` -> `x`.
	 */
	for (uint32_t y = FP_ONE; y <= 2 * FP_ONE - 1; y++) {
		if (random_mask != 0 && (random_mask & random()) != 0)
		    continue;
		/* This is the number of bits required. */
		int RECIP_SHIFT = FP_BITS * 2 - 1;
		uint64_t recip = ((uint64_t)FP_ONE << RECIP_SHIFT) / y;
		for (uint32_t l = 1; l <= FP_ONE; l++) {
			if (random_mask != 0 && (random_mask & random()) != 0)
				continue;
			uint32_t x = ((uint64_t)l * y) >> FP_BITS;
			uint32_t sh = (((uint64_t)x * recip) + (((uint64_t)1 << RECIP_SHIFT) - 1)) >> RECIP_SHIFT;
			uint32_t div = (((uint64_t)x << FP_BITS) + (FP_ONE - 1)) / y;
			if (sh != div) {
				fprintf(stderr, "Note: Recip mismatch, y 0x%08lx x 0x%08lx -> recip 0x%08lx sh 0x%08lx div 0x%08lx\n",
					(unsigned long)y, (unsigned long)x,
					(unsigned long)recip, (unsigned long)sh,
					(unsigned long)div);
			}
			if (sh != l) {
				fprintf(stderr, "Note: Recip mismatch, y 0x%08lx l 0x%08lx -> x 0x%08lx -> recip 0x%08lx sh 0x%08lx div 0x%08lx\n",
					(unsigned long)y, (unsigned long)l,
					(unsigned long)x, (unsigned long)recip,
					(unsigned long)sh, (unsigned long)div);
			}
		}
	}
#endif
}

static void compress_init(Compress *c, FILE *output, StaticModel *static_model)
{
	model_init(&c->model, static_model, false);
	c->output = output;
	c->x = 0;
	c->y = 1;
	c->cbit_count = 0;
	c->carry_out = 0;
	c->carry_digit = 0;
	c->compressed_bytes = 0;
}

static inline void compress_put_byte(Compress *c, byte b)
{
	putc(b, c->output);
	c->compressed_bytes++;
}

static inline void compress_put_word32(Compress *c, uint32_t word)
{
	compress_put_byte(c, (byte)(word >> 24));
	compress_put_byte(c, (byte)(word >> 16));
	compress_put_byte(c, (byte)(word >> 8));
	compress_put_byte(c, (byte)word);
}

static void decompress_init(Decompress *d, FILE *input, StaticModel *static_model)
{
	model_init(&d->model, static_model, true);
	d->input = input;
	d->n_bits_in = 0;
	d->bits_in = 0;
	d->x = 0;
	d->y = 1;
}

static inline byte decompress_get_byte(Decompress *d)
{
	int ch = getc(d->input);
	if (ch == EOF)
		ch = 0;
	return (byte)ch;
}

static inline uint32_t decompress_get_word32(Decompress *d)
{
	uint32_t word = 0;
	word = (word << 8) | decompress_get_byte(d);
	word = (word << 8) | decompress_get_byte(d);
	word = (word << 8) | decompress_get_byte(d);
	word = (word << 8) | decompress_get_byte(d);
	return word;
}

static inline void compress_put_committed_digits(Compress *c)
{
	int cbit = c->cbit_count;
	int keep_bits = count_non_zero_value_width(c->y - 1);

	if (cbit - DIGIT_BITS < keep_bits)
		return;

	do {
		cbit -= DIGIT_BITS;
		xy_extended_type edigit = c->x >> cbit;
		c->x -= edigit << cbit;

		if (edigit >= ((xy_extended_type)1 << DIGIT_BITS)) {
			if (USE_SLOW_CHECKS) {
				if (c->carry_out <= 0)
					abort();
				if (edigit >= (xy_extended_type)1 << (DIGIT_BITS+1))
					abort();
			}
			compress_put_digit(c, c->carry_digit + 1);
			while (--c->carry_out != 0)
				compress_put_digit(c, 0);
			edigit &= DIGIT_HIGH;
		}
		digit_type digit = edigit;

		if (c->carry_out != 0) {
			if (digit == DIGIT_HIGH) {
				c->carry_out++;
				continue;
			}
			compress_put_digit(c, c->carry_digit);
			while (--c->carry_out != 0)
				compress_put_digit(c, DIGIT_HIGH);
		}

		c->carry_digit = digit;
		c->carry_out = 1;
	} while (cbit - DIGIT_BITS >= keep_bits);

	c->cbit_count = cbit;
}

static void compress_finish(Compress *c)
{
	/*
	 * Some bits of the final symbol can be skipped.  Calculate the minimum
	 * number of final bits which identify the correct symbol, assuming the
	 * decompressor pads the bitstream with unlimited zeros.  Note, `cy`
	 * cannot be right-shifted during the loop as its low-order bits still
	 * affect the decision.
	 */
	xy_extended_type cx = c->x, cy = c->y, bit = 1;
	int cbit = c->cbit_count, drop_low_bits = 0;
	while (drop_low_bits < cbit) {
		if ((cx & bit) != 0) {
			if (cx + bit >= cy - 1 + cx)
				break;
			cx += bit;
			cy -= bit;
		}
		bit <<= 1;
		drop_low_bits++;
	}
	/* This shift conveniently forces output padding bits to zero. */
	cx >>= drop_low_bits;
	cbit -= drop_low_bits;

	/* Ensure all pending bits are output to whole final bytes/words. */
	int pad_bits = (-cbit & (DIGIT_BITS - 1));
	c->x = cx << (pad_bits + 1);
	c->cbit_count = cbit + (pad_bits + 1);
	/* Subtle: `count_leading_zeros(y-1)` is undefined for y == 1. */
	c->y = 2;
	compress_put_committed_digits(c);

	if (c->carry_out != 0) {
		compress_put_digit(c, c->carry_digit);
		while (--c->carry_out != 0)
			compress_put_digit(c, DIGIT_HIGH);
	}
}

static inline void compress_symbol(Compress *c, byte symbol)
{
	xy_extended_type cx = c->x;
	xy_type cy = c->y;

	int xy_shift = (FP_BITS + 1) - count_non_zero_value_width(cy);
	cx <<= xy_shift;
	cy <<= xy_shift;
	c->cbit_count += xy_shift;

	lookup_type base = c->model.accum_normalized[(unsigned)symbol];
        lookup_type top = c->model.accum_normalized[(unsigned)symbol+1];

	if (USE_SLOW_CHECKS && (top <= base || top > FQ_ONE))
		abort();

	xy_type xy_base = map_lookup_to_xy(cy, base << (FP_BITS - FQ_BITS));
	xy_type xy_top = map_lookup_to_xy(cy, top << (FP_BITS - FQ_BITS));

	c->x = cx + xy_base;
	c->y = xy_top - xy_base;

	if (USE_TRACE) {
		int symbol_digits       = 3;
		int lookup_digits       = (FQ_BITS + 1 + 3) >> 2;
		int xy_digits           = (FP_BITS + 1 + 3) >> 2;
		fprintf(stderr, "Note: cx/cy 0x%0*llx/0x%0*llx -> symbol %*u -> range [0x%0*lx,0x%0*lx)/0x%0*lx -> xy_range [0x%0*lx,0x%0*lx)/0x%0*lx -> cx/cy 0x%0*llx/0x%0*llx\n",
			xy_digits,     (unsigned long long)cx,
			xy_digits,     (unsigned long long)cy,
			symbol_digits, (unsigned)symbol,
			lookup_digits, (unsigned long)base,
			lookup_digits, (unsigned long)top,
			lookup_digits, (unsigned long)FQ_ONE,
			xy_digits,     (unsigned long)xy_base,
			xy_digits,     (unsigned long)xy_top,
			xy_digits,     (unsigned long)cy,
			xy_digits,     (unsigned long long)c->x,
			xy_digits,     (unsigned long long)c->y);
	}

	/*
	 * New interval must be at least 1, be contained in the previous
	 * interval, and be smaller to ensure progress.
	 */
	if (USE_SLOW_CHECKS && (xy_top <= xy_base || xy_top >= cy))
		abort();

	compress_put_committed_digits(c);
}

static inline unsigned decompress_lookup(Decompress *d, uint32_t lookup)
{
	unsigned symbol;
#if USE_TO_SYMBOL_MAP
	symbol = d->model.to_symbol[lookup];
#elif 1
	symbol = 0;
	if (lookup >= d->model.accum_normalized[symbol + 128])
		symbol += 128;
	if (lookup >= d->model.accum_normalized[symbol + 64])
		symbol += 64;
	if (lookup >= d->model.accum_normalized[symbol + 32])
		symbol += 32;
	if (lookup >= d->model.accum_normalized[symbol + 16])
		symbol += 16;
	if (lookup >= d->model.accum_normalized[symbol + 8])
		symbol += 8;
	if (lookup >= d->model.accum_normalized[symbol + 4])
		symbol += 4;
	if (lookup >= d->model.accum_normalized[symbol + 2])
		symbol += 2;
	if (lookup >= d->model.accum_normalized[symbol + 1])
		symbol += 1;
#else
	for (symbol = 0; symbol < 256; symbol++) {
		if (lookup < d->model.accum_normalized[symbol+1])
			break;
	}
#endif
	return symbol;
}

static inline byte decompress_symbol(Decompress *d)
{
	xy_extended_type edx = d->x, edy = d->y;
	if (edy < FP_ONE) {
		do {
			edx = (edx << DIGIT_BITS) | decompress_get_digit(d);
			edy <<= DIGIT_BITS;
		} while (DIGIT_BITS < FP_BITS && edy < FP_ONE);
	}
	int xy_ignore_bits =
		count_non_zero_value_width(edy) - (FP_BITS + 1);
	xy_type dy = edy >> xy_ignore_bits, dx = edx >> xy_ignore_bits;

	uint32_t lookup = map_xy_to_lookup(dx, dy);
	unsigned symbol = decompress_lookup(d, lookup >> (FP_BITS - FQ_BITS));
	uint32_t base = d->model.accum_normalized[symbol];
	uint32_t top = d->model.accum_normalized[symbol+1];

	if (USE_SLOW_CHECKS && (top <= base || top > FQ_ONE))
		abort();

	xy_type xy_base = map_lookup_to_xy(dy, base << (FP_BITS - FQ_BITS));
	xy_type xy_top = map_lookup_to_xy(dy, top << (FP_BITS - FQ_BITS));

	d->x = edx - ((xy_extended_type)xy_base << xy_ignore_bits);
	d->y = (xy_extended_type)(xy_top - xy_base) << xy_ignore_bits;

	if (USE_TRACE) {
		int symbol_digits       = 3;
		int lookup_digits       = (FQ_BITS + 1 + 3) >> 2;
		int xy_digits           = (FP_BITS + 1 + 3) >> 2;
		fprintf(stderr, "Note: dx/dy 0x%0*lx/0x%0*lx -> lookup 0x%0*lx/0x%0*lx -> symbol %*u -> range [0x%0*lx,0x%0*lx)/0x%0*lx -> xy_range [0x%0*lx,0x%0*lx)/0x%0*lx -> dx/dy 0x%0*lx/0x%0*lx\n",
			xy_digits,     (unsigned long)dx,
			xy_digits,     (unsigned long)dy,
			lookup_digits, (unsigned long)lookup,
			lookup_digits, (unsigned long)FQ_ONE,
			symbol_digits, (unsigned)symbol,
			lookup_digits, (unsigned long)base,
			lookup_digits, (unsigned long)top,
			lookup_digits, (unsigned long)FQ_ONE,
			xy_digits,     (unsigned long)xy_base,
			xy_digits,     (unsigned long)xy_top,
			xy_digits,     (unsigned long)dy,
			xy_digits,     (unsigned long)(dx - xy_base),
			xy_digits,     (unsigned long)(xy_top - xy_base));
	}

	/*
	 * New interval must contain `dx`, be contained in the previous
	 * interval, and be smaller to ensure progress.
	 */
	if (USE_SLOW_CHECKS
	    && (dx < xy_base || dx >= xy_top || xy_top - xy_base >= dy))
		abort();

	return (byte)symbol;
}
