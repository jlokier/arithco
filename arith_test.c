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

#include "arith.c"

static FILE *static_model_calculate(FILE *input, StaticModel *static_model)
{
	char *mem_data = 0;
	size_t mem_size = 0;
	FILE *mem_file = open_memstream(&mem_data, &mem_size);
	if (!mem_file)
		exit_perror("open_memstream");

	static_model_init(static_model);
	while (1) {
		int ch = getc(input);
		if (ch == EOF)
			break;
		putc((byte)ch, mem_file);
		byte b = (byte)ch ^ USE_BYTE_XOR;
		static_model_update(static_model, b);
	}
	fseeko(mem_file, 0, SEEK_SET);
	return mem_file;
}

static count_t compress(FILE *input, FILE *output, StaticModel *static_model,
			int64_t *time_out)
{
	Compress c;
	compress_init(&c, output, static_model);

	int64_t time_before = get_time();
	count_t count_bytes = 0;
	while (1) {
		int ch = getc(input);
		if (ch == EOF)
			break;
		count_bytes++;
		byte b = (byte)ch ^ USE_BYTE_XOR;
		compress_symbol(&c, b);
		model_update(&c.model, b);
	}
	compress_finish(&c);
	*time_out = get_time() - time_before;

	return count_bytes;
}

static void decompress(count_t count_bytes, FILE *input, FILE *output,
		       StaticModel *static_model, int64_t *time_out)
{
	Decompress d;
	decompress_init(&d, input, static_model);
	char *out = xmalloc(count_bytes);

	int64_t time_before = get_time();
	char *end = out;
	while (count_bytes > 0) {
		byte b = decompress_symbol(&d);
#if 0
		if (feof(input) && !ferror(input))
			break;
#endif
		*end++ = b ^ USE_BYTE_XOR;
		model_update(&d.model, b);
		count_bytes--;
	}
	*time_out = get_time() - time_before;

	if (output) {
		for (char *ptr = out; ptr < end; ptr++)
			putc(*ptr, output);
	}
	free(out);
}

int main(int argc, char **argv)
{
	FILE *input = stdin;

	if (argc > 1) {
		input = fopen(argv[1], "r");
		if (input == NULL)
			exit_perror("fopen(input)");
	}

	fprintf(stderr, "Note: Reading input into memory for timing accuracy\n");

	char *input_data = 0;
	size_t input_size = 0;
	FILE *input_mem_file = open_memstream(&input_data, &input_size);
	if (!input_mem_file)
		exit_perror("open_memstream(input_mem_file)");
	while (1) {
		int ch = getc(input);
		if (ch == EOF)
			break;
		putc((byte)ch, input_mem_file);
	}
	if (ferror(input))
		exit_perror("getc(input)");
	if (fflush(input_mem_file) != 0)
		exit_perror("fflush(input_mem_file)");
	if (fseeko(input_mem_file, 0, SEEK_SET) != 0)
		exit_perror("fseeko(input_mem_file, 0)");

	if (input != stdin)
		fclose(input);
	input = input_mem_file;

	char *compressed_data = 0;
	size_t compressed_size = 0;
	FILE *compressed_file = open_memstream(&compressed_data, &compressed_size);
	if (!compressed_file)
		exit_perror("open_memstream(compressed_file)");

#if USE_STATIC_MODEL
	fprintf(stderr, "Note: Starting static_model_calculate\n");
	StaticModel the_static_model, *static_model = &the_static_model;
	input = static_model_calculate(input, static_model);
#else
	StaticModel *static_model = NULL;
#endif

	fprintf(stderr, "Note: Starting recip_table_init\n");
	recip_table_init();

	int64_t time;
	fprintf(stderr, "Note: Starting compression\n");
	count_t count_bytes = compress(input, compressed_file, static_model, &time);
	fprintf(stderr, "Note: Compression time %.6f seconds, %.2f ns per byte\n",
		(double)time / 1e9, (double)time / count_bytes);

	count_t compressed_bytes = ftello(compressed_file);
	fprintf(stderr, "Note: Uncompressed %llu bytes, compressed %llu bytes, removed %.3f %%\n",
		(unsigned long long)count_bytes,
		(unsigned long long)compressed_bytes,
		(double)((int64_t)count_bytes - (int64_t)compressed_bytes) * 100 / count_bytes);

	if (fseeko(compressed_file, 0, SEEK_SET) != 0)
		exit_perror("fseeko(compressed_file, 0)");

	fprintf(stderr, "Note: Starting decompression\n");
	decompress(count_bytes, compressed_file, stdout, static_model, &time);
	fprintf(stderr, "Note: Decompression time %.6f seconds, %.2f ns per byte\n",
		(double)time / 1e9,
		(double)time / count_bytes);

	fclose(compressed_file);
	if (input != stdin)
		fclose(input);
	return EXIT_SUCCESS;
}
