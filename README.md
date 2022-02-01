# R&D in tree-structure friendly arithmetic coding techniques

## Purpose

This is an R&D repo.  The purpose of this **arithmetic code compression**
implementation is to test and explore support compressed-serialisation
primitives for **tree-structured compression** of **balanced search tree data
structures**.

It is not a library, at least not for now.  It's for taking code from the main
source file and using parts in bespoke encoding applications.

This repo contains a simple arithmetic coder (also called a range coder), and a
test and benchmarking program that uses it.  The coding primitives are work in
progress.

Some of the items shown below aren't implemented here but the code is
incorporated into works using them.  Some of those primitives will work their
way back here for general use.

## Key features

Key features (some are work in progress though):

- It is a toolkit of code primitives to crib from, not a library (at least for now).

- Designed to support multiple symbol frequency models interleaved.

- Symbol dictionaries form a syntax, representing stored data structures.

- Fast special encoding of equiprobable byte sequences for things like hashes.

- Final zero bits may be omitted from the compressed bitstream.  This is more
  significant than in regular bulk compression, because these final zeros occur
  at every _fragments_ in the encoding of a fast search tree.

- Static and adaptive models for testing the effects.

- "Rare flags" that take almost _zero bits_ of space in the common case when
  they are not set.  These are very close to 0 or 1 probability and higher
  fractional precision, to interleave metadata flags into data with low
  overhead.  They are an effective alternative to run-length encoding (RLE) in
  many applications, because in high precision arithmetic coding, runs of unset
  rare flags behave like counters, even when interleaved with other symbols.

- Encoding special patch slots and back-patching to insert "forward skips" into
  the stream, where the size of the skip is not known when the slot is first
  encoded.  These are key to supporting fast, balanced search trees and fast
  updates.

- Designed to support "multi-dimensional prefix coding" of related items in a
  stream.  This can be thought of is similar to dictionary compression like LZ
  methods, but supports fast search queries and updates, numerical deltas as
  well as string prefixes, and separate deltas for multiple columns in a tuple.

The use of arithmetic/range encoding for tree data structures that are full of
pointers may be unusual.

There are a number of compile-time options, commented in the code.  They have
been used to explore the design space, test on different architectures, and get
the corner cases of the algorithm working just right on each one.

A number of academic references were studied in the making of this code.  Two
are mentioned in the code but it deserves a more complete list.

## Comparison to ad-hoc, hand-crafted dense byte-codes

Two key benefits over ad-hoc encodings of structured data where the format has
to be optimised for _large-scale compactness_, such as a hand-crafted byte-code
with carefully packed flags, are:

- Arithmetic encoding is _tunable_ with less effort.
- _Simpler code_.  The code which uses the symbol stream actually ends up
  simpler and clearer.
- Adding and removing values that are required, such as flags, pointers, extra
  fields and so on is very easy with automatic low or appropriate overhead
  without requiring the rest of the encoding to be modified to keep it compact.

Tuning can take a lot of time to get right with hand-crafted codes, and when
the data changes, making small changes to the hand-crafted encoding that
remains fast and compact oftens turns into a combinatorial explosion.
Experiments have to be conducted too find which combinations work best, and in
practice there isn't time to run every experiment so only a few combinations
get used.

(For example, changing the code to allocate a bit here, remove it from there,
combine 10 sparse flags into 5 dense ones using a lookup table, add extra
"exceptional" codes to improve a common case, etc.)

It is also difficult or impossible to make such encodings adaptive at run time,
depending how they are implemented.

The arithmetic coding approach has approximately _continuous tuning knobs_,
over a space where hand-crafted byte-code tends to turn into a combinatorial
explosion of trade-offs.  Even the case of adding or removing "exceptional"
byte-codes, or mapping sets of sparse flags into a denser set, are
_automatically_ compressed to a similar or better amount arithmetically.

Therefore it is easier to tune for particular data storage, whether that's by
adjusting knobs at compile time, or it can be made self-tuning.  Self-tuning
works even for fast search tree structures and large databases where the data
is not shaped like a stream, and we don't have to start with good probability
estimates; they can be discovered.

In effect, we exchange some general purpose arithmetic overhead for a toolkit
that greatly simplifies the encoding of complex, compressed or delta-compressed
data structures.

## Unsuitability of ANS with tree-structured compression

[Asymmetric Numeral
Systems](https://en.wikipedia.org/wiki/Asymmetric_numeral_systems) as used in
Zstd compression was examined.

In practice it is faster than arithmetic coding and uses simpler state, but it
requires the compressor to run in the opposite direction to the decompressor.
This is called LIFO, and the regular compression direction is called FIFO.

For the purposes of _tree-structure_ compression, we can't really use LIFO.
It's infeasible to work backwards from leaf nodes to their parents and the
root, because the compression states don't merge.  In FIFO, compression states
in a tree branch rather than merge.

## Alternative techniques

Ad-hoc byte-codes were explored, as described above.  The combinatorial
explosion and difficulty (and time) tuning for large, specialised data sets
motivated further work on arithmetic coding.

I looked into a number of alternative techniques.  ANS (Asymmetric Numeral
Systems), Galois field (carryless) multiplication, carryless (polynomial)
addition, P-adic numbers, other non-linear but monotonic maps, non-monotonic
maps, and LFSR noise to make the non-monotonic probabilities fairer.

In the process I rediscovered something like ANS (though I knew about ANS
already), and may have found hints of an ANS-like method that can be compressed
in the forward direction, as well, but the compressor looks too complex and
slow.  Tree-structure compression for a _fast-updateable_ database requires
compression to be reasonably fast as well, not just decompression.

In the end, for tree-structure applications, I arrived at a new appreciation
for conventional number arithmetic, monotonicity, the _prefix property_ which
says that earlier symbols are directly determined by earlier bits, even though
they affect later bits, and the _bounded compression lookahead_ property which
says that the compressor can emit symbols without having to go through a lot of
state for later symbols before back-patching earlier output.

Bounded compression lookahead is deviated from in this compressor, although it
is a small and controlled deviation.  This is due to the _carry-counting
algorithm_, where compressed digit output is deferred if later arithmetic might
overflow and carry into the deferred digit.  The amount of state needed to do
this in the compressor is small.  The decompressor does not care about it.

There is a different method called _bit stuffing_ (or _digit stuffing_
sometimes) which allows the output to always progress, and no symbol is
withheld.  Bit stuffing uses a little less compression state, is a little
faster to compress, and allows guarantees continuous progress which is
necessary for some streaming applications.  But there is a small compression
ratio penalty and decompression is slightly slower.  So I chose carry-counting,
but bit stuffing may be added and tested.

## License

These files are licensed under either of ["MIT
license"](http://opensource.org/licenses/MIT) or ["Apache License, Version
2.0"](http://www.apache.org/licenses/LICENSE-2.0), at your option.

These files are provided without any warranty.  Use at your own risk.  It is
intended that excerpts be used and changed in other programs, subject to the
terms of the one or both of the above licenses.

You can use them under either one license or both as you wish when
incorporating the code into your own work.  You _may not_ simply copy and paste
the code without attribution.

## History

This is R&D to provide and bring together working primitives for (and from)
other projects.

As with many of my projects the initial Git commit is dated after it's been in
use for some time.  These techniques have been used for years here and there,
but I didn't bring them together and benchmark, optimise sizes and word mode
etc. until 2021.

Compression is not really the point.  It's necessary, but the point is fast,
updatable search trees that are compressed.  This one is related to the `dbio`
fast storage I/O project and tree-structured storage (B-tree designs etc),
going back to Jan 2008.  It was actually part of a distributed database, with
storage being secondary but important.  If we look further my interest in this
area was started with storage integrity on hard drives in Linux:
https://lkml.org/lkml/2008/5/20/320.
