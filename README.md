# gramene-bedify
Convert gramene gene JSON documents into tab delimited bed format. Two parameters, bedFeature and bedCombiner define the behavior.

bedFeature defaults to 'gene', but other options include: 'exon', 'intron', 'CDS', 'utr5', or 'utr3'
bedCombiner defaults to 'canonical', but 'all' is available too

[![Build Status](https://travis-ci.org/warelab/gramene-bedify.svg?branch=master)](https://travis-ci.org/warelab/gramene-bedify)