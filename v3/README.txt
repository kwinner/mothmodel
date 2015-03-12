---
- Distribution Objects
--

Throughout this version, we make an effort to permit arbitrary distributions. Unfortunately, existing Matlab procedures for distribution objects make it difficult to convert to log space computations, as well as being largely non-optimal for use in optimization problems. To address this limitation, we introduce our own probability distribution objects which effectively wrap a set of function pointers for the rng, pdf, cdf, and log pdf/cdf into a single wrapper for readability. Distribution parameters remain separate, since these are often the targets of independent optimization procedures.

Formally, our distribution objects have the following members:
	rand
	pdf
	cdf
	logpdf
	logcdf

More details of the argument specifications for these methods to follow.

Parameters should always be in the form of a cell array (unfortunately, there is no way around this in Matlab)
