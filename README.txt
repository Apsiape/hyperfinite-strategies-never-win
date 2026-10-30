# A finite game that normal hyperfinite strategies approach but never win

Seth Douglas and Nidhal Mghirbi

Version 1.0.0, October 2026. DOI: 10.5281/zenodo.23816141. Concept DOI (all versions): 10.5281/zenodo.23816140. A preprint; it has not been peer reviewed.

## Abstract

We exhibit non-closedness at perfect play for normal injective hosts: an explicit finite two-player game with 0/1 scoring has a perfect commuting-operator correlation approached by finite-dimensional quantum correlations, while no perfect strategy has a normal marginal on an injective von Neumann algebra of either player. The game has 561 Alice questions, 157 Bob questions, 3,696 sampled pairs and at most 20 answers per question. Perfect play is realized on a matrix amplification of the Bellaterra critical-KMS factor of type III with parameter 1/2. After restricting to the support of the chosen normal marginal, every perfect strategy yields a modular Cuntz pair and six self-similar relations that force the canonical trace on the free product of three groups of order two. The resulting noninjective group von Neumann algebra is an expected subalgebra of that support corner, with a faithful normal state-preserving expectation. A two-sided POVM lift through a matrix ultraproduct proves approximation by finite systems. Consequently correlations admitting normal injective hosts are not closed. Under the explicit hypotheses of Halvorson's Proposition 2.31, locally normal quantum-field-theory strategies cannot win exactly. We also formulate a KMS-face criterion for complete modular presentations. Exact arithmetic certificates check every supported losing event; the forcing proof applies to arbitrary normal strategies.

## Source and reproduction

The manuscript is paper/main.tex and the built paper is paper/main.pdf. The standalone source is paper/arxiv/. Both source trees contain the complete ancillary game and certificate programs.

Run `python -u paper/anc/run_checks.py --output-dir CHECK_DIRECTORY` with assertions enabled. Programs use the Python standard library and run sequentially; recorded peak memory stays below one billion bytes.

Build with `python -u paper/build.py --output-dir BUILD_DIRECTORY`. Use an empty build directory outside the source tree. The build uses pdflatex and BibTeX, rejects warnings and box defects, and records a PDF hash.

The full exact certificate and independent sample have different coverage, described in the paper and paper/anc/README.txt. The predicate has 561 and 157 questions, at most 20 answers, and 3,696 uniformly sampled pairs.

## Rights

Text and documentation: CC BY 4.0. Scripts: MIT or Apache 2.0, at your option. See RIGHTS.md and the license files.
