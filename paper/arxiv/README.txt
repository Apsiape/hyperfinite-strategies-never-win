A finite game that normal hyperfinite strategies approach but never win
Seth Douglas and Nidhal Mghirbi
Version 1.0.0, October 2026. A preprint; it has not been peer reviewed.

Abstract

We exhibit non-closedness at perfect play for normal injective hosts: an explicit finite two-player game with 0/1 scoring has a perfect commuting-operator correlation approached by finite-dimensional quantum correlations, while no perfect strategy has a normal marginal on an injective von Neumann algebra of either player. The game has 561 Alice questions, 157 Bob questions, 3,696 sampled pairs and at most 20 answers per question. Perfect play is realized on a matrix amplification of the Bellaterra critical-KMS factor of type III with parameter 1/2. After restricting to the support of the chosen normal marginal, every perfect strategy yields a modular Cuntz pair and six self-similar relations that force the canonical trace on the free product of three groups of order two. The resulting noninjective group von Neumann algebra is an expected subalgebra of that support corner, with a faithful normal state-preserving expectation. A two-sided POVM lift through a matrix ultraproduct proves approximation by finite systems. Consequently correlations admitting normal injective hosts are not closed. Under the explicit hypotheses of Halvorson's Proposition 2.31, locally normal quantum-field-theory strategies cannot win exactly. We also formulate a KMS-face criterion for complete modular presentations. Exact arithmetic certificates check every supported losing event; the forcing proof applies to arbitrary normal strategies.

This folder is a standalone source distribution. main.tex is the main file.
The section files, references.bib, main.bbl and complete anc/ folder are included.
No files from a parent folder are required.

Build with:
  python -u build.py --output-dir BUILD_DIRECTORY
Use an empty directory outside this source tree. pdflatex and BibTeX are needed.
The build rejects warnings, unresolved references and all box defects.

Check source consistency and distributed hashes:
  python -u check_artifacts.py
Reproduce all exact checks, preserving the distributed recorded outputs:
  python -u anc/run_checks.py --output-dir CHECK_DIRECTORY
Choose a fresh directory outside anc/. Assertions must remain enabled.

Text and documentation: CC BY 4.0. Scripts: MIT or Apache 2.0, at your option.
See RIGHTS.md and the three license texts.
