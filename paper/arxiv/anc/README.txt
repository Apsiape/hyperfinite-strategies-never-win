Ancillary finite predicate and exact certificates

The game table lists the valid answer arrays and forbidden pairs of zero-based
answer indices. All 3,696 supported question pairs have weight 1/3696.
Unsupported pairs have weight zero. Invalid answers are rejected.

Run with assertions enabled:
  python -u run_checks.py --output-dir CHECK_DIRECTORY
Choose a fresh check directory outside this source tree. The runner first checks
the distributed hashes, copies the portable inputs, then executes five programs
sequentially with unbuffered output and one numerical thread. It retains logs,
stops on any nonzero exit, and requires every recorded status to be PASS.

compile_and_check.py: exhaustive exact ideal-operator check over Q(sqrt(2));
718 question POVMs, 3,696 pairs, 90,630 forbidden occurrences, 19,014 representative
losing products. Cyclic symmetry gives 86 POVM and 672 pair templates.

independent_predicate_check.py: independent reconstruction of all question and
answer sets and all predicates, plus a 512-pair sample and a mutation check.
The default sample seed is recorded in independent-predicate-results.json;
GAME_CHECK_SEED may select a different seed for a new diagnostic run.

independent_ideal_check.py: separately implemented exact adjacent-generator
rewriting and base effects, on that same 512-pair sample. The recorded sample
covers 12,537 losing events and 485 question POVMs. This is a sample, not an
independent exhaustive ideal certificate.

audit_engine.py: 91,840 direct finite-word action comparisons and structural
checks. audit_forcing_table.py: comparison with the independently enumerated
base-giii-game.json and with the label and conditional requirements.
resource_usage.py only measures memory; it contains no algebra routines.
Each measured peak in the recorded runs is below one billion bytes. The
diagnostics do not establish the infinite-strategy forcing proof.

The *.log and *results.json files record actual outputs. SHA256SUMS uses LF
and covers the programs, input tables, documentation and recorded evidence.
The runner writes new results in the separate directory, leaving this record
unchanged. Text/documentation: CC BY 4.0; scripts: MIT or Apache 2.0.
