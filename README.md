# GC-VOC molecular programs

This package studies public human multiomics associations and molecular programs in gastric cancer. The literature-based C4–C10 aldehyde question is an external anchor. The study does not directly identify or measure VOCs, establish their tissue source, demonstrate causal VOC generation, or provide functional validation.

The working tree contains complete permitted effect, hypothesis, family, uncertainty and diagnostic tables; qualified population and pairing records; 44 figure-source tables; maintained implementations for all seven reported analysis modules; preparation functions; and portable environment specifications. Original registered work (146 hypotheses, 29 families) remains separate from result-informed supplementary work (42 hypotheses, 12 families). Negative, conflicting, uncertain, inactive and non-estimable results are retained.

From the package root:

```sh
./verify
# Set GCVOC_RSCRIPT to an existing qualified Rscript executable.
workflow/environment/run-r workflow/pipelines/check_methods.R
# OUTPUT_DIRECTORY must be absent or empty, outside this release.
workflow/environment/run-r workflow/pipelines/reproduce_core_bulk.R "$OUTPUT_DIRECTORY"
```

The last command independently calculates four prespecified core bulk RNA estimates from the bundled patient scores. It checks 112 GSE29272 pairs and 49 GSE27342 core pairs, the saved tumor-minus-reference differences, and separate two-hypothesis BH families. It does not read the accepted estimates during calculation. This is a tested **patient-score-to-estimate replay**; the raw-accession-to-score route has not been executed for this release.

[REPRODUCIBILITY.md](REPRODUCIBILITY.md) documents every module, explicit input manifests, available preparation code, required external inputs, pinned profiles, synthetic checks and two optional plotting modes. `./verify` needs only Python 3's standard library, works offline and does not depend on Git history. No entry point downloads data or searches another project for missing inputs.

Start with the [original hypothesis ledger](results/integration/original/hypothesis_coverage.tsv), [coverage summary](results/integration/coverage_summary.tsv), [family reconciliation](results/integration/family_reconciliation.tsv), [source dependence](results/integration/source_dependence.tsv), and [publication source index](publication/INDEX.md). Historical identifiers containing `64` or `80` do not specify final biological sample sizes.

Finished manuscript figures and formatted tables are not in the working-tree payload. Reader-generated plots must be written outside the release. LinkedOmics matrices/linked participant values, audited GSE206785 individual/cell derivatives and donor mappings, publisher attachments, controlled-access data, private bindings and opaque analysis containers remain excluded. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for object-specific conditions; no anonymization or blanket third-party license is claimed.

[PUBLIC_MANIFEST.tsv](PUBLIC_MANIFEST.tsv) covers every physical payload file. Original-code licensing remains an owner decision. This is not a claim of finalized licensing, a fully open-source release, publication approval or end-to-end reproducibility. Working-tree exclusions do not remove files from earlier Git commits.
