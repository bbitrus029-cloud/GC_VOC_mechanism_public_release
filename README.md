# GC-VOC molecular programs

This package studies molecular associations in public human gastric-cancer multiomics. It does not measure VOCs, establish their tissue source, demonstrate causal VOC generation, or provide functional validation.

Three scopes remain separate: the original study has **146 hypothesis slots / 29 formal family instances** (136 formal members plus 10 no-formal-family slots); the first post-effect supplement has **42 / 12**; the second has **10 / 2**. The derived inventory has **198 slots / 43 family instances**, with no pooled or study-wide FDR claim. Registrations are [STQJC](https://osf.io/stqjc/), [GE7BK](https://osf.io/ge7bk/) and [C6KE7](https://osf.io/c6ke7/) (DOI 10.17605/OSF.IO/C6KE7). The public C6KE7 adoption cover resolves the frozen plan's historical DRAFT_FOR_REVIEW label without changing its scientific body.

The payload contains permitted full-precision results, uncertainty, diagnostics, population/pairing records, portable methods, source-specific preparation functions and current figure/table source navigation. Negative, conflicting, unavailable and non-estimable results remain explicit. PE2 B is a technical sensitivity outside the ten new slots; A/C/K add no inferential endpoints. F/G remain deferred before registration; H/I/J and K RNA are excluded.

Run from the release root:

```sh
./verify
workflow/environment/run-r workflow/pipelines/check_methods.R
workflow/environment/run-r workflow/pipelines/reproduce_core_bulk.R "$CORE_OUTPUT"
# Header-only external input manifest suffices for each bundled PE2 replay.
workflow/environment/run-r workflow/pipelines/second_post_effect.R B "$INPUT_MANIFEST" "$B_OUTPUT"
workflow/environment/run-r workflow/pipelines/second_post_effect.R D "$INPUT_MANIFEST" "$D_OUTPUT"
workflow/environment/run-r workflow/pipelines/second_post_effect.R E "$INPUT_MANIFEST" "$E_OUTPUT"
```

Use the existing qualified R runtime/libraries and distinct absent or empty output directories outside this tree. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) distinguishes actually tested value-to-estimate routes, synthetic tests, historical evidence and unexecuted external-input routes. `./verify` is offline, read-only and Python-standard-library-only. Entry points do not download or search a scientific source project for missing inputs.

Current Figures 2–8 and table section IDs are mapped in [publication/INDEX.md](publication/INDEX.md). Complete Figures 3, 6 and 7 can be rendered from bundled sources. Figures 2, 4, 5 and 8 support explicitly labelled eligible-panel exports; complete rendering accepts validated scientific intermediates from the documented existing-analysis-output adapters. Those held data remain external, and upstream reconstruction limits remain explicit. The original two plotting modes remain available. Finished figures, formatted tables, previews, manuscripts and editorial masters are excluded from the payload; user-generated outputs belong outside it.

Start with the [original ledger](results/integration/original/hypothesis_coverage.tsv), [historical original/first-supplement coverage](results/integration/coverage_summary.tsv), [PE2 ledger](results/supplementary/second_post_effect/endpoint_ledger.tsv), [current all-scope disposition](publication/RESULT_DISPOSITION.tsv) and [source dependence](results/integration/source_dependence.tsv). Historical `64`/`80` identifiers do not specify final biological n.

LinkedOmics participant-linked data, held GSE206785 individual/cell derivatives and donor mappings, publisher attachments, controlled data, private bindings and opaque containers remain excluded. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Original-code licensing remains an owner decision; this is not unrestricted open source, a licence certification, independent acceptance or an end-to-end reproducibility claim. [PUBLIC_MANIFEST.tsv](PUBLIC_MANIFEST.tsv) covers every physical payload file. Working-tree exclusions do not remove earlier Git history.
