# Current publication source navigation

The current dual-journal mapping is authoritative. Legacy filenames retain historical source identities; they are not current panel labels. [SOURCE_INDEX.tsv](SOURCE_INDEX.tsv) resolves current roles and historical support to canonical public tables. [figure_inputs.tsv](figure_inputs.tsv) records every source dependency, actual public schema and reconstructed scientific identity; [PANEL_MAPPING.tsv](PANEL_MAPPING.tsv) records all 32 current panels and their exact requirements. There is no fixed figure-source file-count expectation.

| Figure | Scientific content | Bundled rendering | Additional input for a complete figure |
|---|---|---|---|
| 2 | Tissue-level program changes | A, C, D | B: LinkedOmics-derived paired protein G6/H4 scores (28/30 patients) |
| 3 | Heterogeneity within G | Complete A–E | None |
| 4 | H expression across molecular layers | A, C, D | B: CPTAC patient protein values; Mun values are bundled but do not make the two-cohort panel complete |
| 5 | Patient-level RNA–protein associations | C–F (PE2 E; four members, same 80 patients) | A/B: held CPTAC paired G/H bridge values (19/21 patients) |
| 6 | Gastric-expression context | Complete A–F; PE2 D in D/E | None |
| 7 | Cellular expression context | Complete A–D; PE2 B in C | None |
| 8 | Captured-cell composition | B/C | A: held GSE206785 donor composition; donor navigation is also external |

Current presentation includes percentage histograms in Figure 2A, native-coefficient bars with stored intervals in Figure 4A, patient boxplots in Figure 4B, and confidence-interval crossbars in Figure 6D. These operations do not change scientific estimates. Existing rank H, historical raw-count standardized H and PE2 normalized H4 remain distinct. K is IHC-only supporting context with no invented main-figure assignment.

[TABLE_DIRECTORY.tsv](TABLE_DIRECTORY.tsv) preserves the current section IDs and **10 November 2026 Directory figure-reference correction** from `layout.json`: 105 workbook sections and three concise reader sections. Project TableS3 is submission Supplementary Table 1; project Tables S1/S2 are submission Supplementary Tables 2/3 and remain outside this source-rendering package; S4–S7 keep their numbers. Section labels do not determine figure references; `related_figures` is carried by stable section ID.

[table_layout.json](table_layout.json) contains eligible current display rows only. Each retained row was screened through every source cell binding; [table_marks.tsv](table_marks.tsv) maps it to public path/row/field tokens. `source_rows`, `eligible_rows` and `held_rows` distinguish full source coverage from the permitted display projection. Held rows are not zero or unavailable scientific outcomes. Full scientific endpoint coverage is separately retained in [RESULT_DISPOSITION.tsv](RESULT_DISPOSITION.tsv); table-layout eligibility is not endpoint selection. No private accepted/display container, omission-alias map, publisher attachment or operational binding container is bundled.

Canonical full-precision scientific tables are the numerical authority. Reused C results add no independent evidence. The source index retains unique diagnostic and supporting rows even when no current main panel draws them. Historical display support is stored in `source/legacy_support/` with its original identity explicitly identified.

The seven held inputs have public scientific projections and adapters from existing analysis outputs; their historical private-view hashes are labelled provenance only. The schema omits private bindings and obsolete coordinates, while independently frozen typed-content checks protect complete populations and values. Full Figures 2/4/5/8 remain conditional on lawfully supplied qualified scientific inputs; only their output-to-panel route has been tested with private fixtures.

Rendering commands and exact prerequisites are in [REPRODUCIBILITY.md](../REPRODUCIBILITY.md). No finished graphics, DOCX/XLSX tables or previews are payload files.

The eligible workbook retains current section notes, including units, denominators, fixed families and dependence qualifications. Held rows include source-label detail from excluded clinical inputs and cells bound to internal qualification records. Aggregate denominators and complete released numerical results remain in the canonical public tables. Native Excel/Calc application validation has not been performed; persisted OOXML cells, source links and font geometry were inspected.
