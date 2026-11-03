# Scientific source navigation

This working-tree package contains 44 retained figure-source tables and complete permitted numerical results, with canonical aliases in [SOURCE_INDEX.tsv](SOURCE_INDEX.tsv). It contains no finished manuscript figures, formatted tables, previews or cached renders. Previously included Figure 2/4/8 PDF/SVG files, Table S3 DOCX and the S4–S7 workbook have been withdrawn from the working-tree payload. Earlier Git commits are not changed by that withdrawal.

| Item | Retained scientific content and boundary |
|---|---|
| Figure 2 | Complete eligible effect/member/site tables and observed patient differences in [Figure2distributions.tsv](Figure2/source/Figure2distributions.tsv), primary paired scores and supplementary rank scores. The derived density grid is omitted. Its descriptive calculation is `stats::density(x, bw=stats::bw.nrd0(x), kernel="gaussian", n=512, from=min(x), to=max(x), cut=0)` per exact analysis ID; only outline height was scaled for display. |
| Figure 3 | Separately eligible aggregate protein/bridge source tables; linked individual inputs and complete exports remain withheld. |
| Figure 4 | Eligible score/member/subprogram/population diagnostic tables; repeated analyses are dependent evidence. |
| Figure 5 | Eligible cross-source reconciliation tables; unresolved individual protein components remain withheld. |
| Figure 6 | Eligible aggregate state/member findings; GSE206785 individual/cell derivatives remain withheld. |
| Figure 7 | Eligible source tables. [Figure7C.tsv](Figure7/source/Figure7C.tsv) contains 20 paired donors for each of `enteroendocrine` and `other_unassigned`, including original counts, denominators, proportions and logit differences. The former C pair/index displays use the enteroendocrine subset; D used both categories. Their separate index, pair-display and repeated individual-display copies are omitted. The source records retain scientific values, not display row positions. GSE206785 individual displays remain excluded. |
| Figure 8 | [Figure8A.tsv](Figure8/source/Figure8A.tsv) contains 80 GSE122401 patients; [Figure8B.tsv](Figure8/source/Figure8B.tsv) contains 33 TCGA patients, with exact paired identities and H/D7 changes. The separate concatenated patient display-order file is omitted. Other fixed-score, incomplete-cohort and influence sources remain. These are associations, not mediation. |
| Table S3 | Essential clinical/reference, eligibility and dependence facts remain in their canonical source/result records, scientific configurations and REPRODUCIBILITY.md. GSE150290 methods retain processing uncertainty, result-informed timing and lineage qualification. No formatted table is recreated. |
| Tables S4–S7 | Full-precision complete effects, family/uncertainty and diagnostic tables under `results/`; no formatted workbook. |
| Figure 1; Main Table 1; S1–S2 | Outside this package; no items imported or invented. |

[The data dictionary](tables/source/data_dictionary.tsv) preserves definitions. Its binding examples describe internal production records; use [SOURCE_INDEX.tsv](SOURCE_INDEX.tsv) for bundled files. The index explicitly distinguishes canonical sources from replacement references to equivalent values; omitted derived grids and display order are not claimed to be byte-identical aliases.

[The plotting interface](../workflow/pipelines/render.R) can generate two tissue-context plots and one protein-member plot from saved permitted sources into an external directory. It does not recreate the complete manuscript packages or bypass holds. No manuscript, final legends or manuscript placement is implied.
