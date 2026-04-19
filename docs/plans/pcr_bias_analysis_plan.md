# PCR Bias Analysis: Measuring Impact at Two Levels

## Problem Definition

PCR amplification of long cassettes (e.g., PEtracer 5mer ~1kb) produces fragments of varying lengths due to amplification bias towards shorter products. This creates:

1. **Over-representation of short fragments** - Most reads cover only the beginning of the cassette
2. **Under-representation of full-length reads** - Few reads span the entire cassette

Observed in allele table (`scc_6_6_alleles_pl.parquet`):
- r1-r4: 95-100% coverage
- r5-r7: ~55% coverage
- r8-r15: 9-19% coverage
- 42% of molecules stop at r4, only 9% cover full cassette

## Level 1 Findings: Reads → Molecules (Segmentation Assembly)

**Conclusion: Segmentation assembly is faithful but input-limited.**

### Read length distribution (scc_6_6: 10.5M reads → 1M contigs)

Bimodal distribution confirmed:
- Dominant peak at 400-450bp (short PCR products)
- Long tail decaying to 2000bp
- p50 = 446bp, p90 = 772bp, p95 = 934bp (cassette ~1kb)

### Per-UMI read composition (4.15M UMIs, median reads/UMI = 1)

| Category | UMIs | % |
|---|---|---|
| All short reads (<500bp) | 1,818,384 | 43.8% |
| Medium only (500-900bp), no long | 1,741,147 | 41.9% |
| At least 1 long read (>=900bp) | 591,657 | 14.3% |

### PCR artifact reads

Checkhealth analysis reveals a large class of reads containing cassette start+end anchors
but zero META sequences (PCR artifacts). These reads:
- Produce **zero segments** in segmentation (edge segments require at least one META)
- Are silently dropped from assembly - they cannot corrupt the contig
- Still count toward the UMI's read count, inflating apparent support
- The only explicit filtering is `flag_aberrant_molecules` which catches tandem duplications
  (any META appearing >1 time), not no-META artifacts

### Assembly fidelity: 100% segment preservation

Every segment type present in raw segments makes it through assembly with zero loss:
- 1,003,508 UMIs compared
- Same segment count in raw vs assembled: **100.0%**
- Lost segments: **0 (0.0%)**
- Mean segment types per UMI: 2.7 (limited by read length, not assembly)

### Key insight: short reads cannot hurt, but cannot help

In segmentation-based assembly:
- Each read is split at META boundaries into independent segments
- Assembly happens per `(umi, start_meta, end_meta)` group
- Short reads contribute segments only for early METAs - they are **absent** from later
  segment groups, not competing with longer reads
- A short read covering META01→META04 adds support for segments (01,02), (02,03), (03,04)
  but is invisible to segments (04,05), (05,06), etc.
- Therefore short fragments cannot "pull down" the contig or corrupt later positions

**The bottleneck is input coverage, not assembly quality.** Contig completeness is purely
determined by whether the UMI received at least one read long enough to cover later METAs.

## Level 2 Findings: Molecules → Alleles (UMI Collapse)

### Allele table structure (scc_6_6: 16,735 molecules, 34 intBCs, 816 cells)

Allele values are complex strings encoding the edit at each cutsite:
- Insertions: `TAGTCATTAC[78:5I:ACTCC]ACTCCGA...`
- Wild-type: contains `[None]`
- Deletions: `GCTCACCTAT[437:581D]` (581bp deletion at position 437)

### Consecutive position equality reveals dominant spanning deletions

| Positions | Identical allele (%) | Interpretation |
|---|---|---|
| r1 == r2 | 2.6% | Independent edits (healthy) |
| r3 == r4 | 1.6% | Independent edits (healthy) |
| r4 == r5 | 35.8% | Spanning deletions begin |
| r7 == r8 | 85.7% | Most molecules share same deletion |
| r10 == r11 | 91.7% | Nearly all molecules identical |
| r14 == r15 | 93.1% | Nearly all molecules identical |

82.6% of molecules have r7==r8==...==r15, all carrying the same allele value.

### Two dominant deletion breakpoints across all intBCs

The same two deletions dominate every intBC in scc_6_6:
- `GCTCACCTAT[437:581D]` - 581bp deletion at position 437
- `GCTCACCTAT[253:765D]` - 765bp deletion at position 253

These are NOT truncation artifacts: 86.5% of molecules carrying these deletions
contain the cassette end anchor in their aligned sequence, confirming the contigs
are full-length with a real internal deletion. Position 437 is a structural hotspot
in the cassette reference — the deletion occurs independently across different intBCs
and samples as a genuine Cas9-mediated resection event.

### Cross-library comparison: deletion enrichment in single-cell

Comparing the same intBC (`TCGAAGCCCTATTGTGGTA`) across libraries:

| Library | Type | [437:581D] | [253:765D] | Wild-type at r7 |
|---|---|---|---|---|
| **F7_1_FACS** | single-molecule | **2.1%** | 0.9% | **72.4%** |
| **Fpp4** | single-molecule | **0.3%** | 0.2% | **70.6%** |
| **scc_6_6** | pooled single-cell | **38.0%** | 32.7% | low |

**Disproven hypothesis: cross-sample PCR contamination.** Early alleles (r1-r4) of
deletion-carrying molecules at F7-specific intBCs match F7 biology (74.8% wild-type
at r1), NOT 4T1_11 biology (89% edited at r1). The intBC is part of the genomic
segment, so cross-intBC contamination cannot occur via standard mechanisms.

### CRITICAL FINDING: Ambient DNA contamination in single-cell library

**The deletion molecules are ambient DNA — free-floating molecules misassigned to
wrong cell barcodes during droplet encapsulation.**

#### Evidence 1: Uniform contamination across cells

96.2% of cells are "mixed" at their intBC — carrying both deleted and non-deleted
molecules. This is biologically impossible since each cell has ONE genotype per
integration site. The deletion fraction is consistent within each intBC (std ~0.14)
but varies across intBCs (6-86%), exactly as expected from a shared ambient pool.

| Metric | Value |
|---|---|
| Cells that are "mixed" (both del + non-del) | 96.2% |
| Cells with 100% deletion | 1.7% (20 cells, 1 intBC) |
| Cells with 0% deletion | 3.1% (mostly 1 intBC) |
| Mean deletion fraction across intBCs | 0.43 +/- 0.18 |

#### Evidence 2: Early allele fingerprint mismatch

Using r1+r2+r3 as a cell-specific fingerprint:
- **Only 4.7%** of deletion molecules match their host cell's fingerprint
- **20.2%** of non-deleted molecules match their own cell's mode (baseline)
- Deleted molecules are overwhelmingly from DIFFERENT cells than assigned

#### Evidence 3: UMI counts are similar

Read counts per molecule are comparable between deleted (median=3, mean=4.2) and
non-deleted (median=4, mean=4.8). UMI deduplication is working — the inflation is
NOT from PCR amplification of the same molecule. Each ambient molecule has a
genuine unique UMI because it IS a real molecule from a different cell.

### Mechanism: pre-capture physical enrichment

1. Free-floating DNA exists in the cell suspension before droplet encapsulation
2. Shorter cassettes (deletion removes 581bp) are physically more stable and more
   easily captured as ambient DNA
3. This enriches the ambient pool from ~2% deletion (bulk rate) to ~38%
4. During droplet encapsulation, ambient molecules get random cell barcodes + unique UMIs
5. UMI deduplication cannot detect this — each is a genuinely distinct molecule

### Impact on collapsed alleles

After mode-voting collapse in scc_6_6:
- r7-r15: **72.4% carry ambient deletion alleles**, 27.6% null
- Zero real allele diversity survives at later positions
- The collapse correctly implements mode voting — but the input is contaminated
  by ambient DNA that overwhelms the genuine signal

### Within 4T1_11: deletions are intBC-specific (real biology)

In the 4T1_11 single-molecule library (no single-cell capture, no ambient DNA issue),
the deletion pattern varies by intBC as expected from independent biology:

| intBC | [437:581D] | [253:765D] | Wild-type |
|---|---|---|---|
| TCCTAAGACATTTTGTGTA | 72.5% | 0.8% | 6.6% |
| TTATACCCCAGGTTTAGTA | 3.9% | 79.3% | 2.3% |
| TCTGAACACAATTCGTGTA | 8.2% | 4.0% | 35.2% |
| TCAGACGGCGAATAGAGTA | 3.7% | 1.7% | 22.4% |

## Implications and Next Steps

### The PCR bias problem has three distinct effects:

1. **Read-level bias (Level 1):** Short PCR fragments limit contig completeness.
   Segmentation assembly handles this correctly — it preserves all available segments
   and cannot be corrupted by short reads. The limitation is purely input coverage.

2. **Ambient DNA enrichment (Level 2):** In single-cell libraries, deletion-carrying
   molecules are physically enriched in the ambient DNA pool because shorter cassettes
   are more stable as free-floating DNA. This creates uniform contamination across ALL
   cells at each intBC. UMI deduplication cannot detect this because each ambient
   molecule has a genuine unique UMI.

3. **Signal destruction (Level 2 consequence):** Ambient deletion molecules dominate
   mode voting at later positions (r7-r15), collapsing all allele diversity to the
   deletion allele. Real lineage signal is destroyed at these positions.

### Recommendations:

1. **Ambient DNA estimation** per intBC: the fraction of "mixed" cells and the
   uniformity of the deletion fraction across cells directly estimates the ambient
   contamination rate
2. **Per-cell fingerprint filtering:** remove molecules whose r1-r3 fingerprint
   doesn't match the cell's consensus (only 4.7% of deletion molecules match,
   vs 20.2% baseline)
3. **SoupX or similar ambient RNA/DNA correction:** apply established single-cell
   ambient correction methods, using the deletion alleles as known ambient markers
4. **Position-aware collapse:** at later positions (r7+), discount or exclude alleles
   that match the known ambient deletion signature before mode voting
5. **Compare with single-molecule libraries:** use F7_1_FACS / Fpp4 as ground truth
   for the expected deletion rate (~2%) to calibrate contamination estimates

## Files Analyzed

All paths relative to `/home/projects/nyosef/pedro/projects/lt/workdir/20260202_pe/`

| File | Sample | Format | Content |
|------|--------|--------|---------|
| `scc_6_6/parsed_reads.arrow` | pooled SC | Arrow IPC | 10.5M raw reads |
| `scc_6_6/intermediate/segments_debug.parquet` | pooled SC | Parquet | Raw segments per read |
| `scc_6_6/intermediate/assembled_debug.parquet` | pooled SC | Parquet | Assembled segments per UMI |
| `scc_6_6/contigs_segmented_valid.arrow` | pooled SC | Arrow IPC | 1M contigs |
| `scc_6_6/alleles_pl.parquet` | pooled SC | Parquet | 16,735 molecules |
| `scc_6_6/alleles_pl_collapsed.parquet` | pooled SC | Parquet | 1,328 collapsed alleles |
| `4T1_11/alleles_pl.parquet` | 4T1_11 SM | Parquet | 747,425 molecules |
| `F7_1_FACS/alleles_pl.parquet` | F7_1_FACS SM | Parquet | 470,853 molecules |
| `Fpp4/alleles_pl.parquet` | Fpp4 SM | Parquet | 342,033 molecules |

## Verification

Run analysis on scc_6_6 sample data:
- 10,480,596 reads → 4,151,188 UMIs → 1,003,508 contigs → 16,735 alleles
- Cross-validated against 4T1_11, F7_1_FACS, Fpp4 single-molecule libraries
