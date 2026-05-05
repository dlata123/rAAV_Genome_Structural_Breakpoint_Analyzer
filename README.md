# AAV Genome Structural Analysis (Snapback + RNAfold)

This repository contains a workflow to analyze AAV genome structures from long-read sequencing data (PacBio). The focus is on identifying snapback breakpoint regions and checking whether these regions are associated with local secondary structure using RNAfold (MFE).

This work is based on approaches described in:

- Structural Analysis of Recombinant AAV Vector Genomes at Single-Molecule Resolution  
---

## Background

AAV genomes produced during manufacturing are often heterogeneous. Along with full-length genomes, truncated forms, snapback structures, and other rearrangements are commonly observed.

With long-read sequencing, these structures can be observed directly at the single-molecule level. To interpret them, tiling and subparser steps are used.

The tiling algorithm aligns each read to reference components such as ITR and payload and breaks it into ordered segments.  
The subparser then simplifies these patterns and extracts breakpoint positions and counts.

The subparser output for snapback and truncation(`*.tile.zmw.counts`) is used as input for this workflow.

---

## What this workflow does

- reads tile count data  
- extracts breakpoint positions in the payload  
- separates plus (coding) and minus (non-coding) strand events  
- generates sequence windows around breakpoints  
- runs RNAfold  
- extracts MFE values  
- combines MFE with breakpoint frequency 
- perform statistical analysis 
- generates summary tables and plots  

---

## Files in this repository

- example.snapback.tile.zmw.counts  
  Example subparser output  

- reference_split.fa  
  Reference FASTA with ITR and payload (synthetic example)

- snapback_analysis.ipynb  
  Main notebook for snapback analysis  

- truncation_analysis.ipynb  
  Workflow for truncation analysis  

- sample_output/  
  Example output files  

---

## Requirements

- Python  
- pandas  
- numpy  
- biopython  
- scipy  
- RNAfold (ViennaRNA)

Optional:
- rpy2  
- R with ggplot2  

---

## How to run

Update paths at the top of the notebook:

```python
from pathlib import Path

PROJECT_DIR = Path(".")

INPUT_TILE_FILE = PROJECT_DIR / "example.snapback.tile.zmw.counts"
REFERENCE_FASTA = PROJECT_DIR / "reference_split.fa"
OUTPUT_FOLDER = PROJECT_DIR / "sample_output"

Then run the notebook step by step.

---
Output

The workflow generates:

breakpoint tables
RNAfold structure output
MFE values
combined datasets
plots showing breakpoint distribution across the payload

MFE values are used in plots, where more negative values indicate more stable secondary structures.

Main idea

This analysis looks at whether snapback breakpoints occur more often in regions with specific secondary structure properties.

Notes
Reference files are synthetic and used for demonstration
RNAfold must be installed and available in PATH

Author
Dimpal Lata
Scientist