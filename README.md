# OBP Gene Family: Genomic Localization & Visualization

Course project for *Biological Omics & Big Data* — genome-wide identification and chromosomal mapping of odorant-binding protein (OBP) genes in the silkworm (*Bombyx mori*).

## What This Does

1. **Data ingestion**: Reads OBP gene family data (gene names, chromosomal locations) and chromosome length information
2. **Intelligent clustering**: Automatically groups genes by chromosome and identifies gene clusters (≥4 genes clustered on the same chromosome)
3. **SCI-style visualization**: Generates a publication-quality chromosome ideogram with gene positions, cluster annotations, and individual gene labels

## Key R Packages

`ggplot2` · `dplyr` · `ggrepel`

## Output

![OBP Chromosome Map](figures/OBP_Genomic_Locations_SCI_Style.png)

- Gene positions marked as red lines across chromosome ideograms
- Gene clusters annotated with brackets and cluster IDs
- Isolated genes labeled with gene names (ggrepel, italicized)
- Exported as both PDF (vector, journal-ready) and PNG (raster preview)

## Files

```
.
├── OBP_location.R              # Main analysis script
├── data/
│   ├── obp_family.csv          # Gene family annotation
│   ├── obp_locations.txt       # Gene chromosomal positions
│   ├── obp_positions.txt       # Detailed positions
│   ├── obp_detailed_positions.txt
│   └── BMSK_chr_lengths.txt    # Chromosome length reference
└── figures/
    ├── OBP_Genomic_Locations_SCI_Style.pdf
    ├── OBP_Genomic_Locations_SCI_Style.png
    └── OBP_Chromosome_Map.png
```

## Usage

```r
source("OBP_location.R")
```

The script will prompt you to select input data files interactively, then generate figures automatically.
