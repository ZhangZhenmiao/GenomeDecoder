
# GenomeDecoder

## Overview

**GenomeDecoder** is a tool for identifying and reconstructing segmental duplications (SDs) within complete genomes. This tool was developed to tackle challenges in modern genomics, including the analysis of SDs in highly complex genomic regions that were previously unresolvable, such as immunoglobulin loci known for rapid evolution.

## Features
- **Identification of SDs**: GenomeDecoder accurately identifies SDs in complete genomes, facilitating research on genomic variations.
- **Evolutionary Analysis of SDs**: Supports the study of SD evolution with a focus on highly duplicated regions.
- **Comparative Genomics**: Enables comparasion of genomic architecture across complete genomes.

## Installation Steps

Download GenomeDecoder:
```bash
git clone https://github.com/ZhangZhenmiao/GenomeDecoder.git
```

Dependencies (except for LJA assembler) can be installed using Conda by creating a new environment `genomedecoder`:

```bash
cd GenomeDecoder
conda env create -f requirements.yml
conda activate genomedecoder
```

Build GenomeDecoder and LJA assembler:
```bash
# build GenomeDecoder
cd src && make
chmod +x consensus_asm edlib_align parse_cigar.py synteny synteny_analysis.py
cd .. && chmod +x GenomeDecoder

# build LJA (development branch required)
cd src/LJA
cmake .
make jumboDBG
cd ../../
```

## Usage

### Sequence preparation:
GenomeDecoder is a tool developed for complete genomes. Each input genome should be in FASTA format, and contain only one contig (in the case of a multi-chromosomal genome, all chromosomes/contigs will be concatenated in a single string using 2000 "N"s as separators automatically by GenomeDecoder). We highly recommend running RepeatMasker before inputting the genome files to GenomeDecoder.

To run RepeatMasker, the species name is required (if not specified, RepeatMasker will treat the sequence as human by default):

```bash
RepeatMasker [-species species_name] input_file.fa
```
This will create `input_file.fa.masked`. See more details in [RepeatMasker](https://www.repeatmasker.org/webrepeatmaskerhelp.html) manual.

### GenomeDecoder Command Syntax:
```bash
GenomeDecoder [-h] -g GENOME [-i ITERATIONS] [-k K_VALUES] [-s SIMPLE] [-c COMPLEX] -o OUTPUT
```
### Parameters:
- **-h, --help**  
  Displays the help message and exits.

- **-g GENOME, --genome GENOME**  
  Specifies the path to the genome file to be analyzed. Each input genome should contain a single contig and be masked using RepeatMasker (in the case of a multi-chromosomal genome, all chromosomes/contigs will be concatenated in a single string using 2000 "N"s as separators automatically by GenomeDecoder). Use `-g` multiple times to compare multiple genomes (generate synteny blocks shared among those genomes).

- **-i ITERATIONS, --iterations ITERATIONS**  
  Specifies the number of iterations to use, applicable for disembroiling more than two input genomes. The default value is 5.

- **-k K_VALUES, --k_values K_VALUES**  
  Defines K-mer values for graph transformations. The default values are `31,51,101,201,401,801,2001`. For highly complex regions, another recommended setting is `-k 21,25,31,51,101,201,401,801,2001` to produce more detailed SD blocks.

- **-s SIMPLE, --simple SIMPLE**  
  Sets the similarity threshold for collapsing simple bubbles. The default value is `0`, which is suitable for most datasets.

- **-c COMPLEX, --complex COMPLEX**  
  Sets the similarity threshold for collapsing complex bubbles. The default value is `0.65`, which is suitable for most datasets.

- **-o OUTPUT, --output OUTPUT**  
  Specifies the output directory where GenomeDecoder will save the results.

## Outputs

In the output directory, GenomeDecoder generates files named `final_blocks_<i>.csv`, where each `i` corresponds to the synteny blocks of the i<sup>th</sup> input genome, with a minimum block length of 2 kb.

Each `final_blocks_<i>.csv` file contains 12 columns:

- **Synteny Block in Transformed Sequence**: Identifies each block-instance with a unique format: `block_ID(multiplicity_in_disembroiled_graph)length_in_kb`.
  
- **Start Vertex in Transformed Graph**: The starting vertex of the block-instance in the disembroiled graph.

- **End Vertex in Transformed Graph**: The ending vertex of the block-instance in the disembroiled graph.

- **Start Pos in Transformed Sequence**: The start coordinate of the block-instance in the disembroiled genome.

- **End Pos in Transformed Sequence**: The end coordinate of the block-instance in the disembroiled genome.

- **Include Left Vertex**: Indicates if the block-instance includes the start vertex (1 for yes, 0 for no).

- **Include Right Vertex**: Indicates if the block-instance includes the end vertex (1 for yes, 0 for no).

- **Is Short Version**: Specifies if the block-instance is a shortened version of the block (1 for yes, 0 for no).

- **Start Pos in Original Sequence**: The start coordinate of the block-instance in the original genome, obtained via Edlib. In case of a multi-chromosomal genome, note that the contigs/chromosomes have been concatnated by 2000 "N"s into a single sequence.

- **End Pos in Original Sequence**: The end coordinate of the block-instance in the original genome, obtained via Edlib. In case of a multi-chromosomal genome, note that the contigs/chromosomes have been concatnated by 2000 "N"s into a single sequence.

- **Similarity of Transformed Block and Original Block**: The percent identity between the block-instance in the disembroiled genome and its corresponding instance in the original genome, calculated by Edlib.

- **Length in Original Sequence**: The length, in base pairs (bp), of this block-instance in the original genome.

### Optional Analysis with `synteny_analysis.py`

Under the `src` directory of this repo, there is an optional script, `synteny_analysis.py`, for transforming CSV outputs into block representations using letters (A-Z). This script is suitable for cases with fewer than three input genomes. **Note:** If the number of blocks exceeds 26, the script will produce an error.

#### Usage:
For two genomes:
```bash
synteny_analysis.py final_blocks_1.csv final_blocks_2.csv
```

For a single genome:
```bash
synteny_analysis.py final_blocks_1.csv
```

The output will be printed to stdout.

### Example Usage

The following command generates synteny blocks for the RepeatMasker-processed files `test_example/HG38.IHG.Masked.Ns_transformed.fa` and `test_example/Orang.IGH.Masked.Ns_transformed.fa`, saving the outputs to `test_example/output/final_blocks_1.csv` and `test_example/output/final_blocks_2.csv`, respectively:

```bash
GenomeDecoder -g test_example/HG38.IHG.Masked.Ns_transformed.fa -g test_example/Orang.IGH.Masked.Ns_transformed.fa -o test_example/output
```

## Contact
For questions or further assistance, please contact **Zhenmiao Zhang** at [zhz142@ucsd.edu](mailto:zhz142@ucsd.edu).
