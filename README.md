
# GenomeDecoder

**GenomeDecoder** is a tool designed for generating synteny blocks in highly repetitive genomic regions. GenomeDecoder can infer synteny blocks in complete genomes; however, for large genomes, the runtime can be substantial.

## Installation

Download GenomeDecoder:
```bash
git clone https://github.com/ZhangZhenmiao/GenomeDecoder.git
```

Dependencies can be installed by creating a new conda environment `genomedecoder`:

```bash
cd GenomeDecoder
conda env create -f requirements.yml
conda activate genomedecoder
```

Build GenomeDecoder:
```bash
# build GenomeDecoder
cd src && make
chmod +x consensus_asm edlib_align parse_cigar.py synteny synteny_analysis.py
cd .. && chmod +x GenomeDecoder

# build LJA
cd src/LJA
cmake .
make jumboDBG
cd ../../
```

## Usage

### Preprocessing
We recommend running RepeatMasker to mask low-complexity repeats. The species name is required (human by default):

```bash
RepeatMasker [-species species_name] input_file.fa
```
This will create `input_file.fa.masked` to be inputted to GenomeDecoder. See more details in [RepeatMasker](https://www.repeatmasker.org/webrepeatmaskerhelp.html) manual.

### Parameters
```bash
GenomeDecoder [-h] -g GENOME [-i ITERATIONS] [-k K_VALUES] [-s SIMPLE] [-c COMPLEX] -o OUTPUT
```
- **-h, --help**  
  Displays the help message and exits.

- **-g GENOME, --genome GENOME**  
  Specifies the path to a genome file to be analyzed. The input genome can be either a complete genome or a genomic subsequence in FASTA format. GenomeDecoder assumes each input contains a single sequence; if the file includes multiple contigs, they will be concatenated using 2,000 “N”s as separators. The `-g` option can be specified multiple times to compare different sequences, or a single time to identify synteny blocks within the sequence.

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

GenomeDecoder generates files named `final_blocks_<i>.csv`. Each `i` corresponds to the synteny blocks of the i<sup>th</sup> input genome, with a minimum block length of 2 kb. It contains 12 columns:

- **Synteny Block in Transformed Sequence**: Provides each block-instance with a unique format: `block_ID(multiplicity_in_disembroiled_graph)length_in_kb`. This is the synteny block.
  
- **Start Vertex in Transformed Graph**: The starting vertex of the block-instance in the disembroiled graph.

- **End Vertex in Transformed Graph**: The ending vertex of the block-instance in the disembroiled graph.

- **Start Pos in Transformed Sequence**: The start coordinate of the block-instance in the disembroiled genome.

- **End Pos in Transformed Sequence**: The end coordinate of the block-instance in the disembroiled genome.

- **Include Left Vertex**: Indicates if the block-instance includes the start vertex (1 for yes, 0 for no).

- **Include Right Vertex**: Indicates if the block-instance includes the end vertex (1 for yes, 0 for no).

- **Is Short Version**: Specifies if the block-instance is a shortened version of the block (1 for yes, 0 for no).

- **Start Pos in Original Sequence**: The start coordinate of the block-instance in the original genome. If the i<sup>th</sup> input genome contains multiple contigs, the coordinate is defined on the concatenated sequence, where contigs are joined by 2,000 “N”s. For example, a coordinate x on the second contig corresponds to (length of the first contig) + 2,000 + x for this column.

- **End Pos in Original Sequence**: The end coordinate of the block instance in the original genome. If the i<sup>th</sup> input genome contains multiple contigs, the coordinate is defined on the concatenated sequence, where contigs are joined by 2,000 “N”s. For example, a coordinate x on the second contig corresponds to (length of the first contig) + 2,000 + x for this column.

- **Similarity of Transformed Block and Original Block**: The percent identity between the block-instance in the disembroiled genome and its corresponding instance in the original genome, calculated by Edlib.

- **Length in Original Sequence**: The length, in base pairs (bp), of this block-instance in the original genome.

If you only need the synteny blocks and their coordinates, focus on the following three columns: **Synteny Block in Transformed Sequence**, **Start Pos in Original Sequence**, and **End Pos in Original Sequence**. The other columns are additional information.

### Transform synteny blocks into alphabets by `synteny_analysis.py`

`synteny_analysis.py` under the `src` can transform CSV outputs into alphabets (A-Z). This script is suitable for fewer than three input genomes. If the number of blocks exceeds 26, the script will produce an error.

#### Usage
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

The following command generates synteny blocks for `test_example/HG38.IHG.Masked.Ns_transformed.fa` and `test_example/Orang.IGH.Masked.Ns_transformed.fa`, saving the synteny blocks to `test_example/output/final_blocks_1.csv` and `test_example/output/final_blocks_2.csv`, respectively:

```bash
GenomeDecoder -g test_example/HG38.IHG.Masked.Ns_transformed.fa -g test_example/Orang.IGH.Masked.Ns_transformed.fa -o test_example/output
```

## Contact
For questions or further assistance, please contact **Zhenmiao Zhang** at [zhz142@ucsd.edu](mailto:zhz142@ucsd.edu).
