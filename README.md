# GenomeDecoder
## Installation
Dependencies: [GCC compiler](https://gcc.gnu.org/), [LJA assembler](https://github.com/AntonBankevich/LJA), [Edlib C/C++ library](https://github.com/Martinsos/edlib), and [Cigar python library](https://pypi.org/project/cigar/).

Download and compile: 
```
git clone https://github.com/ZhangZhenmiao/GenomeDecoder.git
cd GenomeDecoder/src && make && cd ..
chmod +x GenomeDecoder consensus_asm edlib parse_cigar.py
```
Add components to $PATH:
```
# replace /path/to/GenomeDecoder with your directory of GenomeDecoder
export PATH=/path/to/GenomeDecoder:$PATH
```

## Usage
GenomeDecoder can take one or multiple genomes as input. The input genomes should be masked by [RepeatMasker](https://github.com/rmhubley/RepeatMasker), and all the Ns be transformed into random nucleotides.

Usage of GenomeDecoder:

```
# for one genome
GenomeDecoder -g <genome1.fa> -o <out_dir>
# for two genomes
GenomeDecoder -g <genome1.fa> -g <genome2.fa> -o <out_dir>
# for three genomes
GenomeDecoder -g <genome1.fa> -g <genome2.fa> -g <genome3.fa> -o <out_dir>
```

Detailed parameters:
```
usage: GenomeDecoder [-h] -g GENOME [-i ITERATIONS] [-k K_VALUES] [-s SIMPLE] [-c COMPLEX] -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Path to the genome
  -i ITERATIONS, --iterations ITERATIONS
                        The number of iterations for three genomes only
  -k K_VALUES, --k_values K_VALUES
                        K-mer values for graph transformations (default 31,51,101,201,401,801,2001)
  -s SIMPLE, --simple SIMPLE
                        Similarity threshold for collapsing simple bulges (default 0)
  -c COMPLEX, --complex COMPLEX
                        Similarity threshold for collapsing complex bulges (default 0.65)
  -o OUTPUT, --output OUTPUT
                        Output directory of GenomeDecoder
```

## Output
In the output directory, `final_blocks_<number>.csv` are the constructed synteny blocks for corresponding input genomes. For example, `final_blocks_1.csv` is the synteny blocks for the first input genome, `final_blocks_2.csv` is the synteny blocks for the second input genome, etc.

## Running example

The command below will generate synteny blocks for `example/HG38.IHG.Masked.Ns_transformed.fa` and `example/Orang.IGH.Masked.Ns_transformed.fa` in files `example/output/final_blocks_1.csv` and `example/output/final_blocks_2.csv`, respectively.

```
GenomeDecoder -g example/HG38.IHG.Masked.Ns_transformed.fa -g example/Orang.IGH.Masked.Ns_transformed.fa -o example/output
```