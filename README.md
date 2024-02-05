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
GenomeDecoder can take one or two genomes as input.

For one genome:
```
# replace <genome1.fa>, <out_dir> with your paths to the genome and output directory
GenomeDecoder -1 <genome1.fa> -o <out_dir>
```

For two genomes:
```
# replace <genome1.fa>, <genome2.fa>, <out_dir> with your paths to the two genomes and output directory
GenomeDecoder -1 <genome1.fa> -2 <genome2.fa> -o <out_dir>
```

Detailed parameters:
```
usage: GenomeDecoder [-h] -1 GENOME1 [-2 GENOME2] [-k K_VALUES] [-s SIMPLE] [-c COMPLEX] -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -1 GENOME1, --genome1 GENOME1
                        Path to the first genome
  -2 GENOME2, --genome2 GENOME2
                        Path to the second genome (optional)
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
In the output directory, `final_blocks_1.csv` and `final_blocks_2.csv` (if genome2 is provided) are the constructed synteny blocks. 
