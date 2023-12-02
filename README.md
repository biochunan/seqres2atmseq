# seqres2atmseq
Align the SEQRES sequence with the ATMSEQ sequence of a protein chain and output a mask array (0: a missing residue in the ATMSEQ; 1: match)

## Table of Contents
- [seqres2atmseq](#seqres2atmseq)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
  - [Usage](#usage)
    - [Example: use fasta file as input](#example-use-fasta-file-as-input)
      - [FASTA file format](#fasta-file-format)
    - [Example: use SEQRES and ATMSEQ as input](#example-use-seqres-and-atmseq-as-input)

---

## Installation 
If you are using a conda environment, you can install seqres2atmseq by:
```
(base) $ git clone git@github.com:biochunan/seqres2atmseq.git
(base) $ cd seqres2atmseq 
(base) $ pip install . 
```
This will install a command-line tool `seqres2atmseq` in your conda environment.

## Usage
Check help message:
```shell 
(base) $ seqres2atmseq -h  # see below 
```

Help message:
```
usage: seqres2atmseq [-h] [-i FASTA] [-s SEQRES] [-a ATMSEQ] [-n SEQ_NAME] [-c CLUSTAL_OMEGA_EXECUTABLE] [-o OUTPUT] [-v]

Process sequence file.

options:
  -h, --help            show this help message and exit
  -i FASTA, --fasta FASTA
                        Input FASTA file path
  -s SEQRES, --seqres SEQRES
                        SEQRES sequence
  -a ATMSEQ, --atmseq ATMSEQ
                        ATMSEQ sequence
  -n SEQ_NAME, --seq_name SEQ_NAME
                        Sequence name
  -c CLUSTAL_OMEGA_EXECUTABLE, --clustal_omega_executable CLUSTAL_OMEGA_EXECUTABLE
                        Path to clustal omega executable
  -o OUTPUT, --output OUTPUT
                        Output directory or file path
  -v, --verbose         Verbose mode
```

### Example: use fasta file as input 
Run alignment and save the mask json file with a FASTA file as input. We provie a test FASTA file `seq.fasta` in the `test` directory.
Refer to [FASTA file format](#fasta-file-format) for the FASTA file format.

```shell
# current directory: seqres2atmseq
(base) $ seqres2atmseq -i test/test.fasta -o test/mask.json --verbose 
```
- `-i`: input FASTA file path
- `-o`: output file path
  - if the output file path is a directory, the output file will be saved in the directory with name `mask.json`
  - if the output file path is a file, the output file will be saved as the file
- `--verbose`: verbose mode, this will print the alignment result in the terminal

<details>
<summary>stdout:</summary>

```
2023-12-02 21:55:25.118 | DEBUG    | seqres2atmseq.app:main:275 - 
chocolate_A seqres: ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIMRGSSGEGVSCTIRSSLLGLEKTASISIADPFFRSAQ
chocolate_A atmseq: ---QFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIM------GVSCTIRSSLLGLEKTASISIADPFF----
chocolate_A mask  : 0001111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111000000111111111111111111111111110000

2023-12-02 21:55:25.131 | DEBUG    | seqres2atmseq.app:main:275 - 
sunflower_B seqres: ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIMRGSSGEGVSCTIRSSLLGLEKTASISIADPFFRSAQ
sunflower_B atmseq: ---QFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIM------GVSCTIRSSLLGLEKTASISIADPFF----
sunflower_B mask  : 0001111111111111111111111111111111111111111111111111111111111111111111111111111
```
</details>


#### FASTA file format
The FASTA file should be in the following format:
```
>chocolate_A|seqres
ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSS
>chocolate_A|atmseq
QFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSS
>sunflower_B|atmseq
ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSS
>sunflower_B|seqres
ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSS
```
Each protein chain should have a pair of SEQRES and ATMSEQ sequences, as indicated by the suffix `|seqres` and `|atmseq` in the sequence header above. 
Each protein chain is distinguished by the prefix of the sequence header, e.g. `chocolate_A` and `sunflower_B` in the example above. This can be any string, but should be unique for each protein chain. For example, you can use PDB ID and chain id `[pdbID]_[chainID]` e.g. `ABCD_A` as the prefix.

Refer to the example FASTA file `test/test.fasta` in the `test` directory.


### Example: use SEQRES and ATMSEQ as input
Run alignment and save the mask json file with SEQRES and ATMSEQ as input. 
```shell   
# current directory: seqres2atmseq
(base) $ seqres2atmseq \
  -s \
  ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIMRGSSGEGVSCTIRSSLLGLEKTASISIADPFFRSAQ \
  -a \ 
  QFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIMGVSCTIRSSLLGLEKTASISIADPFF \
  -o test/mask.json \
  --verbose
```

stdout: 
```shell 
2023-12-02 21:50:17.716 | DEBUG    | seqres2atmseq.app:main:275 - 
seq seqres: ADLQFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIMRGSSGEGVSCTIRSSLLGLEKTASISIADPFFRSAQ
seq atmseq: ---QFSVLGPSGPILAMVGEDADLPCHLFPTMSAETMELKWVSSSLRQVVNVYADGKEVEDRQSAPYRGRTSILRDGITAGKAALRIHNVTASDSGKYLCYFQDGDFYEKALVELKVAALGSDLHVDVKGYKDGGIHLECRSTGWYPQPQIQWSNNKGENIPTVEAPVVADGVGLYAVAASVIM------GVSCTIRSSLLGLEKTASISIADPFF----
seq mask  : 0001111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111000000111111111111111111111111110000
```
