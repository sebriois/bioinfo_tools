# bioinfo_tools 0.3.1

## Installation

```bash
pip install bioinfo_tools
```

## Parsers

*HEADS UP!* These parsers are still under development and usage is not consistent from one parser to another.

### Fasta parser

```python
from bioinfo_tools.parsers.fasta import FastaParser

fasta_parser = FastaParser()

# by default, sequence IDs are separated by the firstly found '|' or ':'
for seqid, sequence in fasta_parser.read("/path/to/file.fasta"):
    print(seqid, sequence)

# you may specify a specific separator for your sequence ID (e.g white space):
for seqid, sequence in fasta_parser.read("/path/to/file.fasta", id_separator=" "):
    print(seqid, sequence)
```

### GFF parser

```python
from bioinfo_tools.parsers.gff import Gff3

gff_parser = Gff3()
with open("/path/to/file.gff", "r") as fh:
    for gene in gff_parser.read(fh):
        print(gene)

import gzip
with gzip.open("/path/to/file.gz", "rb") as fh:
    for gene in gff_parser.read(fh):
        print(gene)
```

### OBO parser


```python
from bioinfo_tools.parsers.obo import OboParser

obo_parser = OboParser()
with open("/path/to/file.obo") as fh:
    go_terms = obo_parser.read(fh)

for go_term in go_terms.values():
    print(go_term)

    # you may also get the GO term parents via the parser
    parents = obo_parser.get_parents(go_term)
```

## Usage Examples

### Extract all introns sequences by parsing GFF and fasta files

In this example, we focus on a genome assembly. We will first load a GFF file containing gene annotations for this
assembly, then load a fastA file containing the nucleic sequences of each chromosome in the genome.
We will then collect all transcript introns and extract their nucleic sequences.

**__DISCLAIMER__**: for this example to work, your GFF file must expose at least the following feature types in column #3:
 - `gene`
 - one of `transcript|mRNA|RNA` (or lowercased version)


```python
from bioinfo_tools.genomic_features.chromosome import Chromosome
from bioinfo_tools.parsers.gff import Gff3
from bioinfo_tools.parsers.fasta import FastaParser

chromosomes = dict()  # {<chromosome_id>: <bioinfo_tools.genomic_features.Chromosome>}

# start with parsing a GFF file
gff_parser = Gff3()
with open("/path/to/gene_models.gff", "r") as fh:
    for gene in gff_parser.read(fh):
        chromosome = gene['seqid']

        if chromosome not in chromosomes:
            chromosomes[chromosome] = Chromosome(chromosome)  # init a new Chromosome object

        chromosomes[chromosome].add_gene(gene)  # add the current gene to our Chromosome object

# load our chromosome sequences in memory
fasta_parser = FastaParser()
for chromosome, nucleic_sequence in fasta_parser.read("/path/to/genome_chromosomes.fasta"):
    if chromosome not in chromosomes:
        chromosomes[chromosome] = Chromosome(chromosome)
    # attach parsed chromosome sequence to our Chromosome object
    chromosomes[chromosome].attach_nucleic_sequence(nucleic_sequence)

# now, collect introns and extact their nucleic sequence
introns_sequences = dict()  #  {<intron_id>: <intron_sequence>}
for chromosome in chromosomes.values():
    for gene in chromosome.genes:
        for transcript in gene.transcripts:
            for idx, intron in enumerate(transcript.introns):
                intron_id = "%s_intron_%s" % (transcript.transcript_id, idx)
                intron_seq = intron.extract(chromosome.nucleic_sequence)  # that we attached above
                introns_sequences[intron_id] = intron_seq

# from here, you can do what you want with the intron sequences (eg. write them to a fasta file, etc)
# ...
```

__Note:__ when at the transcript level, you can grab its feature types as described in your GFF file by doing so:
```python
for feature in transcript._get_features("exon"):
    print(feature)  # I'm an exon
```
For convenience and clarity, following properties are available on transcript objects:
```python
print(transcript.introns)  # will call transcript._get_features('intron') behind the scenes
print(transcript.exons)  # will call transcript._get_features('exon') behind the scenes
print(transcript.cds)  # will call transcript._get_features('cds') behind the scenes
print(transcript.polypeptide)  # will call transcript._get_features('polypeptide') behind the scenes
print(transcript.five_prime_utr)  # will call transcript._get_features('five_prime_utr') behind the scenes
print(transcript.three_prime_utr)  # will call transcript._get_features('three_prime_utr') behind the scenes
```