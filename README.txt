# bioinfo_tools 0.2

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
with open("/path/to/file.gff") as fh:
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