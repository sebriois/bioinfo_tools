from typing import Tuple, Iterable

import re
from Bio import SeqIO


class FastaParser(object):
    def __init__(self):
        pass
    
    def read(self, fasta_file: str) -> Iterable[Tuple[str,str]]:
        with open(fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                seqid = re.split("\||\:", record.id, 1)[1]
                seq = str(record.seq)
                yield (seqid, seq)
