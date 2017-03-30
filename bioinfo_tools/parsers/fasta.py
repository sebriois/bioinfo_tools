from typing import Tuple, Iterable

from Bio import Seq
from Bio import SeqIO


class FastaParser(object):
    def __init__(self):
        pass
    
    def read(self, fasta_file: str) -> Iterable[Tuple[str,str]]:
        with open(fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                seqid = record.id.split("|")[-1]
                seq = str(record.seq)
                yield (seqid, seq)
