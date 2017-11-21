from typing import Tuple, Iterable

import re
from Bio import SeqIO


class FastaParser(object):
    def __init__(self):
        pass
    
    def read(self, fasta_file: str, id_separator = "\||\:", all_after_sep = True) -> Iterable[Tuple[str,str]]:
        with open(fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                try:
                    if all_after_sep:
                        seqid = re.split(id_separator, record.id, 1)[1]
                    else:
                        seqid = re.split(id_separator, record.id, 1)[-1]
                except IndexError:
                    seqid = record.id
                seq = str(record.seq)
                yield (seqid, seq)
