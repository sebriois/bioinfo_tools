import binascii
import gzip

import os

from bioinfo_tools.parsers.fasta import FastaParser

from bioinfo_tools.genomic_features.chromosome import Chromosome

from bioinfo_tools.parsers.gff import Gff3


class Genome(object):
    def __init__(self, name):
        self.name = name
        self.chromosomes = dict()  # <chromosome_id>: <chromosome.Chromosome>
        self.protein_sequences = dict()  # <protein_id>: sequence
    
    def get_gene(self, gene_id):
        """
        :param gene_id: a gene ID
        :type gene_id: str
        
        :rtype: gene.Gene
        :return: a Gene object
        """
        if not self.chromosomes:
            raise Exception("[%s] gff file not loaded, no gene can be found" % self.name)
        
        return next(
            map(
                lambda chromosome: chromosome.get_gene(gene_id),
                filter(
                    lambda chromosome: gene_id in chromosome._genes,
                    self.chromosomes.values()
                )
            ),
            None
        )
    
    def load_gff(self, filepath):
        """
        Parse the given GFF file and load its features in memory
        :type filepath: str
        :rtype: int
        :return: number of loaded genes
        """
        
        # is filepath a gzip file?
        with open(filepath, 'rb') as test_f:
            is_gzip_file = binascii.hexlify(test_f.read(2)) == b'1f8b'

        if is_gzip_file:
            fh = gzip.open(filepath, "rb")
        else:
            fh = open(filepath, "r")
        
        gene_count = 0
        gff_parser = Gff3()
        for gff_feature in gff_parser.read(fh):
            gene_count += 1
            chromosome_id = gff_feature['seqid']
            
            if chromosome_id not in self.chromosomes:
                chromosome = Chromosome(chromosome_id, self.name)
                self.chromosomes[chromosome_id] = chromosome
            
            self.chromosomes[chromosome_id].add_gene(gff_feature)

        fh.close()
        
        return gene_count
    
    def load_fasta_genome(self, filepath):
        fasta_parser = FastaParser()
        for chromosome_id, nucleic_sequence in fasta_parser.read(filepath):
            if chromosome_id not in self.chromosomes:
                self.chromosomes[chromosome_id] = Chromosome(chromosome_id, self.name)
            
            self.chromosomes[chromosome_id].attach_nucleic_sequence(nucleic_sequence)
    
    def load_fasta_pept(self, filepath):
        fasta_parser = FastaParser()
        for (seqid, sequence) in fasta_parser.read(filepath):
            seqid = seqid.replace("%s|" % os.path.basename(filepath), "")  # TODO: this does not apply widely
            self.protein_sequences[seqid.split("|")[0]] = sequence
    