from Bio.SeqFeature import FeatureLocation

from bioinfo_tools.genomic_features.gene import Gene


class Chromosome(object):
    def __init__(self, chromosome_id):
        self.chromosome_id = chromosome_id
        self.length = 0
        self._genes = {}
    
    def __repr__(self):
        return u"%s (%s genes)" % (self.chromosome_id, len(self.genes))
    
    def add_gene(self, gff_feature):
        """
        :type gff_feature: str
        :rtype: None
        """
        gene = Gene(self, gff_feature)
        self._genes[gene['gene_id']] = gene
    
    @property
    def genes(self):
        return list(self._genes.values())
    
    @property
    def sorted_genes(self):
        if not hasattr(self, "_sorted_genes"):
            self._sorted_genes = sorted(self._genes.values(), key = lambda gene: gene.location.start)
        return self._sorted_genes
    
    def get_gene(self, gene_id):
        return self._genes.get(gene_id, None)

    def attach_nucleic_sequence(self, sequence):
        self.nucleic_sequence = sequence
        self.length = len(sequence)
    
    def extract_sequence(self, start, end):
        location = FeatureLocation(start, end)
        return location.extract(self.nucleic_sequence)
