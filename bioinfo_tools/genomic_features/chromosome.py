from typing import Dict, List

from Bio.SeqFeature import FeatureLocation

from bioinfo_tools.genomic_features.gene import Gene


class Chromosome(object):
    def __init__(self, chromosome_id, assembly_name = None):
        self.chromosome_id = chromosome_id
        self.assembly_name = assembly_name
        self.length = 0
        self._genes = {}
    
    def __repr__(self):
        return u"%s (%s genes)" % (self.chromosome_id, len(self.genes))
    
    def add_gene(self, gff_feature:Dict):
        gene_id = None
        
        if 'gene_id' in gff_feature:
            gene_id = gff_feature.pop('gene_id')
        elif 'attributes' in gff_feature:
            if 'gene_id' in gff_feature['attributes']:
                gene_id = gff_feature['attributes'].pop('gene_id')
            elif 'Name' in gff_feature['attributes']:
                gene_id = gff_feature['attributes']['Name']
            elif 'ID' in gff_feature['attributes']:
                gene_id = gff_feature['attributes']['ID']
        
        if not gene_id:
            raise Exception("gene_id not found in given GFF feature")
        
        # remove all potential duplicated keys
        for key_name in ('chromosome', 'start', 'end', 'strand', 'assembly_name'):
            gff_feature.get('attributes', {}).pop(key_name, None)
        
        gene = Gene(
            gene_id = gene_id,
            chromosome = self,
            start = gff_feature.get('start', 0),
            end = gff_feature.get('end', 0),
            strand = gff_feature.get('strand', None),
            assembly_name = self.assembly_name,
            **gff_feature.get('attributes', {})
        )
        for transcript in gff_feature.get('mRNA', []):
            gene.add_transcript(transcript)
        
        self._genes[gene['gene_id']] = gene
        return gene
    
    @property
    def genes(self) -> List[Gene]:
        return list(self._genes.values())
    
    @property
    def sorted_genes(self) -> List[Gene]:
        if not hasattr(self, "_sorted_genes"):
            self._sorted_genes = sorted(self._genes.values(), key = lambda gene: gene.location.start)
        return self._sorted_genes
    
    def get_gene(self, gene_id) -> Gene:
        return self._genes.get(gene_id, None)

    def attach_nucleic_sequence(self, sequence):
        self.nucleic_sequence = sequence
        self.length = len(sequence)
    
    def extract_sequence(self, start, end):
        location = FeatureLocation(start, end)
        return location.extract(self.nucleic_sequence)
