#!/usr/bin/env python
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from bioinfo_tools.genomic_features.transcript import Transcript


class Gene(object):
    def __init__(self, gene_id, chromosome = None, start = 0, end = 0, strand = None, assembly_name = None, **gff_attributes):
        self.gene_id = gene_id
        self.chromosome = chromosome
        self.attributes = gff_attributes
        self.transcripts = list()
        self.assembly_name = assembly_name
        
        if strand in ("-1", -1, "-"):
            self.strand = -1
        else:
            self.strand = 1
        
        self.location = FeatureLocation(start = start, end = end, strand = self.strand, ref = gene_id)
    
    def __repr__(self):
        return u"%s" % self.gene_id
    
    def __getitem__(self, item):
        return getattr(self, item)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __eq__(self, other):
        if isinstance(other, str):
            return self.gene_id == other
        elif isinstance(other, Gene):
            return self.gene_id == other.gene_id
        return self == other
    
    def get(self, *args):
        return self.__dict__.get(*args)
    
    def add_transcript(self, mRNA_feature, lookup_attributes = ('transcript_id', 'Name', 'ID', 'id')):
        """
        :param mRNA_feature: dictionnary of various gff features describing the transcript
        :type mRNA_feature: dict
        
        :param lookup_attributes: potential attribute names hosting the gene ID
        :type lookup_attributes: tuple
        
        :return: None
        """
        transcript_id = None
        
        if 'transcript_id' in mRNA_feature:
            transcript_id = mRNA_feature.pop('transcript_id')
        elif 'attributes' in mRNA_feature:
            for attribute in lookup_attributes:
                if attribute in mRNA_feature['attributes'] and mRNA_feature['attributes'][attribute].rstrip():
                    transcript_id = mRNA_feature['attributes'][attribute]
                    break

        # remove all potential duplicated keys
        for key_name in ('chromosome', 'start', 'end', 'strand', 'transcript_id'):
            mRNA_feature.get('attributes', {}).pop(key_name, None)

        transcript = Transcript(
            transcript_id = transcript_id,
            chromosome = self.chromosome,
            start = mRNA_feature.get('start', 0),
            end = mRNA_feature.get('end', 0),
            strand = mRNA_feature.get('strand', None),
            features = mRNA_feature.get('features', []),
            **mRNA_feature.get('attributes', {})
        )
        self.transcripts.append(transcript)
    
    def as_fasta(self, **kwargs):
        record = SeqRecord(self.extract_sequence(**kwargs), id = self.gene_id)
        return record.format("fasta")
    
    def extract_sequence(self, upstream = 0, downstream = 0):
        location = FeatureLocation(self.location.start - upstream, self.location.end + downstream, self.strand)
        return location.extract(self.chromosome.nucleic_sequence)

    def get_all_ids(self):
        """
        return all possible IDs for that gene
        :rtype: list[str]
        """
        if not hasattr(self, '_all_ids'):
            all_ids = set()  # all possible IDs for the given gene
            
            gene_name = self.attributes.get('Name', None)
            if gene_name:
                all_ids.add(gene_name)
            
            gene_id = self.attributes.get('ID', None)
            if gene_id:
                all_ids.add(gene_id)
    
            ancestor_identifier = self.attributes.get('ancestorIdentifier', None)
            if ancestor_identifier:
                all_ids.add(ancestor_identifier)
    
            for transcript in self.transcripts:
                all_ids.update(transcript.get_all_ids())
        
            self._all_ids = list(all_ids)
        
        return self._all_ids
