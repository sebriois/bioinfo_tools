#!/usr/bin/env python
from typing import List, Set

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord


class Gene(object):
    def __init__(self, chromosome, gff_feature):
        self.chromosome = chromosome
        self._gff_feature = gff_feature
        self._transcripts = []
        self._all_ids = []
        
        if 'gene_id' in gff_feature:
            self.gene_id = gff_feature['gene_id']
        else:
            self.gene_id = gff_feature['attributes'].get('Name', gff_feature['attributes']['ID'])
        
        strand = gff_feature.get('strand', 0)
        if strand in ("-1", -1, "-"):
            strand = -1
        elif strand in ("+1", "1", 1, "+"):
            strand = 1
        else:
            strand = 0
        
        self.location = FeatureLocation(
            start = gff_feature['start'],
            end = gff_feature['end'],
            strand = strand
        )
        self.attributes = gff_feature.get('attributes', [])
    
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
    
    @property
    def transcripts(self) -> List:
        if not self._transcripts:
            for mRNA in self._gff_feature.get('mRNA', []):
                self._transcripts.append(Transcript(self.chromosome, mRNA))
        return self._transcripts
    
    def as_fasta(self, **kwargs):
        record = SeqRecord(self.extract_sequence(**kwargs), id = self.gene_id)
        return record.format("fasta")
    
    def extract_sequence(self, upstream = 0, downstream = 0):
        location = FeatureLocation(self.location.start - upstream, self.location.end + downstream, self.location.strand)
        return location.extract(self.chromosome.nucleic_sequence)

    def get_all_ids(self) -> List[str]:
        """
        return all possible IDs for that gene
        :rtype: set
        """
        if not self._all_ids:
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


class Transcript(object):
    def __init__(self, chromosome, gff_feature):
        self.chromosome = chromosome
        self.attributes = gff_feature.get('attributes', {})
        self.features = gff_feature.get('features', [])
        
        if 'transcript_id' in gff_feature:
            self.transcript_id = gff_feature['transcript_id']
        else:
            self.transcript_id = gff_feature['attributes'].get('Name', gff_feature['attributes']['ID'])
        
        strand = gff_feature.get('strand', 0)
        if strand in ("-1", -1, "-"):
            strand = -1
        elif strand in ("+1", "1", 1, "+"):
            strand = 1
        else:
            strand = 0
        
        self.location = FeatureLocation(
            start = gff_feature['start'],
            end = gff_feature['end'],
            strand = strand
        )
        
        # cached properties
        self._polypeptide = None
        self._polypeptide_was_updated = False
        self._exons = None
        self._cds = None
        self._five_prime_utr = None
        self._three_prime_utr = None
        
        self._nucleic_coding_sequence = None
        self._protein_sequence = None
    
    def __repr__(self):
        return u"%s" % self.transcript_id
    
    def get(self, attribute):
        return self.__dict__.get(attribute, None)
    
    def get_all_ids(self) -> Set[str]:
        all_ids = set()  # all possible IDs for the given gene
        
        transcript_name = self.attributes.get('Name', None)
        if transcript_name:
            all_ids.add(transcript_name)
        
        transcript_id = self.attributes.get('ID', None)
        if transcript_id:
            all_ids.add(transcript_id)
        
        ancestor_identifier = self.attributes.get('ancestorIdentifier', None)
        if ancestor_identifier:
            all_ids.add(ancestor_identifier)
        
        if self.protein_id:
            all_ids.add(self.protein_id)
        
        return all_ids
    
    @property
    def protein_id(self):
        try:
            polypeptide_feature = next(filter(lambda f: f['type'] == 'polypeptide', self.features))
            return polypeptide_feature['attributes'].get('Name', polypeptide_feature['attributes'].get('ID', None))
        except StopIteration:
            return None
    
    @property
    def protein_sequence(self) -> Seq:
        """
        Return the protein sequence for this transcript using CDS annotations
        :type chromosome_nucleic_sequence: str
        :rtype: str
        """
        if not self._protein_sequence:
            nucleic_sequence = self.nucleic_coding_sequence
            
            if not nucleic_sequence:
                return None
            
            # make sure the sequence length is a multiple of 3
            end_position = len(nucleic_sequence) - len(nucleic_sequence) % 3
            if isinstance(nucleic_sequence, str):
                nucleic_sequence = Seq(nucleic_sequence[:end_position])
            else:
                nucleic_sequence = nucleic_sequence[:end_position]
            
            try:
                self._protein_sequence = nucleic_sequence.translate()
            except TypeError:
                print("Received weird nucleic_sequence type: %s - %s" % (type(nucleic_sequence), nucleic_sequence))
        
        return self._protein_sequence
    
    def nucleic_sequence(self, feature_type):
        """
        :type feature_type: str
        :rtype: Bio.Seq.Seq
        """
        location = self._get_features(feature_type)
        if location:
            return location.extract(self.chromosome.nucleic_sequence)
        return None
    
    @property
    def nucleic_coding_sequence(self) -> Seq:
        if not self._nucleic_coding_sequence:
            if self.cds:
                self._nucleic_coding_sequence = self.cds.extract(self.chromosome.nucleic_sequence)
                if self.location.strand == -1:
                    self._nucleic_coding_sequence = Seq(self._nucleic_coding_sequence).reverse_complement()
        
        return self._nucleic_coding_sequence
    
    @property
    def cds(self) -> CompoundLocation:
        if self._cds is None:
            cds_list = sorted(
                map(
                    lambda exon: FeatureLocation(
                        start = max(exon.start, self.polypeptide.start),
                        end = min(exon.end, self.polypeptide.end)
                    ),
                    filter(
                        lambda exon: exon.start <= self.polypeptide.end and exon.end >= self.polypeptide.start,
                        self.exons.parts
                    ),
                ),
                key = lambda exon: exon.start
            )
            
            if len(cds_list) == 1:
                self._cds = cds_list[0]
            
            if len(cds_list) > 1:
                self._cds = CompoundLocation(cds_list)
        
        return self._cds
    
    @property
    def exons(self):
        if self._exons is None:
            self._exons = self._get_features("exon")
        return self._exons
    
    @property
    def polypeptide(self):
        if self._polypeptide is None:
            self._polypeptide = self._get_features("polypeptide")
        
        return self._polypeptide
    
    @property
    def five_prime_utr(self):
        if self._five_prime_utr is None:
            utr_list = sorted(
                map(
                    lambda exon: FeatureLocation(
                        start = exon.start,
                        end = min(exon.end, self.polypeptide.start)
                    ),
                    filter(
                        lambda exon: exon.start <= self.polypeptide.start,
                        self.exons.parts
                    )
                ),
                key = lambda exon: exon.start
            )
            
            if len(utr_list) == 1:
                self._five_prime_utr = utr_list[0]
            
            elif len(utr_list) > 1:
                self._five_prime_utr = CompoundLocation(utr_list)

        return self._five_prime_utr
    
    @property
    def three_prime_utr(self):
        return self._get_features("three_prime_utr")
    
    @property
    def utr(self, location):
        if int(location) not in (3, 5):
            raise Exception("location must be 3 or 5")
        
        if int(location) == 3:
            return self.three_prime_utr
        if int(location) == 5:
            return self.five_prime_utr
    
    def _get_features(self, feature_type):
        """
        :type feature_type: str
        :rtype: CompoundLocation
        """
        locations = []
        features = sorted(
            filter(lambda f: f['type'].lower() == feature_type, self.features),
            key = lambda f: f['start']
        )
        for feature in features:
            locations.append(FeatureLocation(feature['start'], feature['end']))
        
        if len(locations) == 1:
            return locations[0]
        if len(locations) > 1:
            return CompoundLocation(locations)
        
        return None
