from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, CompoundLocation


class Transcript(object):
    def __init__(self, transcript_id, chromosome = None, start = 0, end = 0, strand = None, features = list(), **gff_attributes):
        self.transcript_id = transcript_id
        self.chromosome = chromosome
        self.attributes = gff_attributes
        self.features = features
        
        if strand in ("-1", -1, "-"):
            strand = -1
        else:
            strand = 1
        
        self.location = FeatureLocation(start = start, end = end, strand = strand, ref = transcript_id)
        
        # cached properties
        self._polypeptide = None
        self._exons = None
        self._five_prime_utr = None
        self._three_prime_utr = None
        
        self._nucleic_coding_sequence = None
        self._protein_sequence = None
    
    def __repr__(self):
        return u"%s" % self.transcript_id
    
    def get(self, attribute):
        return self.__dict__.get(attribute, None)
    
    def get_all_ids(self):
        """
        :rtype: set[str]
        """
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
    def protein_sequence(self):
        """
        Return the protein sequence for this transcript using CDS annotations
        :rtype: Bio.Seq.Seq
        """
        if not self._protein_sequence:
            nucleic_sequence = self.nucleic_coding_sequence
            
            if not nucleic_sequence:
                return Seq("")
            
            # make sure the sequence length is a multiple of 3
            end_position = len(nucleic_sequence) - len(nucleic_sequence) % 3
            if isinstance(nucleic_sequence, str):
                nucleic_sequence = Seq(nucleic_sequence[:end_position])
            else:
                nucleic_sequence = nucleic_sequence[:end_position]
            
            self._protein_sequence = nucleic_sequence.translate()
        
        return self._protein_sequence
    
    def nucleic_sequence(self, feature_type):
        """
        :type feature_type: str
        :rtype: Bio.Seq.Seq
        """
        locations = self._get_features(feature_type)
        if locations:
            if len(locations) > 1:
                feature = CompoundLocation(locations)
            else:
                feature = locations[0]
            
            return feature.extract(self.chromosome.nucleic_sequence)
        return None
    
    @property
    def nucleic_coding_sequence(self):
        """
        :rtype: Bio.Seq.Seq
        """
        if not self._nucleic_coding_sequence:
            if self.cds and len(self.cds) > 1:
                cds = CompoundLocation(self.cds)
            else:
                cds = self.cds[0]
            
            self._nucleic_coding_sequence = cds.extract(self.chromosome.nucleic_sequence)
            if self.location.strand == -1:
                self._nucleic_coding_sequence = Seq(self._nucleic_coding_sequence).reverse_complement()
    
        return self._nucleic_coding_sequence
    
    @property
    def cds(self):
        """
        :rtype: list[FeatureLocation]
        """
        if not hasattr(self, "_cds"):
            self._cds = sorted(
                map(
                    lambda exon: FeatureLocation(
                        start = max(exon.start, self.polypeptide.start),
                        end = min(exon.end, self.polypeptide.end)
                    ),
                    filter(
                        lambda exon: exon.start <= self.polypeptide.end and exon.end >= self.polypeptide.start,
                        self.exons
                    ),
                ),
                key = lambda exon: exon.start
            )
        
        return self._cds
    
    @property
    def exons(self):
        """
        :rtype: list[FeatureLocation]
        """
        if self._exons is None:
            self._exons = self._get_features("exon")
        return self._exons
    
    @property
    def introns(self):
        """
        :type: list[FeatureLocation]
        """
        if not hasattr(self, '_introns'):
            self._introns = []
        
            for i in range(len(self.exons) - 1):
                prev_exon, next_exon = self.exons[i:i + 2]
                start = prev_exon.end + 1
                end = next_exon.start - 1
                if abs(end - start) > 1:
                    self._introns.append(FeatureLocation(start = start, end = end))
        
        return self._introns
    
    @property
    def polypeptide(self):
        """
        :rtype: FeatureLocation
        """
        if self._polypeptide is None:
            self._polypeptide = self._get_features("polypeptide")[0]
        
        return self._polypeptide
    
    @property
    def five_prime_utr(self):
        """
        :rtype: list[FeatureLocation]
        """
        if self._five_prime_utr is None:
            self._five_prime_utr = sorted(
                map(
                    lambda exon: FeatureLocation(
                        start = exon.start,
                        end = min(exon.end, self.polypeptide.start)
                    ),
                    filter(
                        lambda exon: exon.start <= self.polypeptide.start,
                        self.exons
                    )
                ),
                key = lambda exon: exon.start
            )
        
        return self._five_prime_utr
    
    @property
    def three_prime_utr(self):
        """
        :rtype: list[FeatureLocation]
        """
        if not hasattr(self, '_three_prime_utr'):
            self._three_prime_utr = sorted(
                map(
                    lambda exon: FeatureLocation(
                        start = max(exon.start, self.polypeptide.end),
                        end = exon.end
                    ),
                    filter(
                        lambda exon: exon.end >= self.polypeptide.end,
                        self.exons
                    )
                ),
                key = lambda exon: exon.start
            )
        
        return self._three_prime_utr
    
    def _get_features(self, feature_type):
        """
        :type feature_type: str
        :rtype: list[FeatureLocation]
        """
        locations = []
        features = sorted(
            filter(lambda f: f['type'].lower() == feature_type, self.features),
            key = lambda f: f.get('start', 0)
        )
        for feature in features:
            locations.append(FeatureLocation(feature['start'], feature['end']))
        
        return locations
