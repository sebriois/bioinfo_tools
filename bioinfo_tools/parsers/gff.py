import copy
import re


class Gff3(object):
    def __init__(self):
        pass
    
    def read(self, fh):
        """
        :param fh: input file handler
        :return: iterator that yields gene features
        """
        gene = None
        mRNA = None
        
        feature = self._readline(fh)
        while feature:
            if feature['type'] == 'gene':
                if gene is not None:
                    if mRNA is not None:
                        if 'mRNA' not in gene:
                            gene['mRNA'] = []
                        gene['mRNA'].append(mRNA)
                    yield gene
                    mRNA = None
                
                gene = feature
            
            elif feature['type'] in ('mRNA', 'mrna'):
                if mRNA is not None:
                    if 'mRNA' not in gene:
                        gene['mRNA'] = []
                    gene['mRNA'].append(mRNA)
                mRNA = feature
            
            else:
                if 'features' not in mRNA:
                    mRNA['features'] = []
                mRNA['features'].append(feature)
            
            feature = self._readline(fh)
        
        if gene is not None:
            if mRNA is not None:
                if 'mRNA' not in gene:
                    gene['mRNA'] = []
                gene['mRNA'].append(mRNA)
            yield gene
    
    def parse(self, fh):
        genes = []
        for gene in self.read(fh):
            genes.append(gene)
        return genes

    def _readline(self, fh):
        """
        Parse the next line of the file handler
        :rtype: dict
        """
        try:
            line = next(fh)
        except StopIteration:
            return None
        
        if isinstance(line, bytes):
            line = line.decode()

        line = line.strip()
        if not line:
            return self._readline(fh)
        
        if line.startswith("#"):
            return self._readline(fh)
        
        cols = re.split("\s+", line, 8)
        if len(cols) < 9:
            cols = re.split("\t", line, 8)
        
        if len(cols) < 9:
            raise Exception("Missing attributes column in line: ['%s']" % "','".join(cols))

        return {
            'seqid': cols[0],
            'source': cols[1],
            'type': cols[2],
            'start': int(cols[3]) - 1,
            'end': int(cols[4]),
            'score': cols[5] in ('.', '') and '.' or float(cols[5]),
            'strand': cols[6],
            'phase': cols[7] in ('.', '') and '.' or int(cols[7]),
            'attributes': self._parse_attributes(cols[8])
        }

    def _parse_attributes(self, gff_column):
        """
        Parse column #9 of a gff file
        """
        attributes = {}
    
        for attribute in re.split(';', gff_column):
            # handles cases such as ID=Name=toto
            keys = attribute.split("=")[:-1]
            value = attribute.split("=")[-1]
        
            for key in keys:
                attributes[key] = value
    
        return attributes


class Gff(object):
    EXPECTED_TYPES = {'exon', 'intron', 'cds', 'polypeptide', 'five_prime_utr', 'three_prime_utr'}
    
    def __init__(self, fh):
        self.fh = fh
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """
        :return: next gene feature (containing mRNA, CDS, polypeptide, exons)
        """
        line = self._readline()
        while line and line['type'] != 'gene':
            line = self._readline()
        
        if not line:
            raise StopIteration
        
        gene = copy.deepcopy(line)
        
        current_mRNA = None
        
        line = self._readline()
        while line and line['type'].lower() != 'gene':
            if line['type'].lower() == 'mrna':
                if current_mRNA:
                    if 'mRNA' not in gene:
                        gene['mRNA'] = []
                    
                    # self.add_missing_types(current_mRNA)
                    gene['mRNA'].append(current_mRNA)
                
                current_mRNA = copy.deepcopy(line)
            
            elif line['type'].lower() in self.EXPECTED_TYPES:
                if 'features' not in current_mRNA:
                    current_mRNA['features'] = []
                current_mRNA['features'].append(line)
            
            line = self._readline()
        
        if current_mRNA:
            if 'mRNA' not in gene:
                gene['mRNA'] = []
            # self.add_missing_types(current_mRNA)
            gene['mRNA'].append(current_mRNA)
        
        if line:
            self.fh.seek(self.last_fh_position)
        
        return gene
    
    def _readline(self):
        """
        Parse a given line
        """
        self.last_fh_position = self.fh.tell()
        
        # read next line
        line = self.fh.readline()
        line = line.strip()
        
        # skip lines starting with '#'
        while line and line.startswith('#'):
            line = self.fh.readline()
            line = line.strip()
        
        # end of file ?
        if not line:
            return {}
        
        cols = re.split("\s+", line, 8)
        if len(cols) < 9:
            cols = re.split("\t", line, 8)
        
        if len(cols) < 9:
            raise Exception("Missing attributes column in line: ['%s']" % "','".join(cols))
        
        return {
            'seqid': cols[0],
            'source': cols[1],
            'type': cols[2],
            'start': int(cols[3]) - 1,
            'end': int(cols[4]),
            'score': cols[5] in ('.', '') and '.' or float(cols[5]),
            'strand': cols[6],
            'phase': cols[7] in ('.', '') and '.' or int(cols[7]),
            'attributes': self._parse_attributes(cols[8])
        }
    
    def _parse_attributes(self, gff_column):
        """
        Parse column #9 of a gff file
        """
        attributes = {}
        
        for attribute in re.split(';', gff_column):
            # handles cases such as ID=Name=toto
            keys = attribute.split("=")[:-1]
            value = attribute.split("=")[-1]
            
            for key in keys:
                attributes[key] = value
        
        return attributes
    
    def add_missing_types(self, mRNA):
        features = mRNA.get('features', [])
        if not features:
            return
        
        existing_types = {feature['type'].lower() for feature in features}
        if existing_types == self.EXPECTED_TYPES:
            return
        
        missing_features = []
        
        if 'exon' in existing_types and 'polypeptide' in existing_types:
            # get the first polypeptide
            polypeptide = next(filter(lambda feature: feature['type'] == 'polypeptide', features))
            
            # sort exons on their positions
            exons = sorted(
                filter(lambda feature: feature['type'] == 'exon', features),
                key = lambda e: min(e['start'], e['end']),
                reverse = (mRNA['strand'] == '-')
            )
            
            # if 'intron' not in existing_types:
            #     missing_features.extend(self.get_introns(exons))
            if 'cds' not in existing_types:
                missing_features.extend(self.infer_cds(exons, polypeptide))
            if 'five_prime_utr' not in existing_types:
                missing_features.extend(self.infer_five_prime_utr(exons, polypeptide, mRNA['strand']))
            if 'three_prime_utr' not in existing_types:
                missing_features.extend(self.infer_three_prime_utr(exons, polypeptide, mRNA['strand']))
        
        if len(missing_features) > 0:
            for missing_feature in missing_features:
                mRNA['features'].append(missing_feature)
    
    def infer_five_prime_utr(self, exons, polypeptide, strand):
        five_prime_utr = []
        
        if strand == '+':
            for exon in filter(lambda e: e['start'] < polypeptide['start'], exons):
                start = exon['start']
                end = min(exon['end'], polypeptide['start'] - 1)
                
                if abs(end - start) > 1:
                    five_prime_utr.append({'type': 'five_prime_utr', 'start': start, 'end': end})
        
        elif strand == '-':
            for exon in filter(lambda e: e['end'] > polypeptide['end'], exons):
                start = max(exon['start'], polypeptide['end'] + 1)
                end = exon['end']
                
                if abs(end - start) > 1:
                    five_prime_utr.append({'type': 'five_prime_utr', 'start': start, 'end': end})
        
        return five_prime_utr
    
    def infer_three_prime_utr(self, exons, polypeptide, strand):
        three_prime_utr = []
        
        if strand == '+':
            for exon in filter(lambda e: e['end'] > polypeptide['end'], exons):
                start = max(exon['start'], polypeptide['end'] + 1)
                end = exon['end']
                
                if abs(end - start) > 1:
                    three_prime_utr.append({'type': 'three_prime_utr', 'start': start, 'end': end})
        
        elif strand == '-':
            for exon in filter(lambda e: e['start'] < polypeptide['start'], exons):
                start = exon['start']
                end = min(exon['end'], polypeptide['start'] - 1)
                
                if abs(end - start) > 1:
                    three_prime_utr.append({'type': 'three_prime_utr', 'start': start, 'end': end})
        
        return three_prime_utr
    
    def infer_cds(self, exons, polypeptide):
        cds = []
        
        cds_phase = 0
        for exon in filter(lambda e: e['start'] < polypeptide['end'] and e['end'] > polypeptide['start'], exons):
            start = max(exon['start'], polypeptide['start'])
            end = min(exon['end'], polypeptide['end'])
            
            cds.append({
                'type': 'CDS',
                'start': start,
                'end': end,
                'phase': cds_phase
            })
            
            cds_length = (end - start) + 1
            cds_phase = (3 - ((cds_length - cds_phase) % 3)) % 3  # last %3 is to avoid phase = 3 rather than phase = 0
        
        return cds
    
    def get_introns(self, exons):
        introns = []
        
        for i in range(len(exons) - 1):
            prev_exon, next_exon = exons[i:i + 2]
            start = prev_exon['end'] + 1
            end = next_exon['start'] - 1
            if abs(end - start) > 1:
                introns.append({'type': 'intron', 'start': start, 'end': end})
        
        return introns
