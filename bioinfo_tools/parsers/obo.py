import copy
import re

from bioinfo_tools.utils.log import Log


class OboParser(Log):
    DEFAULT_KEPT_FIELDS = ['id', 'name', 'namespace', 'def', 'synonym', 'is_a', 'alt_id']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.terms = {}

    def read(self, file_handler, kept_fields = None, prefix_with = None):
        """
        Parses a OBO file v1.2 format.
        return all Term entries
        """
        self.terms = {}
        kept_fields = kept_fields is not None and kept_fields or self.DEFAULT_KEPT_FIELDS
    
        current_term = None
        for line in file_handler:
            line = line.strip()
            if not line:
                continue  # Skip empty

            if line == "[Term]":
                if current_term:
                    self.terms[current_term['id']] = current_term
                    for alt_id in current_term.get('alt_id', []):
                        self.terms[alt_id] = copy.deepcopy(current_term)
                        self.terms[alt_id]["id"] = alt_id
                current_term = {}

            elif line == "[Typedef]":  # Skip [Typedef] sections
                current_term = None

            elif current_term is not None:  # Only process if we're inside a [Term] environment
                key, value = line.split(": ", 1)
                
                if key not in kept_fields:
                    continue

                value = value.strip()

                m = re.search('"(.*)"', value)  # only keep text wrapped in double quotes (if any)
                if m:
                    value = m.group(1)

                if key == 'is_a':
                    value = value.split(' ! ')[0]

                if key not in current_term:
                    current_term[key] = value
                    if key in ['is_a', 'synonym']:  # force storing of these keys as list
                        current_term[key] = [value]
                else:
                    if not isinstance(current_term[key], list):
                        current_term[key] = [current_term[key]]
                    current_term[key].append(value)

        # Add last term
        if current_term is not None:
            self.terms[current_term['id']] = current_term

        if prefix_with:
            for term_id, term in self.terms.items():
                for key in term.keys():
                    if not key.startswith(prefix_with):
                        self.terms[term_id][prefix_with + key] = term.pop(key, None)

        return self.terms

    def get_parents(self, term, distance_from_child = 1):
        parents = []

        for term_id in term.get('is_a', []):
            parent = copy.deepcopy(self.terms[term_id])
            parent['distance'] = distance_from_child
            parents.append(parent)
            parents.extend(self.get_parents(parent, distance_from_child + 1))

        return parents
