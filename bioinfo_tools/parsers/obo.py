import copy
import re
import ssl

import urllib

# FIXME: bad practice, unsafe, etc.
from bioinfo_tools.utils.log import Log

ssl._create_default_https_context = ssl._create_unverified_context


class OboParser(Log):
    DEFAULT_KEPT_FIELDS = ['id', 'name', 'namespace', 'def', 'synonym', 'is_a']

    def __init__(self, use_proxy = False, **kwargs):
        super().__init__(**kwargs)
        self.terms = {}
        self.use_proxy = use_proxy

    def read(self, uri, kept_fields = None, prefix_with = None):
        """
        Parses a OBO file v1.2 format.
        return all Term entries
        """
        self.terms = {}
        kept_fields = kept_fields is not None and kept_fields or self.DEFAULT_KEPT_FIELDS

        # set proxies as empty
        if not self.use_proxy:
            proxy_handler = urllib.request.ProxyHandler({})
            opener = urllib.request.build_opener(proxy_handler)
            urllib.request.install_opener(opener)

        self.log("Downloading %s" % uri)
        with urllib.request.urlopen(uri) as f:
            current_term = None
            for line in f:
                line = line.decode().strip()
                if not line:
                    continue  # Skip empty

                if line == "[Term]":
                    if current_term:
                        self.terms[current_term['id']] = current_term
                    current_term = {}

                elif line == "[Typedef]":  # Skip [Typedef] sections
                    current_term = None

                elif current_term is not None:  # Only process if we're inside a [Term] environment
                    key, value = line.split(": ", 1)
                    if key == 'is_obsolete':
                        current_term = None
                        continue

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
