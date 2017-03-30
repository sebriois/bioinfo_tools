from urllib import request, parse
import os
import time
import json
import subprocess

from bioinfo_tools.utils.sge import SgeJob


def fetch_taxonomy(species_name):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = parse.urlencode({
        'db': 'taxonomy',
        'term': species_name,
        'retmode': 'json'
    })
    uri = url + '?' + params
    
    with request.urlopen(uri) as r:
        ncbi_response = r.read()
    
    json_response = json.loads(ncbi_response.decode())
    try:
        ncbi_id = json_response.get("esearchresult", {}).get("idlist", [])[0]
    except IndexError:
        print("[%s] Can't figure out taxonomy ID using URL: %s" % (species_name, uri))
        return None, None
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + ncbi_id
    with request.urlopen(url) as r:
        ncbi_response = r.read()
    ncbi_xml = ncbi_response.decode()
    
    return ncbi_id, ncbi_xml


def makeblastdb(fasta_filepath):
    cmd = [
        "makeblastdb",
        "-in", fasta_filepath,
        "-dbtype", "prot"
    ]
    subprocess.check_call(" ".join(cmd), shell = True)


class Blastp(object):
    def __init__(self, use_cluster = False, **kwargs):
        for arg_name in ['query', 'db', 'out']:
            if arg_name not in kwargs:
                raise Exception("'%s' argument is required" % arg_name)

        self.args = kwargs
        if not os.path.exists(self.args['query']):
            raise Exception("query file not found: " + self.args['query'])

        if isinstance(self.args['db'], list):
            self.args['db'] = "'" + " ".join(self.args['db']) + "'"

        blast_cmd = ['blastp']

        for key, value in self.args.items():
            if not key.startswith('-'):
                key = '-' + key
            blast_cmd.append(key + '=' + str(value))
        blast_cmd = " ".join(blast_cmd)

        if use_cluster:
            filename = os.path.basename(self.args['out']).split('.')[0]
            blast_job = SgeJob()
            blast_job.submit(blast_cmd, 'BLASTP_%s' % filename)
        else:
            try:
                subprocess.check_output(blast_cmd, shell = True)
            except subprocess.CalledProcessError as e:
                raise Exception(e.output)

        # wait for output file to be writen on disk
        max_checks = 10
        while not os.path.exists(self.args['out']) and max_checks > 0:
            time.sleep(10)
            max_checks -= 1

        if not os.path.exists(self.args['out']):
            raise Exception("blast output does not exist: " + self.args['out'])


class Blastx(object):
    def __init__(self, use_cluster = False, **kwargs):
        for arg_name in ['query', 'db', 'out']:
            if arg_name not in kwargs:
                raise Exception("'%s' argument is required" % arg_name)

        self.args = kwargs
        if not os.path.exists(self.args['query']):
            raise Exception("query file not found: " + self.args['query'])

        if isinstance(self.args['db'], list):
            self.args['db'] = "'" + " ".join(self.args['db']) + "'"

        blast_cmd = ['blastx']

        for key, value in self.args.items():
            if not key.startswith('-'):
                key = '-' + key
            blast_cmd.append(key + '=' + str(value))
        blast_cmd = " ".join(blast_cmd)

        if use_cluster:
            filename = os.path.basename(self.args['out']).split('.')[0]
            blast_job = SgeJob()
            blast_job.submit(blast_cmd, 'BLASTX_%s' % filename)
        else:
            try:
                subprocess.check_output(blast_cmd, shell = True)
            except subprocess.CalledProcessError as e:
                raise Exception(e.output)

        # wait for output file to be writen on disk
        max_checks = 10
        while not os.path.exists(self.args['out']) and max_checks > 0:
            time.sleep(10)
            max_checks -= 1

        if not os.path.exists(self.args['out']):
            raise Exception("blast output does not exist: " + self.args['out'])
