"""
SIGA.py is a command-line tool to generate Semantically Interoperable Genome Annotations from
GFF files [1,2] according to the RDF specification [3].

References:
[1] Generic Feature Format specification, http://www.sequenceontology.org/
[2] DDBJ/ENA/GenBank Feature Table Definition, http://www.insdc.org/documents/feature-table
[3] Resource Description Framework, https://www.w3.org/TR/rdf11-concepts/

Usage:
  SIGA.py -h|--help
  SIGA.py -v|--version
  SIGA.py db [-ruV] [-d DB_FILE | -e DB_FILEXT] GFF_FILE...
  SIGA.py rdf [-V] [-o FORMAT] [-c CFG_FILE] DB_FILE...

Arguments:
  GFF_FILE...      Input file(s) in GFF version 2 or 3.
  DB_FILE...       Input database file(s) in SQLite.

Options:
  -h, --help
  -v, --version
  -V, --verbose    Show verbose output in debug mode.
  -c FILE          Set the path of config file [default: config.ini]
  -d DB_FILE       Create a database from GFF file(s).
  -e DB_FILEXT     Set the database file extension [default: .db].
  -r               Check the referential integrity of the database(s).
  -u               Generate unique IDs for duplicated features.
  -o FORMAT        Output RDF graph in one of the following formats:
                     turtle (.ttl) [default: turtle]
                     nt (.nt),
                     n3 (.n3),
                     xml (.rdf)
"""

#
# Supported feature types:
#   genome, chromosome, gene, prim_transcript, mRNA, intron, exon, CDS,
#   three_prime_UTR, five_prime_UTR, polyA_site, polyA_sequence, variation
#
# URI data space defined relative to base URI:
#   ../genome/<species name>/<feature type>/<feature ID> + [#<begin|end>|#<start>-<end> chromosome only]
#


from __future__ import print_function
from docopt import docopt
from rdflib import Graph, URIRef, Literal
from rdflib.namespace import Namespace, RDF, RDFS, XSD, DCTERMS
from urllib2 import urlparse, unquote
from datetime import datetime
from ConfigParser import SafeConfigParser, NoSectionError, NoOptionError

import os
import sys
import re
import gffutils as gff
import sqlite3 as sql


__author__ = 'Arnold Kuzniar'
__version__ = '0.4.5'
__status__ = 'alpha'
__license__ = 'Apache License, Version 2.0'


def is_uri(uri):
    """Check if the input string is a URI.
    """
    u = urlparse.urlparse(uri)
    if u.scheme not in ('http', 'https', 'ftp'):
        raise ValueError("Unsupported URI scheme in '{0}'.".format(uri))
    if u.netloc is '':
        raise ValueError("No host specified in '{0}'.".format(uri))
    return u.geturl()


def validate(self):
    """Validate an instance of `SafeConfigParser` class. Check the presence of mandatory sections/options.
    """
    try:
        sections = ('GFF', 'RDF', 'FeatureRewrite', 'FeatureToClass',
                    'DNAstrandToClass', 'Ontologies')
        options = ('download_url', 'base_uri', 'creator', 'license')
        for s in sections:
            for k, v in self.items(s):
                try:
                    if (s in sections[:2] and k in options) or s in sections[3:]:
                        is_uri(v)
                except ValueError, err:
                    print(str(err), file=sys.stderr)
                    sys.exit(1)
    except NoSectionError, err:
        print(str(err), file=sys.stderr)
        sys.exit(1)
    try:
        self.get('GFF', 'species_name')
        self.get('GFF', 'description')
        self.get('GFF', 'download_url')
        try:
            self.getint('GFF', 'ncbitaxon_id')
        except ValueError, err:
            print('NCBI Taxon ID: {0}'.format(err), file=sys.stderr)
            sys.exit(1)
        self.get('RDF', 'version')
        self.get('RDF', 'base_uri')
        self.get('RDF', 'creator')
        self.get('RDF', 'license')
        self.get('FeatureToClass', 'genome')
        self.get('FeatureToClass', 'chromosome')
    except NoOptionError, err:
        print(str(err), file=sys.stderr)
        sys.exit(1)


def normalize_filext(s):
    """Prefix file extension with '.' if needed.
    """
    dot = '.'
    if s.startswith(dot) is False:
        s = dot + s
    return s


def remove_file(fn):
    """Remove file if exists.
    """
    if os.path.exists(fn) is True:
        os.remove(fn)


def normalize_feature_id(id):
    """Ad-hoc function to normalize feature IDs.
    """
    # Note: Sub-optimal URL resolution of genetic features in SGN (https://solgenomics.net/);
    # for example, feature ID 'gene:Solyc00g005000.2' in GFF corresponds to three URLs:
    #   https://solgenomics.net/locus/Solyc00g005000.2/view
    #   https://solgenomics.net/locus/Solyc00g005000/view
    #   https://solgenomics.net/feature/17660839/details
    #
    # Also note the use of 'locus' istead of 'gene' and iternal IDs.
    #
    # Features such as mRNA, exon, intron do not resolve same way:
    #   mRNA:Solyc00g005000.2.1       https://solgenomics.net/feature/17660840/details
    #   exon:Solyc00g005000.2.1.1     https://solgenomics.net/feature/17660841/details
    #   exon:Solyc00g005000.2.1.2     https://solgenomics.net/feature/17660843/details
    #   intron:Solyc00g005000.2.1.1   https://solgenomics.net/feature/17660842/details
    #   five_prime_UTR and three_prime_UTR do not to have corresponding URLs.
    #
    # Moreover, feature IDs are not opaque as these prefixed with type e.g.
    #   gene:Solyc00g005000.2.
    #
    # Normalizing feature IDs by removing the prefixes is reasonable for most
    # feature types except for UTRs because of resulting ambiguity;
    # for example, Solyc00g005000.2.1.0 for both
    #   five_prime_UTR:Solyc00g005000.2.1.0
    #   three_prime_UTR:Solyc00g005000.2.1.0
    return re.sub('gene:|mRNA:|CDS:|exon:|intron:|\w+UTR:', '', id, flags=re.IGNORECASE)


def get_feature_attrs(ft):
    """Concatenate feature attributes (9th column in GFF) into a single string.
    """
    attrs = {}  # selected attributes
    des = []

    for attr in ['Name', 'Note', 'Alias', 'Ontology_term', 'Interpro2go_term', 'Sifter_term']:
        attrs[attr] = None
        attrs[attr.lower()] = None

    for key in ft.attributes.keys():
        if key in attrs:
            val = ft[key]
            val = ', '.join(val) if type(val) == list else str(val)
            des.append('{0}: {1}'.format(key, val.encode('utf-8')))

    if len(des) == 0:
        return None
    else:
        return unquote('; '.join(sorted(des)))


def triplify(self, rdf_format, cfg):
    """Generate RDF triples from `FeatureDB` using Direct Mapping approach.
    """
    base_name, sfx = os.path.splitext(self.dbfn)
    graph = Graph()

    # setup namespace prefixes
    DCMITYPE = Namespace('http://purl.org/dc/dcmitype/')
    graph.bind('dcmitype', DCMITYPE)
    graph.bind('dcterms', DCTERMS)

    for (prefix, uri) in cfg.items('Ontologies'):
        # instantiate dynamically from config
        exec(prefix + " = Namespace('{0}')".format(uri))
        graph.bind(prefix.lower(), eval(prefix))

    # lookup table for RDF mime-types and file extensions
    format_to_filext = dict(
        turtle=['text/turtle', '.ttl'],
        nt=['application/n-triples', '.nt'],
        n3=['text/n3', '.n3'],
        xml=['application/rdf+xml', '.rdf'])

    if rdf_format not in format_to_filext:
        raise IOError(
            "Unsupported RDF serialization '{0}'.".format(rdf_format))

    rdf_mime_type = format_to_filext[rdf_format][0]
    base_uri = URIRef(cfg.get('RDF', 'base_uri'))
    creator_uri = URIRef(cfg.get('RDF', 'creator'))
    license_uri = URIRef(cfg.get('RDF', 'license'))
    version = cfg.get('RDF', 'version')
    download_url = URIRef(cfg.get('GFF', 'download_url'))
    species_name = cfg.get('GFF', 'species_name')
    description = cfg.get('GFF', 'description')
    taxon_id = cfg.getint('GFF', 'ncbitaxon_id')
    taxon_uri = OBO.term('NCBITaxon_{0}'.format(taxon_id))
    genome_uri = URIRef(os.path.join(
        base_uri, 'genome', species_name.replace(' ', '_')))

    genome_type_uri = URIRef(cfg.get('FeatureToClass', 'genome'))

    # add genome info to graph
    graph.add((genome_uri, RDF.type, genome_type_uri))
    graph.add((genome_uri, RDF.type, DCMITYPE.Dataset))
    graph.add((genome_uri, RDFS.label, Literal(
        'genome of {0}'.format(species_name), lang="en")))
    graph.add((genome_uri, RDFS.comment, Literal(description, lang="en")))
    graph.add((genome_uri, DCTERMS.title, Literal(
        'genome of {0}'.format(species_name), lang="en")))
    graph.add((genome_uri, DCTERMS.description, Literal(
        description, lang="en")))
    graph.add((genome_uri, DCTERMS.source, URIRef(download_url)))
    graph.add((genome_uri, DCTERMS.creator, URIRef(creator_uri)))
    graph.add((genome_uri, DCTERMS.created, Literal(
        datetime.now().strftime("%Y-%m-%d"), datatype=XSD.date)))
    graph.add((genome_uri, DCTERMS.license, URIRef(license_uri)))
    graph.add((genome_uri, DCTERMS.term('format'), Literal(
        rdf_mime_type, datatype=XSD.string)))
    graph.add((genome_uri, DCTERMS.hasVersion, Literal(
        version, datatype=XSD.string)))
    # SO predicate no domain/range defined
    graph.add((genome_uri, SO.genome_of, taxon_uri))
    # RO predicate 'in taxon'
    graph.add((genome_uri, OBO.RO_0002162, taxon_uri))
    graph.add((taxon_uri, RDFS.label, Literal(
        'NCBI Taxonomy ID: {0}'.format(taxon_id), datatype=XSD.string)))
    graph.add((taxon_uri, DCTERMS.identifier, Literal(
        taxon_id, datatype=XSD.positiveInteger)))

    # loop through all features in GFF
    for feature in self.all_features():
        try:
            chromosome = str(feature.seqid)
            chromosome_uri = URIRef(os.path.join(
                genome_uri, 'chromosome', chromosome))
            chromosome_type_uri = URIRef(
                cfg.get('FeatureToClass', 'chromosome'))
            feature_id = normalize_feature_id(feature.id)
            feature_type = feature.featuretype
            try:
                feature_type = cfg.get('FeatureRewrite', feature_type)
            except NoOptionError:
                pass
            feature_uri = URIRef(os.path.join(
                genome_uri, feature_type, feature_id))
            feature_type_uri = URIRef(cfg.get('FeatureToClass', feature_type))
            strand_type_uri = URIRef(
                cfg.get('DNAstrandToClass', feature.strand))
            region_uri = URIRef(
                '{0}#{1}-{2}'.format(chromosome_uri, feature.start, feature.end))
            start_uri = URIRef('{0}#{1}'.format(chromosome_uri, feature.start))
            end_uri = URIRef('{0}#{1}'.format(chromosome_uri, feature.end))

            # add chromosome info to graph
            # Note: 'seqid' column refers to chromosome here
            graph.add((chromosome_uri, RDF.type, chromosome_type_uri))
            graph.add((chromosome_uri, RDFS.label, Literal(
                'chromosome {0}'.format(chromosome), datatype=XSD.string)))
            graph.add((chromosome_uri, SO.part_of, genome_uri))
            graph.add((feature_uri, RDF.type, feature_type_uri))
            graph.add((feature_uri, RDFS.label, Literal(
                '{0} {1}'.format(feature_type, feature_id), datatype=XSD.string)))
            graph.add((feature_uri, DCTERMS.identifier, Literal(
                feature_id, datatype=XSD.string)))

            # add feature descriptions (from the attributes field) to graph
            des = get_feature_attrs(feature)
            if des is not None:
                graph.add((feature_uri, RDFS.comment, Literal(
                    des, datatype=XSD.string)))

            # add feature start/end coordinates and strand info to graph
            graph.add((feature_uri, FALDO.location, region_uri))
            graph.add((region_uri, RDF.type, FALDO.Region))
            graph.add((region_uri, RDFS.label, Literal(
                'chromosome {0}:{1}-{2}'.format(chromosome, feature.start, feature.end))))
            graph.add((region_uri, FALDO.begin, start_uri))
            graph.add((start_uri, RDF.type, FALDO.ExactPosition))
            graph.add((start_uri, RDF.type, strand_type_uri))
            graph.add((start_uri, RDFS.label, Literal(
                'chromosome {0}:{1}-*'.format(chromosome, feature.start))))
            graph.add((start_uri, FALDO.position, Literal(
                feature.start, datatype=XSD.positiveInteger)))
            graph.add((start_uri, FALDO.reference, chromosome_uri))
            graph.add((region_uri, FALDO.end, end_uri))
            graph.add((end_uri, RDF.type, FALDO.ExactPosition))
            graph.add((end_uri, RDF.type, strand_type_uri))
            graph.add((end_uri, RDFS.label, Literal(
                'chromosome {0}:*-{1}'.format(chromosome, feature.end))))
            graph.add((end_uri, FALDO.position, Literal(
                feature.end, datatype=XSD.positiveInteger)))
            graph.add((end_uri, FALDO.reference, chromosome_uri))
            # Note: phase is mandatory for CDS but has no corresponding term in
            # ontology

            # add parent-child relationships between features to graph
            for child in self.children(feature, level=1):
                child_feature_id = normalize_feature_id(child.id)
                child_feature_type = child.featuretype
                try:
                    child_feature_type = cfg.get(
                        'FeatureRewrite', child_feature_type)
                except NoOptionError:
                    pass
                child_feature_uri = URIRef(os.path.join(
                    genome_uri, child_feature_type, child_feature_id))
                child_feature_type_uri = URIRef(
                    cfg.get('FeatureToClass', child_feature_type))
                graph.add((feature_uri, SO.has_part, child_feature_uri))
                if feature_type == 'gene' and child_feature_type == 'prim_transcript':
                    graph.add((feature_uri, SO.transcribed_to, child_feature_uri))
        except:
            pass
    rdf_file = base_name + format_to_filext[rdf_format][1]
    with open(rdf_file, 'w') as fout:
        fout.write(graph.serialize(format=rdf_format))

if __name__ == '__main__':
    args = docopt(__doc__, version=__version__)
    debug = args['--verbose']
    sys.tracebacklimit = debug

    if args['db'] is True:  # db mode
        unique_keys = 'create_unique' if args['-u'] is True else 'error'
        fk_check = 'ON' if args['-r'] is True else 'OFF'
        pragmas = dict(foreign_keys=fk_check)

        # populate database(s) from GFF file(s)
        for gff_file in args['GFF_FILE']:
            if os.path.exists(gff_file) is False:
                raise IOError("GFF file '{0}' not found.".format(gff_file))

            if args['-d']:  # one db for all GFF files
                db_file = args['-d']
                if os.path.exists(db_file) is True:
                    db = gff.FeatureDB(db_file)
                    db.update(gff_file)
                else:
                    db = gff.create_db(gff_file, db_file, merge_strategy=unique_keys,
                                       verbose=debug, pragmas=pragmas, force=False)
            else:  # one db per GFF file
                base_name = os.path.splitext(gff_file)[0]
                db_file = base_name + normalize_filext(args['-e'])
                try:
                    db = gff.create_db(gff_file, db_file, merge_strategy=unique_keys,
                                       verbose=debug, pragmas=pragmas, force=False)
                except sql.OperationalError:
                    raise IOError(
                        "Database file '{0}' already exists.".format(db_file))
                except sql.IntegrityError, err:
                    remove_file(db_file)
                    raise IOError(
                        "{0} in database '{1}'.".format(err, db_file))
    else:  # rdf mode
        rdf_format = args['-o']
        cfg_file = args['-c']
        if os.path.exists(cfg_file) is False:
            raise IOError("Config file '{0}' not found.".format(cfg_file))
        SafeConfigParser.validate = validate  # add new method
        cfg = SafeConfigParser()
        cfg.optionxform = str  # option names case sensitive
        cfg.read(cfg_file)
        cfg.validate()
        gff.FeatureDB.triplify = triplify # add new method
        for db_file in args['DB_FILE']:
            db = gff.FeatureDB(db_file)
            db.triplify(rdf_format, cfg)
