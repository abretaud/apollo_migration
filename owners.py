#!/usr/bin/env python

# A script to migrate ownership from an apollo1 gff into a migrated Apollo 2 db
# The script generates sql update queries ready to run on the postgresql db, like this:
#     python owners.py | grep update > owners.sql && psql -h db -Upostgres < owners.sql
# The owners must already exist in the target db (see README)

import re
import sys

from BCBio import GFF
from apollo import ApolloInstance

mrnas_attrs = {}
other_attrs = {}
scaffolds = set()
for rec in GFF.parse('apollo1_dump.gff'):
    for f in rec.features:
        if f.type == 'gene' or f.type == 'pseudogene':
            for sf in f.sub_features:
                sftype = sf.type
                if sftype == 'transcript':
                    sftype = 'mRNA'
                rna_id = "%s %s %s %s %s %s" % (sf.qualifiers['Name'][0], sftype, rec.id, sf.location.strand, sf.location.start, sf.location.end)
                if rna_id in mrnas_attrs:
                    print("Found multiple mrna with same id '%s' !! Aborting" % rna_id)
                    sys.exit()
                mrnas_attrs[rna_id] = f.qualifiers['owner'][0] + "@xxxxxxxxxx"
        else:
            feat_strand = f.location.strand
            if feat_strand is None:
                feat_strand = 0
            other_id = "%s %s %s %s %s %s" % (f.qualifiers['Name'][0], f.type, rec.id, feat_strand, f.location.start, f.location.end)
            if other_id in other_attrs:
                print("Found multiple other feat with same id '%s'!! Aborting" % other_id)
                sys.exit()
            if 'owner' in f.qualifiers:
                f.qualifiers['owner']
            other_attrs[other_id] = f.qualifiers['owner'][0] + "@xxxxxxxxxx"
        scaffolds.add(rec.id)

#print(mrnas_attrs)
#print(other_attrs)

org_common_name = "genus_species"

wa = ApolloInstance('http://localhost:8080', 'admin@xxxxxxxxxx', 'password...')


def print_sql(feat_uid, owner):
    print("update feature_grails_user set user_id = (select id from grails_user where last_name = '%s') where feature_owners_id = (select id from feature where unique_name = '%s');" % (owner, feat_uid))


for sequence in scaffolds:
    print("Processing %s %s" % (org_common_name, sequence))
    # Call setSequence to tell apollo which organism we're working with
    wa.annotations.set_sequence(org_common_name, sequence)
    # Then get a list of features.
    features = wa.annotations.get_features()
    # For each feature in the features
    for feat in features['features']:
        name = feat['name']
        if re.search('-[0-9]{5}$', name):
            name = name[:-6]

        current_id = '%s %s %s %s %s %s' % (name, feat['type']['name'], feat['sequence'], feat['location']['strand'], feat['location']['fmin'], feat['location']['fmax'])
        alt_current_id = '%s %s %s %s %s %s' % (name[:-1], feat['type']['name'], feat['sequence'], feat['location']['strand'], feat['location']['fmin'], feat['location']['fmax'])
        feat_uid = feat['uniquename']
        print("Treating %s..." % (current_id))

        if current_id in mrnas_attrs:
            print("    as an mRNA, %s owned by %s" % (feat_uid, mrnas_attrs[current_id]))
            print_sql(feat_uid, mrnas_attrs[current_id])
            if 'parent_id' in feat:
                print_sql(feat['parent_id'], mrnas_attrs[current_id])
        elif current_id in other_attrs:
            print("    as NOT an mRNA, %s owned by %s" % (feat_uid, other_attrs[current_id]))
            print_sql(feat_uid, other_attrs[current_id])
        elif alt_current_id in mrnas_attrs:
            current_id = alt_current_id
            print("    as an mRNA, alternate, %s owned by %s" % (feat_uid, mrnas_attrs[current_id]))
            print_sql(feat_uid, mrnas_attrs[current_id])
            if 'parent_id' in feat:
                print_sql(feat['parent_id'], mrnas_attrs[current_id])
        elif alt_current_id in other_attrs:
            current_id = alt_current_id
            print("    as NOT an mRNA, alternate, %s owned by %s" % (feat_uid, other_attrs[current_id]))
            print_sql(feat_uid, other_attrs[current_id])
        elif feat['type']['name'] == 'transcript':
            print("    WARNING: Skipping this feature, it's a transcript")
        else:
            print("    WARNING: Skipping this feature, unknown")
