#!/usr/bin/env python

# Takes apollo attributes from Apollo1 gff dump and transfer them to the new Apollo2 server
# Requires latest code from python-apollo github repo
# Does not transfer everything
# there is some bipaa-specific magic, maybe not relevant for other users

import re
import sys

from BCBio import GFF
from apollo import ApolloInstance


mrnas_attrs = {}
other_attrs = {}
scaffolds = set()
for rec in GFF.parse('apollo1_dump.gff'):
    for f in rec.features:
        scaffolds.add(rec.id)
        if f.type == 'gene' or f.type == 'pseudogene':
            for sf in f.sub_features:
                sftype = sf.type
                if sftype == 'transcript':
                    sftype = 'mRNA'
                rna_id = "%s %s %s %s %s %s" % (sf.qualifiers['Name'][0], sftype, rec.id, sf.location.strand, sf.location.start, sf.location.end)
                if rna_id in mrnas_attrs:
                    print("Found multiple mrna with same id '%s' !! Aborting" % rna_id)
                    sys.exit()
                mrnas_attrs[rna_id] = {'gene': f.qualifiers, 'mrna': sf.qualifiers}
        else:
            feat_strand = f.location.strand
            if feat_strand is None:
                feat_strand = 0
            other_id = "%s %s %s %s %s %s" % (f.qualifiers['Name'][0], f.type, rec.id, feat_strand, f.location.start, f.location.end)
            if other_id in other_attrs:
                print("Found multiple other feat with same id '%s'!! Aborting" % other_id)
                sys.exit()
            other_attrs[other_id] = {'gene': f.qualifiers}

#print(mrnas_attrs)
#print(other_attrs)

org_common_name = "genus_species"

wa = ApolloInstance('http://localhost:8080', 'admin@xxxxxxxxxx', 'password...')


def apply_attrs(feature_id, attrs):
    if 'symbol' in attrs:
        wa.annotations.set_symbol(feature_id, attrs['symbol'][0])
    if 'description' in attrs:
        wa.annotations.set_description(feature_id, attrs['description'][0])
    if 'Name' in attrs:
        wa.annotations.set_name(feature_id, attrs['Name'][0])
    if 'Dbxref' in attrs:
        for dbx in attrs['Dbxref']:
            dbx = dbx.split(':')
            if len(dbx) == 2:
                wa.annotations.add_dbxref(feature_id, dbx[0], dbx[1])
    if 'Note' in attrs:
        wa.annotations.add_comment(feature_id, attrs['Note'])

    if 'AnnotGroup' in attrs:
        wa.annotations.add_attribute(feature_id, 'AnnotGroup', attrs['AnnotGroup'][0])


for scaf in scaffolds:
    print("Working on scaffold %s" % scaf)
    data = wa.annotations.set_sequence(org_common_name, scaf)

    data = wa.annotations.get_features()

    for feat in data['features']:
        name = feat['name']
        if re.search('-[0-9]{5}$', name):
            name = name[:-6]

        current_id = '%s %s %s %s %s %s' % (name, feat['type']['name'], feat['sequence'], feat['location']['strand'], feat['location']['fmin'], feat['location']['fmax'])
        alt_current_id = '%s %s %s %s %s %s' % (name[:-1], feat['type']['name'], feat['sequence'], feat['location']['strand'], feat['location']['fmin'], feat['location']['fmax'])
        feat_uid = feat['uniquename']
        print("Treating %s..." % (current_id))

        if current_id in mrnas_attrs:
            print("    as an mRNA")

            # if 'Rep1_Hd35' in current_id or 'U1_Hd19' in current_id:
            apply_attrs(feat_uid, mrnas_attrs[current_id]['mrna'])

            if 'parent_id' in feat:
                apply_attrs(feat['parent_id'], mrnas_attrs[current_id]['gene'])

        elif current_id in other_attrs:
            print("    as NOT an mRNA")
            apply_attrs(feat_uid, other_attrs[current_id]['gene'])
        elif alt_current_id in mrnas_attrs:
            print("    as an mRNA, alternate")
            current_id = alt_current_id
            apply_attrs(feat_uid, mrnas_attrs[current_id]['mrna'])

            if 'parent_id' in feat:
                apply_attrs(feat['parent_id'], mrnas_attrs[current_id]['gene'])
        elif alt_current_id in other_attrs:
            print("    as NOT an mRNA, alternate")
            current_id = alt_current_id
            apply_attrs(feat_uid, other_attrs[current_id]['gene'])
        elif feat['type']['name'] == 'transcript':
            print("    WARNING: Skipping this feature, it's a transcript")
        else:
            print("    WARNING: Skipping this feature, unknown")
