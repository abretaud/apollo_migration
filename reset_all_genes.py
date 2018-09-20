#!/usr/bin/env python

from BCBio import GFF
from apollo import ApolloInstance

# A scary script to remove all annotated genes from an apollo instance
# Will only remove features on scaffolds present in gff file
# (was used after a failed migration)

scaffolds = set()
for rec in GFF.parse('apollo1_dump.gff'):
    for f in rec.features:
        scaffolds.add(rec.id)

wa = ApolloInstance('http://localhost:8080', 'admin@xxxxxxx', 'password...')

org_common_name = "genus_species"

for sequence in scaffolds:
    print("Processing %s %s" % (org_common_name, sequence))
    # Call setSequence to tell apollo which organism we're working with
    wa.annotations.set_sequence(org_common_name, sequence)
    # Then get a list of features.
    features = wa.annotations.get_features()
    # For each feature in the features
    for feature in features['features']:
        wa.annotations.delete_feature(feature['uniquename'])
        print('Deleted %s [type=%s]' % (feature['uniquename'], feature['type']['name']))
