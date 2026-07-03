import sys
import csv

with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type_unambig.csv', 'r') as f:
    csv_reader = csv.DictReader(f)

    for line in csv_reader:
        DMS_id = line['DMS_id']
        target_seq = line['target_seq']
        dms_transform_unambig = line['dms_transform_unambig']

        if dms_transform_unambig != 'yes':
            sys.exit(f'Error: DMS_id {DMS_id} has dms_transform_unambig = {dms_transform_unambig}, expected "yes"')

        print('>' + DMS_id)
        print(target_seq)

