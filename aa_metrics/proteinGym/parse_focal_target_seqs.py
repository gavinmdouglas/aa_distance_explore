import csv

with open('/Users/gavin/Drive/research/aa_distance/aa_distance_zenodo/proteinGym/focal_dms_score_type.csv', 'r') as f:
    csv_reader = csv.DictReader(f)

    for line in csv_reader:
        DMS_id = line['DMS_id']
        target_seq = line['target_seq']
        print('>' + DMS_id)
        print(target_seq)

