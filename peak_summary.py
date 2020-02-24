import pandas as pd
import csv
import argparse


def peak_summary(apeak):
    # Setup dictionary for counting
    annotation_count = {'intergenic': 0, 'exon': 0, 'intron': 0, 'exon-intron': 0,
                        'CDS': 0, 'five_prime_UTR': 0, 'three_prime_UTR': 0}
    f = csv.reader(open(apeak, 'r'))

    # Iterate through file and count
    for ln in f:
        print(ln)
        # Skip header
        if ln[0] == 'peak_id':
            continue
        if ln[10] == '':
            continue

        # Count feature type
        annotation_count[ln[10]] += 1

        # Count exon_type feature
        ex_feature = ln[11].split(';')
        for feat in ex_feature:
            if feat in annotation_count.keys():
                annotation_count[feat] += 1

    # Save summary to csv file
    feat_count_df = pd.DataFrame.from_dict([annotation_count])
    feat_count_df = feat_count_df[['intergenic', 'exon', 'intron', 'exon-intron', 'CDS', 'five_prime_UTR',
                                   'three_prime_UTR']]
    feat_count_df.replace(r'\s+', 0, regex=True, inplace=True)
    feat_count_df.to_csv('feature_count_summary.csv', index=False)


if __name__ == '__main__':
    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Summarizes feature counts from annotated peak file.",
                                     prog="peak_summary.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-ap", type=str, default=None, metavar="annotated_peak file",
                          help="specify path to annotated peaks file.", required=True)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    peak_summary(args.ap)
    print('Summary complete!')
