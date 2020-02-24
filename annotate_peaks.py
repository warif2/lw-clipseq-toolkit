import argparse
import csv
from pybedtools import BedTool
import pandas as pd

def peak_label(peak_file):
    peak_raw = csv.reader(open(peak_file, 'r'), delimiter='\t')
    peak_filtered = csv.writer(open('peaks.bed', 'w'), delimiter='\t')
    peak_num = 1
    for line in peak_raw:
        out = line[0:3] + ['peak_' + str(peak_num)] + ['.'] + [line[5]]
        peak_num += 1
        peak_filtered.writerow(out)


def annotate_peaks(peak_bed, annotation_bed):
    """
    :param annotation_bed: BedTool object of annotation bed file.
    :param peak_bed: BedTool object of peaks bed file.
    :return:
    """
    # Initialize array for output
    output = pd.DataFrame(columns=['peak_id', 'gene_name', 'gene_id', 'gene_type', 'chr', 'start', 'stop',
                                   'strand', 'feature', 'exon_type'])

    # Get intersection of peak bed file with annotation bed file
    # intersect = peak_bed.intersect(annotation_bed, wa=True, wb=True, s=True)
    intersect = peak_bed.intersect(annotation_bed, wa=True, wb=True)

    # Generate a dictionary of annotation result for easy merging
    intersect_dict = {}
    for entry in intersect:

        peak = entry[3]
        if peak not in intersect_dict.keys():
            intersect_dict[peak] = {'feature': [], 'gene_name': [], 'gene_id': [], 'gene_type': [], "exon_type": []}

        feature = entry[13]

        # Fill in gene attributes
        if feature == 'gene':
            attr = entry[15].split(';')
            for ids in attr:
                id_name = ids.split('=')
                if id_name[0] == 'gene_name':
                    intersect_dict[peak]['gene_name'].append(id_name[1])
                if id_name[0] == 'gene_id':
                    intersect_dict[peak]['gene_id'].append(id_name[1])
                if id_name[0] == 'gene_type':
                    intersect_dict[peak]['gene_type'].append(id_name[1])

        # Store exon or intron features in feature
        elif feature == 'exon' or feature == 'intron':
            intersect_dict[peak]['feature'].append(feature)

        # Skip transcript features
        elif feature == 'transcript':
            continue

        # Store other feature types in exon_type
        else:
            if feature == 'UTR':
                utr_type = entry[9].split(':')[0]
                intersect_dict[peak]['exon_type'].append(utr_type)
            else:
                intersect_dict[peak]['exon_type'].append(feature)

    for entry in peak_bed:
        # Fill information into row_dict for output
        row_dict = dict()
        row_dict['peak_id'] = entry[3]
        row_dict['chr'] = entry[0]
        row_dict['start'] = entry[1]
        row_dict['stop'] = entry[2]
        row_dict['strand'] = entry[5]

        # Add annotation from intersect_dict into row_dict
        if row_dict['peak_id'] in intersect_dict.keys():
            row_dict['gene_name'] = ";".join(list(set(intersect_dict[row_dict['peak_id']]['gene_name'])))
            row_dict['gene_id'] = ";".join(list(set(intersect_dict[row_dict['peak_id']]['gene_id'])))
            row_dict['gene_type'] = ";".join(list(set(intersect_dict[row_dict['peak_id']]['gene_type'])))
            feat_list = list(set(intersect_dict[row_dict['peak_id']]['feature']))
            feat_list.sort()
            row_dict['feature'] = "-".join(feat_list)
            row_dict['exon_type'] = ";".join(list(set(intersect_dict[row_dict['peak_id']]['exon_type'])))
        else:
            row_dict['feature'] = 'intergenic'

        output = output.append(row_dict, ignore_index=True)

    return output


if __name__ == '__main__':
    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Annotate peak file with features.",
                                     prog="annotate_peaks.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-p", type=str, default=None, metavar="peak_file", help="specify path to peak file.",
                          required=True)
    required.add_argument("-gbed", type=str, default=None, metavar="gff_bed", help="specify path to gff3 bed file.",
                          required=True)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Label peaks
    peak_label(args.p)

    # Perform annotation intersection
    peak_bt = BedTool('peaks.bed')
    annotation_bt = BedTool(args.gbed)
    output_df = annotate_peaks(peak_bt, annotation_bt)
    output_df.to_csv('annotated_peaks.csv', index=False)
    print('Finished annotating peaks.')