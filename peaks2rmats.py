# Import modules
import pandas as pd
from pybedtools import BedTool
import argparse


def rmats_eclip_intersect(rmats_file, peak_file):
    # Format peak file to allow to input as BedTool object
    ap = pd.read_csv(peak_file)
    ap_input = """"""
    for index, row in ap.iterrows():
        try:
            line = ' '.join([row['chr'], str(row['start']), str(row['stop']), row['peak_id'], '.', row['strand'],
                             row['feature'], str(row['exon_type']), str(row['gene_id'])])
            ap_input = ap_input + line + '\n'
        except:
            continue

    peaks = BedTool(ap_input, from_string=True)

    # Format rMATS output to allow input as BedTool object
    rm = pd.read_csv(rmats_file)
    rm_input = """"""
    for index, row in rm.iterrows():
        if row['Event'] == 'SE' or row['Event'] == 'RI':
            if row['strand'] == '+':
                upper_array = [row['chr'], str(row['upstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-upper', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['exonEnd']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                lower_array = [row['chr'], str(row['exonEnd']), str(row['downstreamEE']), str(row['ID']) + '-lower',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
            if row['strand'] == '-':
                lower_array = [row['chr'], str(row['upstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-lower', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['exonEnd']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                upper_array = [row['chr'], str(row['exonEnd']), str(row['downstreamEE']), str(row['ID']) + '-upper',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]

        if row['Event'] == 'MXE':
            if row['strand'] == '+':
                upper_array = [row['chr'], str(row['downstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-upper', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamEE']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                lower_array = [row['chr'], str(row['upstreamEE']), str(row['Unnamed: 13']), str(row['ID']) + '-lower',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
            if row['strand'] == '-':
                lower_array = [row['chr'], str(row['downstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-lower', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamEE']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                upper_array = [row['chr'], str(row['upstreamEE']), str(row['Unnamed: 13']), str(row['ID']) + '-upper',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]

        if row['Event'] == 'A5SS':
            if row['strand'] == '+':
                upper_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamEE']),
                               str(row['ID']) + '-upper', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['upstreamEE']), str(row['exonEnd']), str(row['ID']) + '-exon', '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                lower_array = [row['chr'], str(row['exonEnd']), str(row['downstreamEE']), str(row['ID']) + '-lower',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
            if row['strand'] == '-':
                upper_array = [row['chr'], str(row['upstreamES']), str(row['exonEnd']), str(row['ID']) + '-upper', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamES']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                lower_array = [row['chr'], str(row['downstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-lower', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]

        if row['Event'] == 'A3SS':
            if row['strand'] == '-':
                lower_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamEE']),
                               str(row['ID']) + '-lower', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['upstreamEE']), str(row['exonEnd']), str(row['ID']) + '-exon', '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                upper_array = [row['chr'], str(row['exonEnd']), str(row['downstreamEE']), str(row['ID']) + '-upper',
                               '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
            if row['strand'] == '+':
                lower_array = [row['chr'], str(row['upstreamES']), str(row['exonEnd']), str(row['ID']) + '-lower', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]
                exon_array = [row['chr'], str(row['exonStart_0base']), str(row['upstreamES']), str(row['ID']) + '-exon',
                              '.',
                              row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                              row['GeneID']]
                upper_array = [row['chr'], str(row['downstreamES']), str(row['exonStart_0base']),
                               str(row['ID']) + '-upper', '.',
                               row['strand'], str(row['IncLevelDifference']), str(row['FDR']), row['Event'],
                               row['GeneID']]

        rm_input = rm_input + ' '.join(exon_array) + '\n'
        rm_input = rm_input + ' '.join(upper_array) + '\n'
        rm_input = rm_input + ' '.join(lower_array) + '\n'

    rmats = BedTool(rm_input, from_string=True)

    # Perform intersection of rmats with peaks
    intersect = rmats.intersect(peaks, wa=True, wb=True, s=True)
    # intersect.saveas('output.bed')

    intersect_dict = {}
    for line in intersect:
        rid = line[3].split('-')[0]
        region = line[3].split('-')[1]
        event = line[8]
        gene_id = line[9]
        peak_id = line[13]
        region_type = line[16].split(';')
        feature_type = line[17].split(';')
        code = rid + '-' + event
        if code not in intersect_dict.keys():
            intersect_dict[code] = {'ID': rid, 'GeneID': gene_id, 'Event': event, 'exon_count': 0, 'lower_count': 0,
                                    'upper_count': 0,
                                    'total': 0, 'exon': [], 'upper': [], 'lower': [], 'Region_type': region_type,
                                    'Feature_type': feature_type}
            intersect_dict[code][region].append(peak_id)
            intersect_dict[code][region + '_count'] += 1
            intersect_dict[code]['total'] += 1
        else:
            for elm in region_type:
                if elm not in intersect_dict[code]['Region_type']:
                    intersect_dict[code]['Region_type'].append(elm)
            for elm in feature_type:
                if elm not in intersect_dict[code]['Feature_type']:
                    intersect_dict[code]['Feature_type'].append(elm)
            intersect_dict[code][region].append(peak_id)
            intersect_dict[code][region + '_count'] += 1
            intersect_dict[code]['total'] += 1

    intersect_list = []
    columns = ['exon', 'lower', 'upper', 'Region_type', 'Feature_type']
    for entry in intersect_dict.keys():
        for col in columns:
            intersect_dict[entry][col] = ';'.join(intersect_dict[entry][col])
            if intersect_dict[entry][col] == '':
                intersect_dict[entry][col] = 'NA'
        # intersect_dict[entry]['exon'] = ';'.join(intersect_dict[entry]['exon'])
        # intersect_dict[entry]['lower'] = ';'.join(intersect_dict[entry]['lower'])
        # intersect_dict[entry]['upper'] = ';'.join(intersect_dict[entry]['upper'])
        # intersect_dict[entry]['Region_type'] = ';'.join(intersect_dict[entry]['Region_type'])
        # intersect_dict[entry]['Feature_type'] = ';'.join(intersect_dict[entry]['Feature_type'])
        intersect_list.append(intersect_dict[entry])

    # Create pandas dataframe of the intersection_list
    intersect_df = pd.DataFrame.from_dict(intersect_list)
    intersect_df = intersect_df[['ID', 'GeneID', 'Event', 'Region_type', 'Feature_type', 'exon', 'lower',
                                 'upper', 'exon_count', 'lower_count', 'upper_count', 'total']]

    # Setup output_df from rm dataframe for merge with intersect_df
    rm.ID = rm.ID.astype(str)
    output_df = rm[['ID', 'GeneID', 'Event', 'geneSymbol', 'PValue', 'FDR', 'IncLevelDifference']]
    output_df = pd.merge(output_df, intersect_df, on=['ID', 'GeneID', 'Event'], how='left').fillna('NA')

    # Return table
    return output_df


if __name__ == '__main__':
    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Finds peaks on .",
                                         prog="peak_summary.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-r", type=str, default=None, metavar="rmats_file",
                          help="specify path to filtered rmats file.", required=True)
    required.add_argument("-ap", type=str, default=None, metavar="annotated_peak file",
                          help="specify path to annotated peaks file.", required=True)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Perform clip intersection
    intersection_df = rmats_eclip_intersect(args.r, args.ap)
    intersection_df.to_csv('rmats_with_clip.csv', index=False)

    print('Summary complete!')