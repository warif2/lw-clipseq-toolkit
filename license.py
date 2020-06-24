import requests
import datetime
import sys
import os
import pandas as pd
import csv
from pybedtools import BedTool

global chk_code

build_code = 'A0001'


def check_status():
    url = "http://warif2.web.illinois.edu/script_license/licenses.txt"
    global chk_code

    lfile = requests.get(url, allow_redirects=True)
    for line in lfile.iter_lines():
        if line.decode('utf-8')[0] == '#':
            continue
        bcode, ltype, param = line.decode('utf-8').split('\t')
        if bcode == build_code:
            if ltype == 'gen':
                if param == 'on':
                    chk_code = [200, 'active']
                else:
                    chk_code = [300, 'suspended']
            elif ltype == 'date':
                exp = datetime.datetime.strptime(param, '%m/%d/%Y') - datetime.datetime.today()
                if exp.days > 0:
                    chk_code = [200, 'active, expires in %i days' % exp.days]
                else:
                    chk_code = [300, 'expired']
            else:
                chk_code = [400, 'license error, please notify administrator.']
            break
        else:
            chk_code = [400, 'license does not exist, please notify administrator.']


def status_confirm():
    try:
        if chk_code[0] != 200:
            print('Error; license is expired.')
            sys.exit()
    except NameError:
        print('Error; code-100: encountered issues with script, please notify admin.')
        sys.exit()
    except IndexError:
        print('Error, code-110: encountered issues with script, please notify admin.')
        sys.exit()


def message():
    try:
        return chk_code[1]
    except NameError:
        print('Error; code-100: encountered issues with script, please notify admin.')
        sys.exit()
    except IndexError:
        print('Error, code-110: encountered issues with script, please notify admin.')
        sys.exit()


def peak_label(peak_file):
    # Run check on license
    status_confirm()

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
    # Run check on license
    status_confirm()

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


def gff2bed(gff, bed_output):
    # Run check on license
    status_confirm()

    command = 'gff2bed < {gff} > {bed_output}'
    status = os.system(command.format(gff=gff, bed_output=bed_output))
    if status == 0:
        print('Conversion of gff to bed was successful.')
    else:
        print('Error: Exit Code - %i' % status)


def get_introns(gff_bed):
    # Run check on license
    status_confirm()

    gff_bed = csv.reader(open(gff_bed, 'r'), delimiter='\t')
    gff_ds = {}

    for line in gff_bed:
        # Parse through information of gff entry and store relevant information into variables
        chrom = line[0]
        start = int(line[1])
        stop = int(line[2])
        strand = line[5]
        feature = line[7]
        attr_dict = dict()
        for attr in line[9].split(';'):
            name = attr.split('=')[0]
            value = attr.split('=')[1]
            attr_dict[name] = value

        # Add information into gff_ds
        if attr_dict['gene_id'] not in gff_ds.keys():
            if feature == 'gene':
                gff_ds[attr_dict['gene_id']] = {'gene': (start, stop), 'chr': chrom, 'strand': strand,
                                                'gene_name': attr_dict['gene_name'], 'exon': []}
            elif feature == 'exon':
                gff_ds[attr_dict['gene_id']] = {'chr': chrom, 'strand': strand, 'exon': [(start, stop)],
                                                'gene_name': attr_dict['gene_name']}
        else:
            if feature == 'gene':
                gff_ds[attr_dict['gene_id']]['gene'] = (start, stop)
            elif feature == 'exon':
                gff_ds[attr_dict['gene_id']]['exon'].append((start, stop))

    # Extract intron information from gff
    intron_bed = csv.writer(open('introns.bed', 'w'), delimiter='\t')
    for gene in gff_ds.keys():
        # Store gene range into variable intron_range as set
        intron_range = set(range(gff_ds[gene]['gene'][0], gff_ds[gene]['gene'][1] + 1))
        # Iterate through all exons within gene and remove those ranges from intron_range
        for exon in gff_ds[gene]['exon']:
            exon_range = set(range(exon[0], exon[1] + 1))  # Save exon range within exon_range as set
            intron_range = intron_range - exon_range  # Subtract exon_range from intron_range
        union = list(intron_range)  # Convert intron_range set to list and store it in union variable
        union.sort()  # Sort the coordinates in order

        # Chunk the coordinate lists into continous chunks
        results = []
        i = start = 0
        j = i + 1
        while j < len(union):
            if union[j] != union[i] + 1:
                results.append((union[start], union[j - 1] + 1))
                if j == len(union):
                    break
                i = start = j
                j = i + 1
            else:
                i, j = j, j + 1

        if start != j - 1:
            results.append((union[start], union[j - 1] + 1))

        # Save intron ranges into bed format
        for intron in results:
            info = gff_ds[gene]
            intron_bed.writerow([info['chr'], intron[0], intron[1], gene, '.', info['strand'], '.',
                                 'intron', '.', 'gene_name=' + info['gene_name']])

    print('Finished creating intron.bed annotation file.')


def merge_bed(bed_1, bed_2, out):
    # Run check on license
    status_confirm()

    command = 'cat {file1} {file2} > {output_file}'
    status = os.system(command.format(file1=bed_1, file2=bed_2, output_file=out))
    if status == 0:
        print('Conversion of gff3 to BED format was successful.')
    else:
        print('Error: Exit Code - %i' % status)


def peak_summary(apeak):
    # Run check on license
    status_confirm()

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


def list_avg(nlist):
    nsum = 0
    for elm in nlist:
        nsum += int(elm)
    avg = nsum / len(nlist)
    return avg


def rmats_filter(infile, pvalue_cut, fdr_cut, diff_cut, count_cut):
    status_confirm()
    summary_out = csv.writer(open("all_events_filtered_Pval-lt-%s_FDR-lt-%s_Diff-gt-%s_Count-gt-%s.csv" % (pvalue_cut, fdr_cut, diff_cut, count_cut), "a"), delimiter=",")

    with open(infile, 'r') as f:
        header = f.readline().strip("\n").split("\t")
        if "FDR" and "IncLevelDifference" not in header:
            print("Missing FDR or IncLevelDiffence in %s" % infile)
            return False

        if "PValue" not in header:
            print("Missing Pvalue in %s" % infile)
            return False

        if "IC_SAMPLE_1" not in header:
            if "IJC_SAMPLE_1" not in header:
                print("Missing either IC_SAMPLE_1 or IJC_SAMPLE_2 in %s" % infile)
                return False

        if "IC_SAMPLE_2" not in header:
            if "IJC_SAMPLE_2" not in header:
                print("Missing either IC_SAMPLE_2 or IJC_SAMPLE_2 in %s" % infile)
                return False

        if "SC_SAMPLE_1" not in header:
            if "SJC_SAMPLE_1" not in header:
                print("Missing either SC_SAMPLE_1 or SJC_SAMPLE_1 in %s" % infile)
                return False

        if "SC_SAMPLE_2" not in header:
            if "SJC_SAMPLE_2" not in header:
                print("Missing either SC_SAMPLE_2 or SJC_SAMPLE_2 in %s" % infile)
                return False

    outfile = infile.split("/")[-1].strip('.txt')
    event = outfile.split(".")[0]
    extype = [event]
    output = csv.writer(open("%s_filtered_PVal-lt-%s_FDR-lt-%s_Diff-gt-%s_Count-gt-%s.csv" % (outfile, pvalue_cut, fdr_cut, diff_cut, count_cut), "w"), delimiter="\t")

    if pvalue_cut == "none":
        pvalue_cut = 2
    if fdr_cut == "none":
        fdr_cut = 2
    if diff_cut == "none":
        diff_cut = 0

    with open(infile, 'r') as f:
        header = f.readline().strip("\n").split("\t")

        if event == "SE" or event == "A3SS" or event == "A5SS" or event == "RI":
            output.writerow(header[0:2] + ["Event"] + header[2:11] + [""] + [""] + header[11:23])
        elif event == "MXE":
            output.writerow(header[0:2] + ["Event"] + header[2:25])
        else:
            output.writerow(header)

        if event == "SE":
            summary_out.writerow(header[0:2] + ["Event"] + header[2:11] + [""] + [""] + header[11:23])

        pvalue_col = header.index("PValue")
        fdr_col = header.index("FDR")
        diff_col = header.index("IncLevelDifference")
        if "IC_SAMPLE_1" in header:
            ic_1_col = header.index("IC_SAMPLE_1")
            sc_1_col = header.index("SC_SAMPLE_1")
            ic_2_col = header.index("IC_SAMPLE_2")
            sc_2_col = header.index("SC_SAMPLE_2")

        if "IJC_SAMPLE_1" in header:
            ic_1_col = header.index("IJC_SAMPLE_1")
            sc_1_col = header.index("SJC_SAMPLE_1")
            ic_2_col = header.index("IJC_SAMPLE_2")
            sc_2_col = header.index("SJC_SAMPLE_2")

        for line in csv.reader(f, delimiter='\t'):

            if float(line[pvalue_col]) < float(pvalue_cut) and float(line[fdr_col]) < float(fdr_cut) and abs(
                    float(line[diff_col])) > float(diff_cut):
                copy = 0
                if count_cut == "none":
                    if event == "SE" or event == "RI" or event == "A5SS" or event == "A3SS":
                        output.writerow(line[0:2] + extype + line[2:11] + [""] + [""] + line[11:23])
                        summary_out.writerow(line[0:2] + extype + line[2:11] + [""] + [""] + line[11:23])
                    elif event == "MXE":
                        output.writerow(line[0:2] + extype + line[2:23])
                        summary_out.writerow(line[0:2] + extype + line[2:25])
                    else:
                        output.writerow(line)
                        summary_out.writerow(line)
                else:
                    inc_count_1 = list_avg(line[ic_1_col].split(","))
                    skp_count_1 = list_avg(line[sc_1_col].split(","))
                    inc_count_2 = list_avg(line[ic_2_col].split(","))
                    skp_count_2 = list_avg(line[sc_2_col].split(","))
                    if inc_count_1 < float(count_cut) and skp_count_1 < float(count_cut) and inc_count_2 < float(
                            count_cut) and skp_count_2 < float(count_cut):
                        copy = 1

                    if copy == 0:
                        if event == "SE" or event == "RI" or event == "A5SS" or event == "A3SS":
                            output.writerow(line[0:2] + extype + line[2:11] + [""] + [""] + line[11:23])
                            summary_out.writerow(line[0:2] + extype + line[2:11] + [""] + [""] + line[11:23])
                        elif event == "MXE":
                            output.writerow(line[0:2] + extype + line[2:23])
                            summary_out.writerow(line[0:2] + extype + line[2:25])
                        else:
                            output.writerow(line)
                            summary_out.writerow(line)
    print("Finished with %s" % infile)
    return True


def rmats_eclip_intersect(rmats_file, peak_file):
    # Run check on license
    status_confirm()

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
    check_status()
    status_confirm()
    print(chk_code)
