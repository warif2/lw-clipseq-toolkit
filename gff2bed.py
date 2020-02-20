import argparse
import os
import csv


def gff2bed(gff, bed_output):
    command = 'gff2bed < {gff} > {bed_output}'
    status = os.system(command.format(gff=gff, bed_output=bed_output))
    if status == 0:
        print('Conversion of gff to bed was successful.')
    else:
        print('Error: Exit Code - %i' % status)


def get_introns(gff_bed):
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
    command = 'cat {file1} {file2} > {output_file}'
    status = os.system(command.format(file1=bed_1, file2=bed_2, output_file=out))
    if status == 0:
        print('Conversion of gff3 to BED format was successful.')
    else:
        print('Error: Exit Code - %i' % status)


if __name__ == '__main__':
    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Convert gff3 file to BED format.",
                                     prog="gff2bed.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-gff", type=str, default=None, metavar="GFF", help="specify path to gff3 file.",
                          required=True)
    required.add_argument("-out", type=str, default=None, metavar="OUTPUT BED",
                          help="specify desired label for output.", required=True)
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Convert gff to bed
    gff2bed(args.gff, 'features.bed')

    # Create introns.bed
    get_introns('features.bed')

    # Merge beds
    merge_bed('features.bed', 'introns.bed', args.out)
