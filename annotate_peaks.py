import sys
import argparse
import license
from pybedtools import BedTool

if __name__ == '__main__':
    # Check license
    license.check_status()

    # Setup of argparse for script arguments
    class licenseAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print("License status: %s" % license.message())
            sys.exit()

    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Annotate peak file with features.",
                                     prog="annotate_peaks.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-pk", type=str, default=None, metavar="peak_file", help="specify path to peak file.",
                          required=True)
    required.add_argument("-gbed", type=str, default=None, metavar="gff_bed", help="specify path to gff3 bed file.",
                          required=True)
    required.add_argument("-out", type=str, default=None, metavar="output_dir",help="specify desired output path.",
                          required=True)
    optional.add_argument("-l", "--license", action=licenseAction, metavar="", nargs=0,
                          help='show license status and exit')
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Label peaks
    license.peak_label(args.pk, args.out)

    # Perform annotation intersection
    peak_bt = BedTool(args.out + '/peaks.bed')
    annotation_bt = BedTool(args.gbed)
    output_df = license.annotate_peaks(peak_bt, annotation_bt)
    output_df.to_csv(args.out + '/annotated_peaks.csv', index=False)
    print('Finished annotating peaks.')
