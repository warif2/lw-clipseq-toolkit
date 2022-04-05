import argparse
import license
import sys

if __name__ == '__main__':
    # Check license
    license.check_status()

    # Setup of argparse for script arguments
    class licenseAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print("License status: %s" % license.message())
            sys.exit()

    # Setup of argparse for script arguments
    parser = argparse.ArgumentParser(description="Convert gff3 file to BED format.",
                                     prog="gff2bed.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-gff", type=str, default=None, metavar="GFF", help="specify path to gff3 file.",
                          required=True)
    required.add_argument("-out", type=str, default=None, metavar="OUTPUT DIR",
                          help="specify desired output path.", required=True)
    optional.add_argument("-l", "--license", action=licenseAction, metavar="", nargs=0,
                          help='show license status and exit')
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Convert gff to bed
    license.gff2bed(args.gff, args.out + '/features.bed')

    # Create introns.bed
    license.get_introns(args.out + '/features.bed', args.out)

    # Merge beds
    license.merge_bed(args.out + '/features.bed', args.out + '/introns.bed', args.out + '/annotation.bed')
