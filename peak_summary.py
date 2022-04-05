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
    parser = argparse.ArgumentParser(description="Summarizes feature counts from annotated peak file.",
                                     prog="peak_summary.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-ap", type=str, default=None, metavar="annotated_peak file",
                          help="specify path to annotated peaks file.", required=True)
    required.add_argument("-out", type=str, default=None, metavar="output_dir", help="specify desired output path.",
                          required=True)
    optional.add_argument("-l", "--license", action=licenseAction, metavar="", nargs=0,
                          help='show license status and exit')
    parser._action_groups.append(optional)
    args = parser.parse_args()

    license.peak_summary(args.ap, args.out)
    print('Summary complete!')
