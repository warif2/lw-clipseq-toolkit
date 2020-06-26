# Import modules
import sys
import argparse
import license


if __name__ == '__main__':
    # Check license
    license.check_status()

    # Setup of argparse for script arguments
    class licenseAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print("License status: %s" % license.message())
            sys.exit()

    parser = argparse.ArgumentParser(description="Finds peaks on exons found in rmats output.",
                                         prog="peaks2rmats.py")
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument("-rf", type=str, default=None, metavar="rmats_file",
                          help="specify path to filtered rmats file", required=True)
    required.add_argument("-ap", type=str, default=None, metavar="annotated_peak file",
                          help="specify path to annotated peaks file", required=True)
    optional.add_argument("-l", "--license", action=licenseAction, metavar="", nargs=0,
                          help='show license status and exit')
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # Perform clip intersection
    intersection_df = license.rmats_eclip_intersect(args.rf, args.ap)
    intersection_df.to_csv('rmats_with_clip.csv', index=False)

    print('Summary complete!')