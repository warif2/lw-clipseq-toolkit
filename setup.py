import sys, os
binpath = '/'.join(os.path.realpath(__file__).split("/")[:-1])+'/'
sys.path.insert(1, binpath)
import license_setup

if __name__ == '__main__':
    license_setup.license_setup()