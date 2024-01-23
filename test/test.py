import argparse
import sys
import os

###############################################################################
# Setup the parser
###############################################################################
parser = argparse.ArgumentParser(description="Run FEN test cases utility.")
parser.add_argument("-ts", "--test-suite", type=str, 
    default="", help="Select the test case suite to run.",)
parser.add_argument("-tc", "--test-case", type=str, 
    default="", help="Select the test case to run.",)
parser.add_argument("-s", "--silent", default=False, action="store_true",
    help="Do not pop up figures.")
parser.add_argument("-cp", "--check-plot", default=False, action="store_true",
    help="Clean the test case directories.")
parser.add_argument("-c", "--clean", default=False, action="store_true",
    help="Clean the test case directories.")

# Collect args
args = parser.parse_args()

###############################################################################
# Clean
###############################################################################
# If clean is selected, perform cleaning operations and exit
if args.clean:
    os.system("cd small_test && bash clean_all.sh")
    sys.exit()

###############################################################################
# List of available run functions
###############################################################################
# Run all small cases test suite
def run_small_test_suite(s):
    os.system("cd small_test && bash run_all.sh" + s)

# Run only Navier-Stokes test cases
def run_navier_stokes_suite(s):
    os.system("cd small_test && cd navier_stokes && bash run.sh" + s)

# Run only Volume-of-Fluid test cases
def run_volume_of_fluid_suite(s):
    os.system("cd small_test && cd volume_of_fluid && bash run.sh" + s)

# Run only Immersed-Boundary test cases
def run_eulerian_solid_suite(s):
    os.system("cd small_test && cd eulerian_solid && bash run.sh" + s)

# Run only Immersed-Boundary test cases
def run_deformable_solid_suite(s):
    os.system("cd small_test && cd deformable_solid && bash run.sh" + s)

# Run single test case
def run_single_test_case(test_case: str, s: str):
    for root, dirs, files in os.walk('.', topdown=True):
        for d in dirs:
            if test_case in d:
                dirname = os.path.join(root, d)
    os.system("cd " + dirname + " && bash run.sh" + s)

###############################################################################
# Based on input select the run function
###############################################################################
def execute():
    # Optional argument for silent run
    s = ''
    if args.silent:
        s = ' silent'

    if args.test_suite == '':
        exit
    elif args.test_suite == 'navier_stokes' or args.test_suite == 'NS':
        run_navier_stokes_suite(s)
    elif args.test_suite == 'volume_of_fluid' or args.test_suite == 'VoF':
        run_volume_of_fluid_suite(s)
    elif args.test_suite == 'eulerian_solid' or args.test_suite == 'ES':
        run_eulerian_solid_suite(s)
    elif args.test_suite == 'deformable_solid' or args.test_suite == 'DF':
        run_deformable_solid_suite(s)
    elif args.test_suite == 'small_test_suite' or args.test_suite == 'STS':
        run_small_test_suite(s)
    else:
        raise Exception('Wrong test case suite selected. Available are:\n\
                        navier_stokes (NS) \n\
                        volume_of_fluid (VoF) \n\
                        eulerian_solid (ES)\n\
                        small_test_suite (STS).')

    if args.test_case == '':
        exit
    else:
        run_single_test_case(args.test_case, s)

    if args.check_plot:
        os.system('eom small_test/*/*/*.png')

if __name__ == '__main__':
    execute()
