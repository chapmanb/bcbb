"""Script to extract the version of software that is run in the analysis pipeline.

This information is meant to be included in the report file and is crucial for
e.g. reproducing and publishing the results.

Usage:
    software_versions.py <post-processing config file>

The software used in the pipeline and which will be included in the report are specified in 
a post-processing configuration file, under the 'program' section. In case the program is a 
python script rather than an executable, the git commit hash will be reported.
"""
import os
import sys
import yaml
from optparse import OptionParser
from bcbio.templates.version import get_versions

def main(config_file):

    # Parse the config yaml file
    if not os.path.exists(config_file):
        log.info("Could not find specified yaml configuration file %s." % config_file)
        return
    
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    for name,ver in get_versions(config).items():
        print "%s: %s" % (name,ver)

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)

