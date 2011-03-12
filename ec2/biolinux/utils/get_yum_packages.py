"""Convert list of apt packages to matching yum packages.

This needs to run on a machine with yum to check for the existance of
package names.
"""
import os
import re
import sys
import subprocess
import platform
from contextlib import nested
import StringIO

def main(orig_file):
    new_file = "%s-yum%s" % os.path.splitext(orig_file)
    with nested(open(orig_file), open(new_file, "w")) as \
               (orig_handle, new_handle):
        for line in orig_handle:
            if line.lstrip().startswith("- "):
                base, orig_package = line.split("- ")
                yum_package = get_yum_package(orig_package.strip())
                if yum_package:
                    new_handle.write("%s- %s\n" % (base, yum_package))
            else:
                new_handle.write(line)

def get_yum_package(pname):
    print 'In', pname
    # hacks for package names that cause it to hang
    if pname in ["ri"]:
        return None
    elif pname in ["perl"]:
        return pname
    cl = subprocess.Popen(["yum", "search", pname], stdout=subprocess.PIPE)
    cl.wait()
    arch_pname = "%s.%s" % (pname, platform.machine())
    for line in cl.stdout.read().split("\n"):
        if line.startswith(arch_pname):
            return pname
    return None

if __name__ == "__main__":
    main(*sys.argv[1:])
