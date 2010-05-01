"""Scrape the Biolinux website to retrieve a list of packages they install.

http://www.jcvi.org/cms/research/projects/jcvi-cloud-biolinux/included-software

This needs to run on a machine with an apt system to check for the existance of
package names.
"""
import sys
import urllib2
import re
import subprocess
import StringIO

from BeautifulSoup import BeautifulSoup

def main():
    url = "http://www.jcvi.org/cms/research/projects/jcvi-cloud-biolinux/included-software"
    in_handle = urllib2.urlopen(url)
    soup = BeautifulSoup(in_handle)
    tables = soup.findAll("table", {"class": "contenttable"})
    to_check = []
    for t in tables:
        for row in soup.findAll("tr", {"class" : re.compile("tableRow.*")}):
            for i, item in enumerate(row.findAll("p", {"class": "bodytext"})):
                if i == 0:
                    to_check.append(str(item.contents[0]))
    to_check = list(set(to_check))
    packages = [get_package(n) for n in to_check]
    not_ported = [to_check[i] for i, p in enumerate(packages) if p is None]
    packages = [p for p in packages if p]
    print len(to_check), len(packages)
    with open("biolinux-packages.txt", "w") as out_handle:
        out_handle.write("\n".join(sorted(packages)))
    with open("biolinux-missing.txt", "w") as out_handle:
        out_handle.write("\n".join(sorted(not_ported)))

def get_package(pname):
    """Try and retrieve a standard or biolinux package for the package name.
    """
    # custom hacking for painfully general names that take forever
    if pname in ["act", "documentation"]:
        pname = "bio-linux-%s" % pname
    print 'In', pname
    cl = subprocess.Popen(["apt-cache", "search", pname], stdout=subprocess.PIPE)
    cl.wait()
    for line in cl.stdout.read().split():
        package = line.split()[0]
        if package == pname or package == "bio-linux-%s" % pname:
            print 'Out', package
            return package
    return None

if __name__ == "__main__":
    main(*sys.argv[1:])
