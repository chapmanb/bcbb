"""Install instructions for non-packaged java programs.
"""
import os

from fabric.api import *
from fabric.contrib.files import *

from shared import _if_not_installed

@_if_not_installed("lein")
def install_leinengin(env):
    """Standard clojure build tool: http://github.com/technomancy/leiningen
    """
    run("wget http://github.com/technomancy/leiningen/raw/stable/bin/lein")
    run("chmod a+rwx lein")
    sudo("mv lein %s" % os.path.join(env.system_install, "bin"))
    run("lein self-install")

def install_incanter(env):
    """Clojure based statistics and graphics environment.

    http://github.com/liebke/incanter
    """
    # Can we handle requirements cleaner and keep YAML config style?
    install_leinengin(env)

    clojure_dir = os.path.join(env.local_install, "clojure")
    install_dir = os.path.join(clojure_dir, "incanter")
    if not exists(install_dir):
        if not exists(clojure_dir):
            run("mkdir %s" % clojure_dir)
        with cd(clojure_dir):
            run("git clone git://github.com/liebke/incanter.git")
        with cd(install_dir):
            run("lein deps")
