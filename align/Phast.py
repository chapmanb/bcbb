"""Wrappers for PHAST applications: conservation, alignments and phylogeny.

http://compgen.bscb.cornell.edu/phast/
"""
import types

from Bio.Application import _Option, _Argument, _Switch, AbstractCommandline

class PhastConsCommandline(AbstractCommandline):
    def __init__(self, cmd="phastCons", **kwargs):
        self.parameters = [
            _Argument(["alignment"], ["input"],
                    None, True, ""),
            _Argument(["models"], ["input"],
                    None, True, ""),
            _Option(["--target-coverage"], ["input"],
                    None, False, "", False),
            _Option(["--expected-length"], ["input"],
                    None, False, "", False),
            _Option(["--rho"], ["input"],
                    None, False, "", False),
            _Option(["--msa-format"], ["input"],
                    None, False, "", False),
            _Option(["--estimate-trees"], ["input"],
                    None, False, "", False),
            _Switch(["--no-post-probs"], ["input"]),
            _Option(["--most-conserved"], ["input"],
                    None, False, "", False),
            _Option(["--estimate-rho"], ["input"],
                    None, False, "", False),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
