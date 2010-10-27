This directory contains details of the software installed with
[CloudBioLinux][1]. This is the right place to dig around if you are interested
in adding packages to the image, or would like to get an overview of what is
installed. The configuration files are written in easily readable [YAML format][2].

* `main.yaml` --  High level category view of packages and libraries that are
  installed with CloudBioLinux.

* `packages.yaml` -- A full list of operating system packages that are included,
  organized by category. The names are standard [Ubuntu APT package names][3].

* `custom.yaml` -- Installed software that is not included in the standard
  package repository. These are often specialized biological packages that have
  not yet been cleanly packaged. Actual installation code is in the `custom`
  sub-directory.

* `python-libs.yaml`, `r-libs.yaml`, `perl-libs.yaml`, `ruby-libs.yaml`,
  `clojure-libs.yaml` -- Libraries installed for a number of programming
  languages. These are installed by the language specific library managers
  (easy\_install for Python, cpan for Perl, gem for Ruby).


[1]: http://cloudbiolinux.com/
[2]: http://en.wikipedia.org/wiki/YAML
[3]: https://help.ubuntu.com/community/AptGet/Howto
