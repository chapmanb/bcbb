#!/bin/bash

#  See this blog:
#  http://alestic.com/2009/08/runurl
#
#  The basic Alestic images support execution of the user data, and also
#  come with "runurl" pre-installed.  Therefore we should be able to simply stick the
#  following into user data:
#  
#  #!/bin/bash -ex
#  runurl cloudbiolinux.com/bootstrap
#
#  Of to support non-alestic images:
#  
#  #!/bin/bash -ex
#  which runurl || wget -O- run.alestic.com/install/runurl | bash
#  runurl cloudbiolinux.com/bootstrap

#  Now the URL given would actually be a redirect to:
#  http://github.com/chapmanb/bcbb/raw/master/ec2/biolinux/bootstrap.sh

if [ `id -un` != root ] ; then
	echo "This script needs to be run as root."
        exit 1
fi
HOME=~root; USER=root

#  So, here we go...
tmpdir=/tmp/bio-linux-bootstrap
#Skip this bit if the directory was found - eases debugging
if [ ! -e "$tmpdir" ] ; then

	#  Install git and fabric and python stuff and ssh
	export DEBIAN_FRONTEND=noninteractive
	apt-get update
	apt-get -y install git-core fabric python-setuptools python-yaml openssh-server

	# Make a directory to work in
	mkdir $tmpdir ; cd $tmpdir

	# Pull the GIT stuff - this is wasteful but I don't know a better way.
	git clone --depth=0 http://github.com/chapmanb/bcbb.git .
	ls -A | grep -v ec2 | xargs rm -r

	#Some fixups to the fabfile - these can probably be changed in the upstream but
        #I want my script to work now...
	if ! which curl ; then
		echo "Removing use of curl from fabfile"
		sed -i -e 's/curl -s /wget -q -O- /' ec2/biolinux/fabfile.py
	fi
	if ! id -un ubuntu ; then
		echo "Removing env.user=ubuntu from fabfile"
		sed -i -e '/env.user = "ubuntu"/d' ec2/biolinux/fabfile.py
	fi
fi

cd "$tmpdir"/ec2/biolinux/

# Fabric relies on an ssh connection, but this script is local and running as root,
# so I need to allow root logins with SSH and get myself a password-less login.
# What a faff
ssh-keygen -t rsa -N '' -C 'localhost passwordless login key' -f ~/.ssh/id_local_rsa
( echo -n 'from="localhost,::1,127.0.0.0/24" ' ; cat ~/.ssh/id_local_rsa.pub ) >> ~/.ssh/authorized_keys2
/usr/sbin/sshd -o PermitRootLogin=without-password -o ListenAddress=localhost:2232 -o AllowGroups=root -o PidFile=/var/run/sshd_local.pid -p 2232

# Run Fabric on localhost
#fab -i ~/.ssh/id_local_rsa -H localhost:2232 uname_a
fab -i ~/.ssh/id_local_rsa -H localhost:2232 install_biolinux

# Clean up keys and sshd
rm ~/.ssh/id_local_rsa.pub ~/.ssh/id_local_rsa ~/.ssh/authorized_keys2
xargs kill -TERM < /var/run/sshd_local.pid

# Say what we did
echo "All finished.  You probably want to reboot now."
