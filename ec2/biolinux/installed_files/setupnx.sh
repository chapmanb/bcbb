#!/bin/bash

# Set up NX to work with password access.
# This can be run manually or triggered from ~/.bash_login

if [ `id -u` != 0 ] ; then
    echo "This script must be run as root."
    exit 1
fi

CURRENTUSER=${1:-$USER}

cat <<BLURB
The NX remote desktop software needs password-based authentication enabled.
This functionality will now be activated.  While you can still log in using your private
key file, please note that password authentication will now be possible both via NX and 
regular ssh for any user account which has a password set.

You can choose to set a password for $CURRENTUSER or give a new user name to create a
new non-privileged account.

BLURB

read -p "User to set password for [$CURRENTUSER]: " USERTOMAKE
if [ -z "$USERTOMAKE" ] ; then USERTOMAKE="$CURRENTUSER" ; fi
if getent passwd "$USERTOMAKE" > /dev/null ; then
    passwd "$USERTOMAKE"
else
    useradd -m -s /bin/bash "$USERTOMAKE"
    passwd "$USERTOMAKE"
fi

#force SSH to allow password logins
echo "Reconfiguring SSH to allow password login."
sed -i 's/^PasswordAuthentication .*/PasswordAuthentication yes/' /etc/ssh/sshd_config
/etc/init.d/ssh reload

echo "Setting up NX."
#force freenx-server to fix itself in any case
if ! dpkg -s freenx-server >&/dev/null ; then
    apt-get -y install freenx-server
fi
dpkg-reconfigure -pcritical freenx-server >&/dev/null

echo "DONE!  You can find clients for NX here: "
echo "               http://www.nomachine.com/download.php"

