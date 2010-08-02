#!/bin/sh

firstboot_file=/var/log/firstboot.done 
if [ -e $firstboot_file ]; then
    return 1;
fi
#Find out if a password was given
USERINFO=`wget -q -O- http://169.254.169.254/2009-04-04/user-data | grep ^USER1PASS=`
USER1PASS=${USERINFO#USER1PASS=}

#or, more trusting
#wget -q -O- http://169.254.169.254/2009-04-04/user-data > /tmp/user-data
#source /tmp/user-data

if [ -n "$USER1PASS" ]; then
    useradd -m -s /bin/bash user1
    chpasswd "user1:$USER1PASS" user1

    #force SSH to allow password logins
    sed -i 's/^PasswordAuthentication .*/PasswordAuthentication yes/' /etc/ssh/sshd_config
fi

#force freenx-server to fix itself in any case
dpkg-reconfigure -pcritical freenx-server

echo 'firstboot' > $firstboot_file
