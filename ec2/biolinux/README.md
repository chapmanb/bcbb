CloudBioLinux provides a [Fabric build file][3] which will install a
large selection of Bioinformatics and machine learning libraries on a
bare machine. This is ideally designed for Amazon EC2 or VirtualBox,
where you start with a bare bones system and bootstrap to a ready to
go instance. Packages to install are fully customizable, and by
default include a large suite of bioinformatics tools and libraries.

# Using an instance

## Amazon

See the 'Getting Started with CloudBioLinux' guide on the
[CloudBioLinux website][1] for a detailed description. The short
version for users familiar with Amazon is:

* Login to the [Amazon EC2 console][2].
* Click Launch Instance, and choose the latest CloudBioLinux AMI from
  the [website][1] in the community AMI section.
* After finishing the configuration wizard, find the host details of
  your running instance from the Instances section.
* Connect to your machine via ssh or VNC.

## VirtualBox with vagrant

Install [VirtualBox 4.0][v2] and [vagrant][v1]. Then add the pre-built
VirtualBox and start it up:

        vagrant box add biolinux_version https://s3.amazonaws.com/chapmanb/biolinux_version.box
        mkdir tmp/biolinux
        cd tmp/biolinux
        vagrant init biolinux_version
        vagrant up

You now have a running virtual machine with CloudBioLinux, and can
connect to it with:

        vagrant ssh

# Building an image from scratch

## Amazon

Install [Fabric][3]:

        sudo apt-get install python-setuptools
        sudo easy_install fabric

Start an Amazon EC2 ubuntu machine using an [Alestic Ubuntu image][4]
as the base, then install CloudBioLinux on it:

        fab -f fabfile.py -H hostname -i private_key_file install_biolinux

## VirtualBox with vagrant

Add a base image and boot it up. Ideally this will eventually be a
64-bit Maverick box, but those aren't available for the latest
Vagrant and VirtualBox:

        vagrant box add lucid32 http://files.vagrantup.com/lucid32.box
        mkdir tmp/biolinux
        cd tmp/biolinux
        vagrant init lucid32
        vagrant up

Run the fabfile, building CloudBioLinux:

        fab -H vagrant -f /path/to/bcbb/ec2/biolinux/fabfile.py install_biolinux

Then build the box, renaming package.box to `cloudbiolinux_date` and
move it to a public webserver, like Amazon S3:

        vagrant package
        mv package.box biolinux_20110122.box
        s3cmd put --acl-public --guess-mime-type biolinux_20110122.box
              s3://chapmanb/biolinux_20110122.box

[1]: http://cloudbiolinux.com/
[2]: https://console.aws.amazon.com/ec2/home
[3]: http://fabfile.org/
[4]: http://alestic.com/
[v1]: http://vagrantup.com/
[v2]: http://digitizor.com/2011/01/07/virtualbox-4-0-install-ubuntu/

# Technical details for using build scripts

## Targets for local builds

The Fabric build files are useful for automating installation of
scientific software on local systems as well as Amazon cloud
servers. There are a number of build targets that can be used to
install a subset of the total available packages.

The main build command is used as described above, but the host
can point to a local machine or any server on your network:

      fab -f fabfile.py -H localhost -c config/fabricrc.txt install_biolinux

The `config/fabricrc.txt` configuration file specifies install
directories and other server specific details.

With this basic command, you can substitute `install_biolinux` with
serveral more specific targets:

* `install_biolinux:packages` -- Install all of the defined system
  packages.
* `install_biolinux:libraries` -- Install all libraries for various
  programming languages.
* `install_libraries:language` -- Install libraries for a specific
  language.
* `install_biolinux:custom` -- Install all custom programs.
* `install_custom:your_package_name` -- Install a specific custom
   program.

### Custom package installs

The custom directory contains installation instructions for programs that are
not available from standard package repositories. These instructions are written
in Python using the [Fabric][3] remote deployment tool and can also be used for
installing individual packages locally on your machine. To do this, run:

      fab -f fabfile.py -H localhost install_custom:your_package_name

To build and install `your_package_name` on the local machine. We welcome
additional custom bioinformatics package definitions for inclusion in
CloudBioLinux. `custom/shared.py` contains a number of higher level functions
which make it easier to write installation instructions.

## EC2 quickstart
This provides a quick cheat sheet of commands for getting up and running on EC2 using
Amazon's command line tools.

### Initial set up

The first time using EC2, you'll need to install the toolkit and credentials
for connecting. Follow these basic directions:
<http://docs.amazonwebservices.com/AWSEC2/latest/GettingStartedGuide/>

Login to your Amazon EC2 account (<http://aws.amazon.com/account/>) and go to
Security Credentials/X.509. Create a new certificate and download the public
`cert-*.pem` and `private pk-*.pem` files. Put these in `~.ec2`.

Download and unzip the ec2 api tools, which require java:
<http://developer.amazonwebservices.com/connect/entry.jspa?externalID=351&categoryID=88>

Set up .zshrc/.bashrc:

       export EC2_HOME=$HOME/install/ec2/ec2-api-tools
       export EC2_PRIVATE_KEY=~/.ec2/pk-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export EC2_CERT=~/.ec2/cert-UBH43XTAWVNQMIZRAV3RP5IIBAPBIFVP.pem
       export PATH=$PATH:$EC2_HOME/bin

To test, you should be able to run the command:

       % ec2-describe-regions

Now generate a privatekey for logging in:

       % ec2-add-keypair yourmachine-keypair

This will produce an RSA private key. You should copy and paste this to your
.ec2 directory for future use:

       % vim ~/.ec2/id-yourmachine.keypair
       % chmod 600 ~/.ec2/id-yourmachine.keypair

Allow ssh and web access to your instances:

       % ec2-authorize default -p 22
       % ec2-authorize default -p 80

### Starting an instance

Each time you'd like to use EC2, you need to create a remote instance to work
with. This can be done nicely via the AWS console, but if you'd like
to use the commandline, these quick start docs follow
<http://docs.amazonwebservices.com/AWSEC2/latest/GettingStartedGuide/running-an-instance.html>.

Pick an AMI, start an instance and ensure that it is running:

       % ec2-run-instances ami-1ad03273 -k sobchak-keypair
       RESERVATION	r-0a7af462	678711657553	default
       INSTANCE	i-0ca39764	ami-1ad03273			pending

       % ec2-describe-instances i-0ca39764
       RESERVATION	r-0a7af462	678711657553	default
       INSTANCE	i-0ca39764	ami-1ad03273 ec2-174-129-68-135.compute-1.amazonaws.com

Now you can ssh in using the key you created:

       % ssh -i ~/.ec2/id-sobchak.keypair root@ec2-174-129-68-135.compute-1.amazonaws.com

You're in and paying per hour. When done:

       % ec2-terminate-instances i-0ca39764

### Creating EBS stores

Create a 200Gb store:

       ec2-create-volume -z us-east-1c -s 200
       ec2-describe-volumes vol-7568e61c

Establish a store from an existing snapshot:

       ec2-create-volume --snapshot snap-10dbab78 -z us-east-1c

Attach a store to a running instance:

       ec2-attach-volume -d /dev/sdh -i i-351d225e vol-7568e61c

Only for a new instance -- ssh to the machine and create a
filesystem:

       grep -q xfs /proc/filesystems || sudo modprobe xfs
       sudo mkfs.xfs /dev/sdh

Mount the directory, need to do this every time on the machine:

       sudo mkdir -p /mnt/biodata
       sudo mount -t xfs -o noatime /dev/sdh /mnt/biodata/

## Eucalyptus Public Cloud Quickstart

* Sign up for an account at: <https://ecc.eucalyptus.com>

* Download credentials from: <https://ecc.eucalyptus.com/#credentials>.
  Unpack all files to ~/.euca directory

* Include on your .bashrc file so that the environment variables are
  loaded: `source .euca/eucarc`
  NOTE: Eucalyptus environment variables override those from EC2. AWS will
  not work properly.

* Follow the starting an instance directions above with euca2ools:

       euca-add-keypair biolinux
       euca-run-instances emi-CBEA100C -k biolinux
       euca-authorize -P tcp -p 22 -s 0.0.0.0/0 default


## Bundling a new EMI in Eucalyptus

Eucalyptus provides some prebuilt images under the "Extra" panel on the default installation. Those
images should be packaged and uploaded to the Eucalyptus cloud you're working on. The following
is just a packaging example with one of the images, further automation should be done here.

        export PREFIX="/scratch/euca-fedora-11-i386"

        export IMAGE=${PREFIX}/fedora.11.x86.img
        export KERNEL="${PREFIX}/kvm-kernel/vmlinuz-2.6.28-11-server"
        export KERNEL_B="vmlinuz-2.6.28-11-server"

        export INITRD="${PREFIX}/kvm-kernel/initrd.img-2.6.28-11-server"
        export INITRD_B="initrd.img-2.6.28-11-server"

        euca-bundle-image -i $KERNEL --kernel true
        euca-upload-bundle -b kernels -m /tmp/${KERNEL_B}.manifest.xml
        export EKI=`euca-register kernels/${KERNEL_B}.manifest.xml | cut -f1 -d" "`

        euca-bundle-image -i $INITRD --ramdisk true
        euca-upload-bundle -b ramdisks -m /tmp/`basename ${INITRD}`.manifest.xml
        export ERI=`euca-register ramdisks/${INITRD_B}.manifest.xml | cut -f1 -d" "`

        euca-bundle-image -i $IMAGE --kernel $EKI --ramdisk $ERI
        euca-upload-bundle -b cloudbiolinux -m /tmp/`basename ${KERNEL}`.manifest.xml
        euca-register cloudbiolinux/`basename ${KERNEL}`.manifest.xml
