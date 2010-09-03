import boto
import collections

OWNER = '678711657553' # Brad's owner ID

def images_and_snapshots(owner):
    """Retrieve Biolinux image and snapshot information.
    """
    conn = boto.connect_ec2()
    images = conn.get_all_images(owners=[owner])
    images32 = _sorted_images(images, "CloudBioLinux 32")
    images64 = _sorted_images(images, "CloudBioLinux 64")
    datalibs = _data_libraries(conn, owner)
    print images32
    print images64
    print datalibs

def _data_libraries(conn, owner):
    library_types = collections.defaultdict(list)
    snaps = conn.get_all_snapshots(owner=owner)
    for snap in snaps:
        if snap.description.startswith("CloudBioLinux Data"):
            # the type is everything except the start and date
            data_type = " ".join(snap.description.split()[2:-1])
            library_types[data_type].append(snap)
    final = dict()
    for name, snaps in library_types.iteritems():
        snaps = [(s.description, s) for s in snaps]
        snaps.sort(reverse=True)
        final[name] = [(s.id, d) for (d, s) in snaps]
    return final

def _sorted_images(images, start_name):
    """Retrieve a sorted list of images with most recent first.
    """
    images = [(i.name, i) for i in images if i.name.startswith(start_name)]
    images.sort(reverse=True)
    return [(i.id, name) for (name, i) in images]

images_and_snapshots(OWNER)
