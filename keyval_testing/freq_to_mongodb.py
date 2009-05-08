#!/usr/bin/env python
"""Write a file of read frequencies to a MongoDB database.

./mongod --dbpath /store3/alt_home/chapman/dbs/mongodb --quiet

Original stats:
Total time  : 5:13.66s
Memory      : 4954
Percent CPU : 99.9%

    -rwxrwxr-x 1 chapman users    0 2009-05-05 17:58 mongod.lock
    -rw------- 1 chapman users  64M 2009-05-05 17:58 reads_090504.0
    -rw------- 1 chapman users 128M 2009-05-05 17:59 reads_090504.1
    -rw------- 1 chapman users 256M 2009-05-05 17:59 reads_090504.2
    -rw------- 1 chapman users 512M 2009-05-05 18:01 reads_090504.3
    -rw------- 1 chapman users 1.0G 2009-05-05 18:03 reads_090504.4
    -rw------- 1 chapman users  16M 2009-05-05 17:58 reads_090504.ns

With _id change:

Total time  : 2:58.90s
Memory      : 4988
Percent CPU : 99.9%

-rw------- 1 chapman users  64M 2009-05-08 10:19 reads_090504.0
-rw------- 1 chapman users 128M 2009-05-08 10:20 reads_090504.1
-rw------- 1 chapman users 256M 2009-05-08 10:20 reads_090504.2
-rw------- 1 chapman users 512M 2009-05-08 10:21 reads_090504.3
-rw------- 1 chapman users  16M 2009-05-08 10:19 reads_090504.ns
"""
import sys
from pymongo.connection import Connection
from pymongo import ASCENDING

def main(in_file):
    conn = Connection("mothra")
    db = conn["reads_090504"]
    col = db["read_to_freq"]

    with open(in_file) as in_handle:
        for read_index, freq in enumerate(in_handle):
            col.insert(dict(_id=read_index, freq=int(freq)))
    #col.create_index("_id", ASCENDING)

if __name__ == "__main__":
    main(sys.argv[1])
