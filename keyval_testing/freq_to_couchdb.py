#!/usr/bin/env python
"""Write a file of read frequencies to a CouchDB database.

Did not finish:

    Total time  : 22:15:19.10s
    Memory      : 4896
    Percent CPU : 2.6%

    -rw-rw-r-- 1 chapman users 5.7G 2009-05-07 07:52 read_to_freq.couch

    1842694 documents loaded
"""
import sys
import couchdb.client

def main(in_file):
    db_name = "reads_090504/read_to_freq"
    server = couchdb.client.Server("http://mothra:5984/")
    if db_name in server:
        db = server[db_name]
    else:
        db = server.create(db_name)

    with open(in_file) as in_handle:
        for read_index, freq in enumerate(in_handle):
            db[str(read_index)] = dict(read_id=read_index, frequency=freq)

if __name__ == "__main__":
    main(sys.argv[1])
