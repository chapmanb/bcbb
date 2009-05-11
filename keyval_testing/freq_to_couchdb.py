#!/usr/bin/env python
"""Write a file of read frequencies to a CouchDB database.

Did not finish on first pass with naive code:

    Total time  : 22:15:19.10s
    Memory      : 4896
    Percent CPU : 2.6%

    -rw-rw-r-- 1 chapman users 5.7G 2009-05-07 07:52 read_to_freq.couch

    1842694 documents loaded

Re-done with improvements from Chris and Paul for bulk loading:

    Total time  : 5:31.01s
    Memory      : 7490
    Percent CPU : 10.4%

    -rw-rw-r-- 1 chapman users 236M 2009-05-11 08:01 read_to_freq.couch
"""
import sys
import couchdb.client

def main(in_file):
    db_name = "reads_090504/read_to_freq"
    server = couchdb.client.Server("http://mothra:5984/")
    # tune this somewhere between 500-5000 depending on your doc size
    bulk_size = 5000
    bulk_docs = []
    if db_name in server:
        db = server[db_name]
    else:
        db = server.create(db_name)

    with open(in_file) as in_handle:
        for read_index, freq in enumerate(in_handle):
            bulk_docs.append(dict(_id=str(read_index), frequency=int(freq)))
            if len(bulk_docs) >= bulk_size:
                db.update(bulk_docs)
                bulk_docs = []
        db.update(bulk_docs)

if __name__ == "__main__":
    main(sys.argv[1])
