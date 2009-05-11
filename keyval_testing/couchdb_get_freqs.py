#!/usr/bin/env python
"""Get frequencies from a CouchDB remote database.

Total time  : 14:22.19s
Memory      : 8154
Percent CPU : 63.8%
"""
import sys
import random
import couchdb.client

def main():
    db_name = "reads_090504/read_to_freq"
    server = couchdb.client.Server("http://mothra:5984/")
    db = server[db_name]

    max_records = 2810717
    num_trials = 500000
    for index in range(num_trials):
        read_id = str(random.randint(0, max_records))
        doc = db[read_id]
        freq = int(doc["frequency"])
        if index % 10000 == 0:
            print index, read_id, freq

if __name__ == "__main__":
    main()
