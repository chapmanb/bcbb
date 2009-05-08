#!/usr/bin/env python
"""Get frequencies from a CouchDB remote database.

Total time  : 14:08.47s
Memory      : 8112
Percent CPU : 63.2%
"""
import sys
import random
import couchdb.client

def main():
    db_name = "reads_090504/read_to_freq"
    server = couchdb.client.Server("http://mothra:5984/")
    db = server[db_name]

    max_records = 2810718
    num_trials = 500000
    misses = 0
    for index in range(num_trials):
        read_id = str(random.randint(0, max_records))
        try:
            doc = db[read_id]
        except:
            misses += 1
        freq = int(doc["frequency"])
        if index % 10000 == 0:
            print index, misses, read_id, freq

if __name__ == "__main__":
    main()
