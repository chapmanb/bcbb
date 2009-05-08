#!/usr/bin/env python
"""Get frequencies from a MongoDB remote database.

Total time  : 3:57.29s
Memory      : 7848
Percent CPU : 34.6%
"""
import sys
import random
from pymongo.connection import Connection

def main():
    conn = Connection("mothra")
    db = conn["reads_090504"]
    print db.validate_collection("read_to_freq")
    col = db["read_to_freq"]
    print col.index_information()
    print col.options()

    max_records = 2810718
    num_trials = 500000
    for index in range(num_trials):
        read_id = random.randint(0, max_records)
        doc = col.find_one(dict(_id=read_id))
        if doc is None:
            print index, read_id
        freq = int(doc["freq"])
        if index % 10000 == 0:
            print index, freq, doc

if __name__ == "__main__":
    main()
