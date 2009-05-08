#!/usr/bin/env python
"""Get frequencies from a Tokyo Tyrant remote database.

Total time  : 3:20.37s
Memory      : 6706
Percent CPU : 35.4%
"""
import sys
import random
import pytyrant
import json

def main():
    db = pytyrant.PyTyrant.open("mothra", 1978)

    max_records = 2810718
    num_trials = 500000
    for index in range(num_trials):
        read_id = str(random.randint(0, max_records))
        freq = int(json.loads(db[read_id])['frequency'])
        if index % 10000 == 0:
            print index, freq, read_id

if __name__ == "__main__":
    main()
