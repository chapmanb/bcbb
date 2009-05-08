"""Write file of read frequencies to a Tokyo Tyrant server.

ttserver test.tcb#opts=ld#bnum=1000000#lcnum=10000

Total time  : 11:58.90s
Memory      : 2536
Percent CPU : 39.2%

-rw-r--r-- 1 chapman users 24M 2009-05-07 09:01 test.tcb
"""
import sys
import pytyrant
import json

def main(in_file):
    db = pytyrant.PyTyrant.open("mothra", 1978)

    with open(in_file) as in_handle:
        for read_index, freq in enumerate(in_handle):
            db[str(read_index)] = json.dumps(dict(frequency=freq))

if __name__ == "__main__":
    main(sys.argv[1])
