#!/usr/bin/python
""" Reads a Illumina sample sheet CSV file, modify values and write a new file.
This is intended for quick samplesheet sanitizing and creation via cmdline.

- The 48 index sequences are extracted from a file provided by Illumina.

- The indexes must be given by 'X_indexN', where 'X' means any string, and 'N'
  means any integer in the range from 1 to 48.

- By default, an extra 'A' is added to the index sequence. Use option '-n'
  to switch this off.

- The options allow setting values for all rows of the columns FCID, SampleRef
  Recipe and Operator.

- If no option is given, any empty value of a column will be replaced by the
  following values: FCID (value in first row), SampleRef ('unknown'),
  Description ('Test'), Control ('N'), Recipe ('R1'), Operator ('NN').

Per Kraulis
"""

import sys, csv

VERSION = '1.1'


INDEX_LOOKUP = dict(index1='ATCACG',
                    index2='CGATGT',
                    index3='TTAGGC',
                    index4='TGACCA',
                    index5='ACAGTG',
                    index6='GCCAAT',
                    index7='CAGATC',
                    index8='ACTTGA',
                    index9='GATCAG',
                    index10='TAGCTT',
                    index11='GGCTAC',
                    index12='CTTGTA',
                    index13='AGTCAA',
                    index14='AGTTCC',
                    index15='ATGTCA',
                    index16='CCGTCC',
                    index17='GTAGAG',
                    index18='GTCCGC',
                    index19='GTGAAA',
                    index20='GTGGCC',
                    index21='GTTTCG',
                    index22='CGTACG',
                    index23='GAGTGG',
                    index24='GGTAGC',
                    index25='CTAGCT',
                    index26='CTAGCT',
                    index27='CTATAC',
                    index28='CTCAGA',
                    index29='CTGCTG',
                    index30='TAATCG',
                    index31='TACAGC',
                    index32='TATAAT',
                    index33='TCATTC',
                    index34='TCCCGA',
                    index35='TCGAAG',
                    index36='TCGGCA',
                    index37='GCCGCG',
                    index38='GCCTTA',
                    index39='GCTCCA',
                    index40='GGCACA',
                    index41='GGCCTG',
                    index42='TCTACC',
                    index43='TGAATG',
                    index44='TGCCAT',
                    index45='TGCTGG',
                    index46='TGGCGC',
                    index47='TTCGAA',
                    index48='TTCTCC')
INDEX_LOOKUP.update(dict([(k.replace('index', ''), v)
                          for k,v in INDEX_LOOKUP.items()]))


class Samplesheet(object):
    """Read a sample sheet CSV file.
    Process the contents to fill in values.
    Write the modified sample sheet CSV file.
    """

    HEADER = ('FCID',
              'Lane',
              'SampleID',
              'SampleRef',
              'Index',
              'Description',
              'Control',
              'Recipe',
              'Operator')

    def __init__(self, silent=True, verbose=False, fcid=None,
                 sampleref=None, extra_a=True, recipe=None, operator=None):
        self.silent = silent
        self.verbose = verbose
        self.fcid = fcid
        self.sampleref = sampleref
        self.extra_a = extra_a
        self.recipe = recipe
        self.operator = operator

    def message(self, message):
        "Output message to stdout."
        if self.silent: return
        if self.verbose:
            print message

    def warning(self, message):
        "Output warning to stdout."
        if self.silent: return
        print 'WARNING:', message

    def error(self, message, quit=True):
        "Output a message to stderr, and optionally quit."
        if not self.silent:
            sys.stderr.write(message)
            sys.stderr.write('\n')
        if quit:
            sys.exit(1)

    def read(self, filename):
        "Read the CSV file with the given name."
        try:
            infile = open(filename, mode='rU')
        except IOError, msg:
            self.error("could not open file '%s': %s\n" % (filename, msg))
        dialect = csv.Sniffer().sniff(infile.read(1024), delimiters=',;\t')
        infile.seek(0)
        reader = csv.reader(infile, dialect)
        self.rows = list(reader)
        infile.close()
        if not self.rows:
            self.error('no rows in input file')
        self.message("read %s rows from %s" %(len(self.rows), filename))

    def process(self):
        "Process the CSV table contents, checking and modifying as required."
        self.process_header()
        self.process_row_length()
        self.process_fcid()
        self.process_lane()
        self.process_sampleref()
        self.process_index()
        self.process_description()
        self.process_control()
        self.process_recipe()
        self.process_operator()

    def process_header(self):
        "Remove the header, if any, from the rows. Set the output header."
        self.header = self.HEADER
        header = self.rows.pop(0)
        if header[0].strip().upper() == 'FCID':
            self.header = header
        else:
            self.rows.insert(0, header)
        if not self.rows:
            self.error('no data rows in input file')

    def process_row_length(self):
        "Make sure that rows are at least 9 items long."
        for row in self.rows:
            if len(row) < 9:
                row.extend([''] * (9 - len(row)))

    def process_fcid(self):
        """Set the FCID column values, if not specified.
        Set the value if FCID was specified.
        Else replace empty values by the value of the first row.
        """
        if self.fcid:
            for row in self.rows:
                row[0] = self.fcid
        else:
            fcid = self.rows[0][0]
            if not fcid:
                self.error('no FCID specified')
            for row in self.rows[1:]:
                if not row[0]:
                    row[0] = fcid

    def process_lane(self):
        "Convert the lane values to integer."
        self.lanes = dict()
        for i, row in enumerate(self.rows):
            try:
                lane = int(row[1])
            except ValueError, msg:
                self.error("row %i: lane '%s' not a number" % (i+1, row[1]))
            else:
                row[1] = lane
                self.lanes.setdefault(lane, set())

    def process_sampleref(self):
        "Set the SampleRef column values, or 'unknown' if empty value"
        if self.sampleref:
            for row in self.rows:
                row[3] = self.sampleref
        else:
            for row in self.rows:
                if not row[3]:
                    row[3] = 'unknown'

    def process_index(self):
        "Set the index sequence; the index label from the SampleID column."
        for i, row in enumerate(self.rows):
            label = row[2].split('_')[-1].lower()
            if label in self.lanes[row[1]]:
                self.warning("row %i: index '%s' used multiple times, lane '%i'"
                             % (i+1, label, row[1]))
            else:
                self.lanes[row[1]].add(label)
            try:
                sequence = INDEX_LOOKUP[label]
            except KeyError:
                self.error("row %i: no such index '%s'" % (i+1, label))
            else:
                if self.extra_a:
                    sequence += 'A'
                if row[4] != sequence:
                    if row[4]:
                        self.warning("row %i: replacing index sequence '%s'"
                                     % (i+1, row[4]))
                    row[4] = sequence

    def process_description(self):
        "Set the description to 'Test' if empty value."
        for row in self.rows:
            if not row[5]:
                row[5] = 'Test'

    def process_control(self):
        "Set the control to 'N' if empty value."
        for row in self.rows:
            if not row[6]:
                row[6] = 'N'

    def process_recipe(self):
        "Set the Recipe column values, or 'R1' if empty value."
        if self.recipe:
            for row in self.rows:
                row[7] = self.recipe
        else:
            for row in self.rows:
                if not row[7]:
                    row[7] = 'R1'

    def process_operator(self):
        "Set the Operator column values, or 'NN' if empty value."
        if self.operator:
            for row in self.rows:
                row[8] = self.operator
        else:
            for row in self.rows:
                if not row[8]:
                    row[8] = 'NN'

    def write(self, filename):
        "Write the CSV file with the given name."
        try:
            outfile = open(filename, mode='wb')
        except IOError, msg:
            self.error("could not create file '%s': %s\n" % (filename, msg))
        writer = csv.writer(outfile, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(self.header)
        writer.writerows(self.rows)
        outfile.close()
        self.message("wrote %s rows to %s" %(len(self.rows)+1, filename))


if __name__ == '__main__':
    import optparse
    usage = '%prog [options] infile outfile'
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % VERSION)
    parser.add_option('-s', '--silent', action='store_true', default=False,
                      help='no information output at all')
    parser.add_option('-v', '--verbose', action='store_true', default=False,
                      help='output more information')
    parser.add_option('-f', '--fcid', action='store', type='string',
                      help='set FCID for all rows')
    parser.add_option('-r', '--sampleref', action='store', type='string',
                      help='set SampleRef for all rows')
    parser.add_option('-a', '--extra_a', action='store_true', default=True,
                      help="add an extra A to sequence (default)")
    parser.add_option('-n', '--no_extra_a', action='store_false',dest='extra_a',
                      help="do not add an extra A to sequence")
    parser.add_option('-p', '--recipe', action='store', type='string',
                      help='set Recipe for all rows')
    parser.add_option('-o', '--operator', action='store', type='string',
                      help='set Operator for all rows')
    options, args = parser.parse_args()
    if len(args) != 2:
        parser.print_usage(sys.stderr)
        sys.exit(1)
    samplesheet = Samplesheet(silent=options.silent,
                              verbose=options.verbose,
                              fcid=options.fcid,
                              sampleref=options.sampleref,
                              extra_a=options.extra_a,
                              recipe=options.recipe,
                              operator=options.operator)
    samplesheet.read(args[0])
    samplesheet.process()
    samplesheet.write(args[1])
