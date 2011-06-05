[Cascalog][1] based approach for [BioStar question][2] about
comparing intervals from two BED files stored in HDFS. To run:

        % hadoop fs -mkdir /tmp/bed-hadoop
        % hadoop fs -mkdir /tmp/bed-hadoop/bed-1
        % hadoop fs -mkdir /tmp/bed-hadoop/bed-2
        % hadoop fs -put one.bed /tmp/bed-hadoop/bed-1
        % hadoop fs -put two.bed /tmp/bed-hadoop/bed-2
        % lein deps
        % lein run :bed-overlap /tmp/bed-hadoop/bed-1 /tmp/bed-hadoop/bed-2
        RESULTS
        ----------
        chr1 20 30
        chr2 40 50
        ----------

[1]: http://github.com/nathanmarz/cascalog
[2]: http://biostar.stackexchange.com/questions/8821/hadoop-genomic-segments-and-join
