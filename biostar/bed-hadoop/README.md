[Cascalog][1] based approach for [BioStar question][2] about
comparing intervals from two BED files stored in HDFS. Install:

* [Hadoop][3]
* [Leiningen][4]

Then run:

        % lein deps
        % lein uberjar
        % hadoop fs -mkdir /tmp/bed-hadoop/bed-1
        % hadoop fs -mkdir /tmp/bed-hadoop/bed-2
        % hadoop fs -put test/one.bed /tmp/bed-hadoop/bed-1
        % hadoop fs -put test/two.bed /tmp/bed-hadoop/bed-2
        % hadoop jar bed-hadoop-0.0.1-SNAPSHOT-standalone.jar
                     bed_hadoop.core /tmp/bed-hadoop/bed-1 /tmp/bed-hadoop/bed-2
        RESULTS
        ----------
        chr1 20 30
        chr2 40 50
        ----------

[1]: http://github.com/nathanmarz/cascalog
[2]: http://biostar.stackexchange.com/questions/8821/hadoop-genomic-segments-and-join
[3]: http://www.cloudera.com/hadoop/
[4]: https://github.com/technomancy/leiningen#readme
