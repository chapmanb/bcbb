(defproject bed-hadoop "0.0.1-SNAPSHOT"
  :description "BED processing with cascalog on Hadoop"
  :dependencies [[org.clojure/clojure "1.2.1"]
                 [cascalog "1.8.0-SNAPSHOT"]]
  :dev-dependencies [[org.apache.hadoop/hadoop-core "0.20.2-dev"]
                     [swank-clojure "1.3.1"]]
  :run-aliases {:bed-overlap bed-hadoop.core}
  :aot [bed-hadoop.core])
