(defproject findorf "0.0.1-SNAPSHOT"
  :description "ORF finder for BioStar question 5902"
  :dependencies [[org.clojure/clojure "1.2.0"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [org.biojava/biojava3-core "3.0.1"]]
  :repositories {"biojava" "http://www.biojava.org/download/maven/"}
  :run-aliases {:findorf findorf.core})
