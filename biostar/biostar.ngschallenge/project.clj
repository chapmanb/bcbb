(defproject biostar.ngschallenge "0.0.1-SNAPSHOT"
  :description "Custom code for BioStar NGS challenge"
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [bcbio.variation "0.0.1-SNAPSHOT"]]
  :run-aliases {:find biostar.ngschallenge.core})
