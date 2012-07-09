(ns clj-genomespace.test-core
  "Basic usage tests for GenomeSpace integration.
  Requires setting GS_USERNAME and GS_PASSWORD environmental variables
  for login."
  (:use [clojure.java.io]
        [clojure.test])
  (:require [clj-genomespace.core :as gs]))

(deftest genomespace-files 
  (let [client (gs/get-client (System/getenv "GS_USERNAME") :password (System/getenv "GS_PASSWORD"))
        test-out-fname "gstest.vcf"]
    (when (.exists (file test-out-fname))
      (.delete (file test-out-fname)))
    (gs/upload client "cdk-test" (str "test/data/" test-out-fname))
    (gs/download client "cdk-test" "gstest.vcf" ".")
    (is true (.exists (file test-out-fname)))
    (println (gs/list-dirs client "."))
    (println (gs/list-files client "cdk-test" "vcf"))))