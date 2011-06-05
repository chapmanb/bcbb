(ns bed-hadoop.core
  (:use [clojure.string :only [split]]
        [cascalog.api])
  (:gen-class))

;; The base query -- define overlap function and query from BED inputs

(deffilterop overlaps [s1 e1 s2 e2]
  "Pass filter if s1,e1 start end pair overlaps s2,e2 pair."
  (or (and (>= s1 s2) (< s1 e2))
      (and (>= e1 s2) (< e1 e2))))

(defn run-bed-overlap [bed-one bed-two]
  "Cascalog query to pull data from two bed files and intersect."
  (?<- (stdout) [?chr ?s1 ?e1]
       (bed-one ?chr ?s1 ?e1)
       (bed-two ?chr ?s2 ?e2) ; Matches chromosome with bed-one
       (overlaps ?s1 ?e1 ?s2 ?e2)))

;; Define example BED intervals for demonstration.
(def test-bed-one
  [
   ["chr1" 20 30]
   ["chr1" 40 50]
   ["chr2" 30 40]
   ["chr2" 40 50]
   ])

(def test-bed-two
  [
   ["chr1" 15 35]
   ["chr2" 45 55]
   ])

(defn memory-bed-overlap []
  (run-bed-overlap test-bed-one test-bed-two))

;; Pull intervals from BED file on HDFS.

(defmapop parse-bed-line [line]
  (take 3 (split line #"\t")))

(defn to-int [x] (Integer/parseInt x))

(defn bed-from-hfs [dir]
  "With a BED file from HDFS, produce chromsome, start, end."
  (let [source (hfs-textline dir)]
    (<- [?chr ?s ?e]
        (source ?line)
        (parse-bed-line ?line :> ?chr ?s-str ?e-str)
        (to-int ?s-str :> ?s)
        (to-int ?e-str :> ?e))))

(defn bed-hfs-overlap [dir1 dir2]
  (run-bed-overlap (bed-from-hfs dir1) (bed-from-hfs dir2)))

(defn -main [dir1 dir2]
  (bed-hfs-overlap dir1 dir2))

