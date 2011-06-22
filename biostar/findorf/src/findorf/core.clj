(ns findorf.core
  (:use [clojure.java.io])
  (:require [clojure.contrib.str-utils2 :as str2])
  (:import [org.biojava3.core.sequence.io FastaReaderHelper]
           [org.biojava3.core.sequence.transcription Frame
            TranscriptionEngine$Builder]))

(defn longest-region [sqn-coll]
  "Retrieve start and end of the longest sequence region in a collection."
  (->>
   (loop [sqn-coll sqn-coll cur-pos 0 coords []]
     (if (empty? sqn-coll) coords
         (let [cur-len (.length (first sqn-coll))
               cur-len-w-stop (+ cur-len (if (empty? (rest sqn-coll)) 0 1))
               next-pos (+ cur-pos cur-len)]
           (recur (rest sqn-coll) next-pos (cons [cur-pos cur-len-w-stop] coords)))))
   (sort-by second >)
   first))

(defn frame-coords [[start end] size frame]
  "Orient coordinates relative to the provided frame."
  (let [map-frame (fn [[start end] frame]
                    (cond
                     (.endsWith frame "TWO") [(+ 1 start) (+ 1 end)]
                     (.endsWith frame "THREE") [(+ 2 start) (+ 2 end)]
                     :else [start end]))
        map-dir (fn [[start end] frame size]
                  (if (.startsWith frame "REVERSED")
                    [(- size end) (- size start)]
                    [start end]))]
    (-> [start end]
        (map-frame frame)
        (map-dir frame size)
        (conj frame))))

(defn tx-sqn-frame [dna frame]
  "Longest transcribed region within a given translation frame."
  (let [tx-engine (-> (TranscriptionEngine$Builder.)
                      (.translateNCodons true)
                      (.trimStop false)
                      (.build))
        dna-coords (fn [coords]
                     (let [dna-start (* 3 (first coords))]
                       [dna-start (+ dna-start (* 3 (second coords)))]))]
    (-> (.getRNASequence dna tx-engine frame)
        (.getProteinSequence tx-engine)
        .toString
        (str2/partition #"\*")
        longest-region
        dna-coords
        (frame-coords (.getLength dna) (str frame)))))

(defn longest-tx-sqn [sqn]
  "Transcribe sequence, returning the longest frame."
  (->> (map #(tx-sqn-frame sqn %) (Frame/getAllFrames))
       (sort-by #(- (second %) (first %)) >)
       first))

(defn fasta-seqs [in-file]
  "Lazy steam of names and DNA sequences from input FASTA file."
  (let [seq-map (FastaReaderHelper/readFastaDNASequence (file in-file) true)]
    (for [key (.keySet seq-map)]
      [key (.get seq-map key)])))

(defn -main [fasta-file]
  (doseq [[name seq] (fasta-seqs fasta-file)]
    (println name (longest-tx-sqn seq))))
