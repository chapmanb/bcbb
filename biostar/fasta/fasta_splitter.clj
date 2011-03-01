(import '(org.utgenome.format.fasta FASTAPullParser))

(use '[clojure.java.io]
     '[clojure.contrib.duck-streams :only (with-out-writer)]
     '[clojure.contrib.str-utils2 :only (join)]
     '[clojure.contrib.seq :only (indexed)])

(defn seq-iterator [parser]
  "Lazily retrieve sequence lines one at a time from the FASTA parser."
  (lazy-seq
    (when-let [line (.nextSequenceLine parser)]
      (cons line (seq-iterator parser)))))

(defn partition-seq [parser piece-size]
  "Break up a sequence collection into indexed pieces, pulling lines as needed."
  (->> (seq-iterator parser)
    (mapcat seq)
    (partition-all piece-size)
    indexed))

(defn write-fasta-pieces [in-file piece-size out-file]
  (with-out-writer out-file
    (let [parser (new FASTAPullParser (file in-file))
          piece-size (. Integer parseInt piece-size)]
      (loop [work-parser parser]
        (when-let [cur-id (.nextDescriptionLine work-parser)]
          (doseq [[i s] (partition-seq work-parser piece-size)]
            (println (format ">%s_%s\n%s" cur-id (+ i 1) (join "" s)))
              )
          (recur work-parser))))))

(when *command-line-args*
  (apply write-fasta-pieces *command-line-args*))
