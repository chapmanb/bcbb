(import '(org.biojava.bio.seq.io SeqIOTools))
(use '[clojure.java.io])

(defn seq-lengths [seq-iter]
  "Produce a lazy collection of sequence lengths given a BioJava StreamReader"
  (lazy-seq
    (if (.hasNext seq-iter)
      (cons (.length (.nextSequence seq-iter)) (seq-lengths seq-iter)))))

(defn fasta-to-lengths [in-file seq-type]
  "Use BioJava to read a Fasta input file as a StreamReader of sequences"
  (seq-lengths (SeqIOTools/fileToBiojava "fasta" seq-type (reader in-file))))

(defn lazy-avg [coll]
  "Collect the count and number of values in a collection in one pass.

  1 define the function that increments the sum and counts.
  2 reduce the collection, updating the sum and count, and assign to variables
  3 divide the sum by count, keeping it safe in case of 0 division
  "
  (let [f (fn [[s c] val] [(+ s val) (inc c)])
        [sum cnt] (reduce f [0 0] coll)]
    (if (zero? cnt) 0 (/ sum cnt))))

(when *command-line-args*
  (println
    (lazy-avg (apply fasta-to-lengths *command-line-args*))))
