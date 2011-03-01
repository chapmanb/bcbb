(use '[clojure.java.io])
(use '[clojure.contrib.str-utils2 :only (join)])

(defn fasta-lengths [in-file]
  "Generate collection of FASTA record lengths, splitting at '>' delimiters"
  (->> (line-seq (reader in-file))
    (partition-by #(.startsWith ^String % ">"))
    (filter #(not (.startsWith ^String (first %) ">")))
    (map #(join "" %))
    (map #(.length ^String %))))

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
    (lazy-avg (fasta-lengths (first *command-line-args*)))))
