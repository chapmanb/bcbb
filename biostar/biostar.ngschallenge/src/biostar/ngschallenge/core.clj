(ns biostar.ngschallenge.core
  (:use [bcbio.variation.variantcontext :only [parse-vcf]]))

(defn find-interest-vrns [fname wanted-genotype]
  "Check for variations of interest with the given genotypes.
  Requires that a call is:
   - Not filtered
   - Has a high or moderate predicted impact
   - Is in the desired genotype (heterozygous or homozygous)"
  (letfn [(is-problem? [vc]
            (and (= 0 (count (:filters vc)))
                 (= (-> vc :genotypes first :type) wanted-genotype)
                 (contains? #{"HIGH" "MODERATE"}
                            (get (:attributes vc) "SNPEFF_IMPACT"))))]
    (filter is-problem? (parse-vcf fname))))

(defn combine-interest-vrns [father-fname mother-fname child-fname]
  "Combine variations of interest from family variant calls."
  (letfn [(get-unique-positions [xs]
            (vec (set (map (juxt :chr :start) xs))))
          (count-by-position [coll x]
            (assoc coll x (+ 1 (get coll x 0))))]
    (reduce count-by-position {}
            (->>
             (map get-unique-positions
                  [(find-interest-vrns father-fname "HET")
                   (find-interest-vrns mother-fname "HET")
                   (find-interest-vrns child-fname "HOM_VAR")])
             flatten
             (partition 2)))))

(defn filter-interest-vrns [pos-map]
  "Print positions of interest present more than once in father/mother/child."
  (println (filter (fn [[_ count]] (> count 1)) pos-map)))

(defn -main [father-fname mother-fname child-fname]
  (->> [father-fname mother-fname child-fname]
       (apply combine-interest-vrns)
       filter-interest-vrns))
