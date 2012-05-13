(ns clj-genomespace.core
  (:import [org.genomespace.client GsSession])
  (:use [clojure.java.io]))

;; ## API for accessing GenomeSpace

(defprotocol GsAccess
  "Provide API for accessing GenomeSpace through CDK."
  (gs-upload [this dirname local-file])
  (gs-download [this dirname fname]))

;; ## Helper functions

(defn- gs-mkdir [dm dirname]
  (.createDirectory dm
                    (-> dm .listDefaultDirectory .getDirectory)
                    dirname))

(defn- gs-remote-file [dm dirname fname]
  (->> (.list dm dirname)
       .findFiles
       (filter #(= fname (.getName %)))
       first))

;; Implementation and factory

(defrecord GsClient [session gsuser dm]
  GsAccess
  (gs-upload [this dirname local-file]
    (.uploadFile dm (file local-file)
                 (gs-mkdir dm dirname)))
  (gs-download [this dirname fname]
    (.downloadFile dm (gs-remote-file dm dirname fname)
                   (file fname) false)))

(defn get-gs-client
  "Retrieve a GenomeSpace client given the login username and password."
  [user passwd]
  (let [session (GsSession.)
        user (.login session user passwd)]
    (GsClient. session user (.getDataManagerClient session))))
