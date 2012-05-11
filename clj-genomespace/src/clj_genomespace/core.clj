(ns clj-genomespace.core
  (:import [org.genomespace.client GsSession])
  (:use [clojure.java.io]))

(defn- gs-mkdir [dm-client dirname]
  (.createDirectory dm-client
                    (-> dm-client .listDefaultDirectory .getDirectory)
                    dirname))

(defn- gs-login [user passwd]
  (doto (GsSession.)
    (.login user passwd)))

(defn gs-upload [user passwd dirname local-file]
  (let [session (gs-login user passwd)
        dm-client (.getDataManagerClient session)]
    (.uploadFile dm-client (file local-file)
                 (gs-mkdir dm-client dirname))))

(defn- gs-remote-file [dm dirname fname]
  (->> (.list dm dirname)
       .findFiles
       (filter #(= fname (.getName %)))
       first))

(defn gs-download [user passwd dirname fname]
  (let [dm (.getDataManagerClient (gs-login user passwd))]
    (.downloadFile dm (gs-remote-file dm dirname fname)
                   (file fname) false)))
