(ns clj-genomespace.core
  (:import [org.genomespace.client GsSession])
  (:use [clojure.java.io]))

;; ## API for accessing GenomeSpace

(defprotocol GsAccess
  "Provide API for accessing GenomeSpace through CDK."
  (gs-upload [this dirname local-file])
  (gs-download [this dirname fname])
  (get-user-token [this]))

;; ## Helper functions

(defn- gs-mkdir [dm dirname]
  (.createDirectory dm
                    (-> dm .listDefaultDirectory .getDirectory)
                    dirname))

(defn- gs-full-path
  "Convert relative directory name into full GenomeSpace directory."
  [dm dirname]
  (str "/users/"
       (-> dm .listDefaultDirectory .getDirectory .getName)
       "/" dirname))

(defn- gs-remote-file
  "Retrieve GenomeSpace reference to remote file."
  [dm dirname fname]
  (->> (gs-full-path dm dirname)
       (.list dm)
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
                   (file fname) false))
  (get-user-token [this]
    (.getToken gsuser)))

(defmulti get-gs-client
  "Retrieve a GenomeSpace client given username and password or token."
  (fn [_ method _] method))

(defmethod get-gs-client :password
  [user _ passwd]
  (let [session (GsSession.)
        gsuser (.login session user passwd)]
    (GsClient. session gsuser (.getDataManagerClient session))))

(defmethod get-gs-client :token
  [user _ token]
  (let [session (GsSession. token)
        gsuser (-> (.getUserManagerClient session)
                   (.getUser user))]
    (GsClient. session gsuser (.getDataManagerClient session))))
