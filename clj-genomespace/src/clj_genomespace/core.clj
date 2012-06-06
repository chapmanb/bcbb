(ns clj-genomespace.core
  (:import [org.genomespace.client GsSession]
           [org.genomespace.client.exceptions AuthorizationException])
  (:use [clojure.java.io]))

;; ## API for accessing GenomeSpace

(defprotocol GsAccess
  "Provide API for accessing GenomeSpace through CDK."
  (upload [this dirname local-file])
  (download [this dirname fname out-dirname])
  (get-user-token [this])
  (list-dirs [this base-dir])
  (list-files [this dirname ftype])
  (logged-in? [this]))

;; ## Helper functions

(defn- gs-user-path
  "Convert relative directory name into full GenomeSpace directory.
  XXX This should use .listPersonalDirectory when available from API."
  ([dm gsuser]
     (gs-user-path dm gsuser nil))
  ([dm gsuser dirname]
     (let [base (str "/" (-> dm .listDefaultDirectory .getDirectory .getName)
                     "/" (.getUsername gsuser))]
       (cond
        (or (nil? dirname) (= dirname ".")) base
        (.startsWith dirname "/") dirname
        :else (str base "/" dirname)))))

(defn- gs-mkdir [dm gsuser dirname]
  (.createDirectory dm
                    (gs-user-path dm gsuser)
                    dirname))

(defn- gs-remote-file
  "Retrieve GenomeSpace reference to remote file."
  [dm gsuser dirname fname]
  (->> (gs-user-path dm gsuser dirname)
       (.list dm)
       .findFiles
       (filter #(= fname (.getName %)))
       first))

(defn- gs-get-dirs
  "Retrieve list of directories relative to the base directory"
  [dm gsuser dirname]
  (let [base (gs-user-path dm gsuser dirname)]
    (map #(str base "/" (.getName %))
         (.findDirectories (.list dm base)))))

(defn- gs-list-files
  "Retrieve files of a specific filetype in a directory."
  [dm gsuser dirname ftype]
  (letfn [(name-and-ftype [gs-file-meta]
            {:name (.getPath gs-file-meta)
             :ftype (.getName (.getDataFormat gs-file-meta))}
            )]
    (let [base (gs-user-path dm gsuser dirname)]
      (->> (.findFiles (.list dm base))
           (map name-and-ftype)
           (filter #(= (:ftype %) ftype))
           (map :name)))))

;; Implementation and factory

(defrecord GsClient [session gsuser dm]
  GsAccess
  (upload [_ dirname local-file]
    (.uploadFile dm (file local-file)
                 (gs-mkdir dm gsuser dirname)))
  (download [_ dirname fname out-dirname]
    (.downloadFile dm (gs-remote-file dm gsuser dirname fname)
                   (file out-dirname fname) false))
  (get-user-token [_]
    (.getToken gsuser))
  (logged-in? [_]
    (.isLoggedIn session))
  (list-dirs [_ base-dir]
    (gs-get-dirs dm gsuser base-dir))
  (list-files [_ dirname ftype]
    (gs-list-files dm gsuser dirname ftype)))

(defmulti get-client
  "Retrieve a GenomeSpace client given username and password or token."
  (fn [_ method _] method))

(defmethod get-client :password
  [user _ passwd]
  (let [session (GsSession.)
        gsuser (try (.login session user passwd)
                    (catch AuthorizationException e nil))]
    (when gsuser
      (GsClient. session gsuser (.getDataManagerClient session)))))

(defmethod get-client :token
  [user _ token]
  (let [session (GsSession. token)
        gsuser (-> (.getUserManagerClient session)
                   (.getUser user))]
    (GsClient. session gsuser (.getDataManagerClient session))))
