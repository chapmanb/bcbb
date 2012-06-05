(defproject clj-genomespace "0.1-SNAPSHOT"
  :description "Access GenomeSpace data integration platform with simple Clojure API"
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [org.clojars.chapmanb/genomespace-cdk "0.1-SNAPSHOT"]
                 [com.sun.jersey/jersey-client "1.11"]
                 [com.sun.jersey.contribs/jersey-apache-client "1.11"
                  :exclusions [commons-httpclient]]
                 [com.sun.jersey/jersey-json "1.11"]
                 [org.apache.commons/com.springsource.org.apache.commons.httpclient "3.1.0"
                  :exclusions [org.apache.commons/com.springsource.org.apache.commons.codec]]
                 [commons-io "2.0.1"]
                 [org.apache.commons/com.springsource.org.apache.commons.codec "1.4.0"]
                 [commons-lang "2.5"]
                 [org.apache.log4j/com.springsource.org.apache.log4j "1.2.16"]
                 [org.apache.servicemix.bundles/org.apache.servicemix.bundles.jets3t "0.8.1_1"]
                 [org.apache.servicemix.bundles/org.apache.servicemix.bundles.aws-java-sdk
                  "1.3.0_1" :exclusions [org.apache.httpcomponents/httpclient commons-logging
                                         commons-codec org.apache.httpcomponents/httpcore]]
                 [org.apache.httpcomponents/com.springsource.org.apache.httpcomponents.httpclient
                  "4.1.1"]])
