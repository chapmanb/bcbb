# GenomeSpace with Clojure

This is a simple API to access [GenomeSpace][1] from Clojure using the Java
CDK. This allows upload and download of files to GenomeSpace. GenomeSpace
makes these files available to Galaxy, GenePattern and other tools.

The library is available from [Clojars][2] for inclusion in [Leiningen][3]
managed projects.

[1]: http://www.genomespace.org/
[2]: https://clojars.org/clj-genomespace
[3]: http://leiningen.org/

## Usage

Download Clojure libraries and the GenomeSpace CDK and start a REPL:

    $ lein deps
    $ lein repl

Usage, from the REPL:

    user> (require '[clj-genomespace.core :as gs])
    user> (def client (gs/get-client "chapmanb" :password "password"))
    user> (gs/upload client "cdk-test" "/path/to/yourfile.vcf")
    user> (gs/download client "cdk-test" "yourfile.vcf" ".")
    user> (gs/list-dirs client ".")
    user> (gs/list-files client "cdk-test" "vcf")
    

## License

The code is freely available under the [MIT license][l1].

[l1]: http://www.opensource.org/licenses/mit-license.html
