Clojure solution for BioStar code golf:

http://biostar.stackexchange.com/questions/5902/code-golf-finding-orf-and-corresponding-strand-in-a-dna-sequence

With a file of FASTA sequences, it prints the coordinates of the longest
translation product for each sequence. All 6 frames are searched: forward and
reverse in 3 frames.

It uses BioJava to handle FASTA parsing and translation:

http://www.biojava.org

and uses Leiningen to make the dependencies:

https://github.com/technomancy/leiningen#readme

Usage:

    % lein deps
    % lein run :findorf test.fa
    t-1 [0 15 ONE]
    t-2 [1 16 TWO]
    t-3 [2 17 THREE]
    t-r1 [0 15 REVERSED_ONE]
    t-r2 [0 18 REVERSED_TWO]
    t-r3 [0 15 REVERSED_THREE]
    t-1-inner [3 21 ONE]
    EX720612.1 [0 252 ONE]
