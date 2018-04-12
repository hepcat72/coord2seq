# coord2seq
Extract subsequences from a larger sequence, given a set of coordinates.

## WHAT IS THIS:

This script takes fasta files and coordinate files (or coordinates supplied on the command line) and extracts the sequences from the fasta file corresponding to the supplied coordinates (inclusive and indexed from a start coordinate of 1).

## INSTALLATION

    perl Makefile.PL
    make
    sudo make install

## USAGE

By example (relates to the example INPUT below):

    coord2seq.pl -i fasta_file -c '301 499' > fasta_file

    coord2seq.pl -i fasta_file -f coord_file -b 1 -e 2 -s '5 1 2' -d 4 -p 5 -i .subseqs

## INPUT

Example (tab delimited):

    1	300	Orf1	+	chrom1
    301	499	Orf2	-	chrom1
    907	634	Orf3		chrom1

## OUTPUT

Fasta format files.
