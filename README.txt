For information about purposes and methods of ETHA, see Dadra et al 2016 "Reconstruction of full-length Plasmodium falciparum var exon 1 sequences reveals severe malaria and pregnancy-associated malaria vars in uncomplicated malaria infections in Malian children".

To run ETHA to reconstruct var exon 1 sequences, you will need a set of Illumina reads and a pre-existing whole genome assembly of the same isolate.  You will also need access to three software dependencies:

* Jellyfish (tested with jellyfish-2.0.0beta6.1) http://www.cbcb.umd.edu/software/jellyfish/
* Glimmer (testedw with glimmer-3.02) https://ccb.jhu.edu/software/glimmer/
* MUMmer (tested with version 3.06) http://mummer.sourceforge.net/

You will need to make sure that the executables for each of these packages are available in you path. Edit the primary driver script "run_exon_1" to assign the PATH variable appropriately to include the correct paths on your system.

Running the pipeline consists of three steps:

1) Running Jellyfish on the Illumina reads to get counts of all observed 71mers. See the Jellyfish documentation for instructions for this step.

2) Setting up the working directory with four input files. These should be copied or symlinked to these exact names:
** asm.seq.fa, the whole genome assembly
** 71.mer_counts, the output of step 1
** exon1.all.tail.fa, which lists the tail ends of known exon 1 sequences. A version is included with this code. Augmenting the included version with sequences likely to be similar to those of the target strain may improve sensitivity
**reference.fa, which lists previously known exon1 sequences

3) Running the main driver script:

run_exon_1 $etha $working_directory $lower_kmer_bound $upper_kmer_bound

Here, $etha is the full path of the directory containing the code and this README file, $working_directory is the path of the directory created in step 2, and the kmer bounds are numbers indicating the minimum and maximum numbers of times a 71mer must be seen in the Illumina data to be used. These should be set to reflect the reasonable variation in read depth that characterizes the particular dataset. Note that if you know what value you will be using for the lower bound, you can save storage space by asking Jellyfish to keep only those kmers above that value.

ETHA will run for some hours, putting all of its intermediate and output files in the same working directory. For most purposes, the files you will be most interested in will be these:

* finish/results.deduplicated.fa, the output of ETHA proper, high confidence var sequences
* finish/union.fa, the output of ETHA proper, plus var-like sequences identified in the whole genome assembly which are not accounted for in the ETHA output.

If you run into any difficulty or question that is not addressed here, please email elliott.drabek@gmail.com or jcsilva@som.umaryland.edu
