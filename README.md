# abui-network
Bayesian phylogenetic network inference in BEAST, using the example of the Timor-Alor-Pantar language Abui

This method was presented at SLE2019 in Leipzig.

## Data Source
The data underlying this study was taken from [LexiRumah 2.2.3](https://github.com/lessersunda/lexirumah-data/commit/67a89d9dc733c068d78288878482ff56eb206e5a), and automatically cognate-coded using LexStat (shipped with [LingPy 2.6.4](https://doi.org/10.5281/zenodo.1544172)) using the pylexirumah interface

    python -m pylexirumah.autocode ~/lexirumah-data/cldf/cldf-metadata.json --soundclass asjp --threshold 0.55 --initial-threshold 0.5
    
The resulting `tsv` file containing alignments contains unescaped quotation marks, which is mitigated by

    sed -e 's/"/”/g' aligned.tsv > beastling/aligned_fixed.tsv
    
The result of that command is version-controlled in the `data/` directory, together with a small (and technically invalid, but working) wrapper to turn it into a CLDF file for use by BEASTling.

## Configuration file generation

The basic configuration is contained in `beastling/spdc.conf`, which describes a tree analysis using a stochastic pseudo-Dollo covarion model of a subset of the data described above. Generate a the corresponding BEAST XML file by

    cd beast
    beastling ../beastling/spdc.conf
    
and apply the necessary changes to run a network inference using the SpeciesNetwork package by running the

    python make_network_config.py
    
script. The grouping of dialects into network taxa and the number of different gene trees to infer are currently hard-coded at the top of that Python script.
This produces a BEAST XML file `network-abui-neighbours-spdc.xml`, which can be run using BEAST, eg. asjp

    mkdir -p ../results
    cd ../results
    beast ../beast/network-abui-neighbours-spdc.xml

## Analysis

The output of the BEAST run consists of several log files. The file `abui-neighbours-spdc.log` contains information about the general progress of the Markov chain, in particular the likelihood values and posterior probabilities of the samples visited. The `genetrees*.trees` files are Nexus Newick files containing the “gene” tree samples. Concepts are associated to specific trees according to the indices logged in `indexes.log`. The file `networks.newick` lists the sampled phylogenetic network structures in eXtended Newick format in a Nexus wrapper.

To normalize the γ annotations of the networks such that the duplicate leaf, which represents the contact edge, is always the node with the minor contribution (γ<0.5), use [my fork of](https://github.com/Anaphory/treesed/) Luke Maurits' [`phyltr`](https://pypi.org/project/phyltr/) tool like

    phyltr cat -n -b 50 networks.newick > normalized.newick
    
If, by unlikely chance, the first network sampled contains no contact edges, you may have to change the burnin (`-b`) parameter or find a workaround. The resulting normalized networks can be inspected in Tim Vaughan's [IcyTree](https://icytree.org). The distance histograms were generated using

    phyltr cat -b 50 networks.newick | phyltr distancehistogram -n --max 6 
    
and the SplitsTree network was generated from the output of
    
    phyltr cat -b 50 networks.newick | phyltr distances -n > distances.nex
    
using SplitsTree 4.

## Adaptation

To vary this analysis, change options in the following places.

1. To change the data, modify the `languages` or `data` entry in `beastling/spdc.conf`. You can also try to vary the other model options defined there, but that might break the network inference script.
2. If you change the set of langugages, you will also have to change their grouping, which is defined in the `species` dictionary in `beast/make_network_config.py`.
3. To change the number of different trees inferred, change the `n_trees` parameter in `beast/make_network_config.py`.

## Requirements

 - BEAST with the following packages:
   - Babel
   - speciesnetwork [with my modifications to prevent NullPointer exceptions, and to provide TreeSelectors]
   in addition to the requirements listed by BEASTling
 - phyltr [with my additions for network support]
 - beastling
 - SplitsTree 4
