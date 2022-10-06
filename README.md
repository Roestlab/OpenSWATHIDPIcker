# IDpicker for OpenSWATH

This is a package that is an implementation
of a protein inference algorithm called IDpicker[1].
It is an algorithm that reads a protein-peptide
mapping as graph, merges protein that map to an identical set
of peptides, separate the graph into components,
then greedily choose the protein that corresponds
to most peptide for each component, until all 
peptides are covered by a chosen protein.

This implementation can only read the output of openSWATH [2]. 
The advantage of this algorithm compare to other ones
for openSWATH is that, this one can handle peptide that 
map to multiple proteins. 

There is only 1 function (just 'main') for this package
with options on
input_file path, q_limit_pep, context, run_id

main(input_file: str, q_limit_pep: str, context: str, run_id: str)

input_file is a path to the osw file

q_limit_pep is a number represent the maximum qvalue the user wants

context is either 'global' or 'run-specific'

run_id is either the id of the run that user want to use,
or if the context is global, it can be anything, it will not be used.

We strongly advice to install PyProphet in a Python virtualenv.

citations

1. Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic parsimony through bipartite graph analysis improves accuracy and transparency. Journal of proteome research, 6(9), 3549–3557. https://doi.org/10.1021/pr070230d
2. Röst, H. L., Rosenberger, G., Navarro, P., Gillet, L., Miladinović, S. M., Schubert, O. T., Wolski, W., Collins, B. C., Malmström, J., Malmström, L., & Aebersold, R. (2014). OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data. Nature biotechnology, 32(3), 219–223. https://doi.org/10.1038/nbt.2841
