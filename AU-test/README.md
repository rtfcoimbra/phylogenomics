Topology testing 
============================

# TODO

 [] automated AU plotting 
 [] automated species tree estimation
# Pipeline

1. Create GF alignments (previous)
2. Provide alternative trees for AU test.
3. Edit `run-phylo.sh` to collect `.fa` or `.fasta` files and adjust number of threads.
4. Run
```
bash run-phylo.sh dir-with-alignments/
```
5. Continue working with gene trees in `phylo-gftrees.tree` and AU values in `phylo-gftrees.au`.


