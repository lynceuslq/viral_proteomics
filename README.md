# viral_proteomics
pipelines to build sequential and functional profiles for RNA viruses

![rdrp-discovery-04-06-2021](https://user-images.githubusercontent.com/55744039/125015736-02047800-e068-11eb-8b32-96563288a7fa.png)

Here is an example of using the pipelines for RNA-dependent-RNA-polymerase discovery,
[rdrp-discovery-04-06-2021-liqian.pdf](https://github.com/lynceuslq/viral_proteomics/files/6788370/rdrp-discovery-04-06-2021-liqian.pdf)


Dependencies and databses 
1. Blast suite
2. Diamond (v2.9 and above) and non-redundant protein database for Diamond blastp (nr.dmnd) # here I use nr as databse for protein database, but namely you can use others such as Uniprot, as long as generated by Diamond
3. CD-Hit suite
4. InterproScan & its database

Iput files
1. a fasta file containing seed protein sequences (to define sequential similarity)
2. a list of conserved domains from CDD,separated by "\n" (to define functional similarity)
3. a list of viral families you want to look at, separated by "\n"

