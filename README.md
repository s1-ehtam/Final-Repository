# Final-Repository

Welcome to my final repository! Here I will be explaining the analytic methods that I used for my final term paper involving D-2-hydroxy-acid dehydrogenases.

### Background Information
Many gene families have undergone rounds of gene duplication and gene loss in the formation of their lineages from a common ancestral lineage. These modifications may lead to functional diversification among multiple lineages within the gene family where each gene may vary in their biochemical role.

D-2-hydroxy-acid dehydrogenases contain a wide range of oxidoreductases that are responsible for catalyzing the reversible reduction of 2-keto acids to 2-hydroxy acids through the oxidation of electron carriers like NAD+ (Matelska 2018). D-2-hydroxy-acid dehydrogenases possess various substrate specificities resulting in the activation of different cellular mechanisms. The D-2-hydroxy-acid dehydrogenases seem to possess similar biochemical functions and contain similar binding substrates and characteristics. It is likely assumed that there were duplicatin events primarily in the common ancestor or early evolutionary period rather than any recent mass duplication and loss events occurring.

We will be using methods of Bioinformatic Searches, Gene Tree-Species Tree Reconciliation, and Domain Identification to determine the evolutionary history of this gene family.

![image](https://github.com/user-attachments/assets/09d59c09-ca11-488b-9513-16c33fb51f20)


## Find and align D-2-hydroxy-acid dehydrogenases using StartingGene LDHD as the query sequence

First Download the query protein 
```
ncbi-acc-download -F fasta -m protein "NP_919417.1"
```
Peform a blast search usng the query protein
```
blastp -db ../allprotein.fas -query NP_919417.1.fa -outfmt 0 -max_hsps 1 -out LDHD.blastp.typical.out (perform blast search)
```

We will be exploring the blast results. Use the less command to view the output from ```LDHD.blastp.typical.out```
```
less LDHD.blastp.typical.out (look at the output of the blast)
```

Let's create a more detailed and easier-to-process output of the same analysis. The -outfmt flag specifies a particular output format that will be useful for our analysis.
```
blastp -db ../allprotein.fas -query NP_919417.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out LDHD.blastp.detail.out (creates a more detailed and easier to process output of the same analysis)
```

Look at output to find the total human hits using less.
```
less -S LDHD.blastp.detail.out 
```

Use grep to determine the total human hits without having to count by hand.
```
grep -c H.sapiens LDHD.blastp.detail.out (tells the number of hits with the certain species without having to count by hand)
```

Filtering the BLAST output for high-scoring putative homologs. Next we need to choose which putative homologs to include. We only want to include high-scoring matches for our analysis. In addition, it will help to minimize the possibility that we include false homologs. Let's require the e-value to be less than 1e-30.
```
awk '{if ($6< 1e-30)print $1 }' LDHD.blastp.detail.out > LDHD.blastp.detail.filtered.out (require the e-value to be less than 1e-30.)
```

Count the total number of hits in the BLAST results after the filter. We will use the wc command.
```
wc -l LDHD.blastp.detail.filtered.out (count the total number of hits in the BLAST results after the filter.)
```

Use this command below to determine the number of paralogs in each species.
```
grep -o -E "^[A-Z].[a-z]+" LDHD.blastp.detail.filtered.out | sort | uniq -c (tells us the number of paralogs for each of the 11 species)
```

How many homologs does each of the 11 species have?
```
  2 C.carcharias
  2 C.mydas
  2 D.rerio
  2 E.caballus
  2 F.catus
  2 G.aculeatus
  2 G.gallus
  2 H.sapiens
  2 S.salar
  2 S.townsendi
  2 X.laevis
```

## Estimate a maximum-likelihood, midpoint rooted phylogeny for D-2-hydroxy-acid dehydrogenases from Sequence Data
To infer the optimal phylogenetic tree based on a sequence alignment, we will use the software (IQ-TREE)(http://www.iqtree.org/). 

```
sed 's/ /_/g' ~/lab04-$MYGIT/LDHD/LDHD.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas
```
Use IQ-TREE to find the maximum likehood tree estimate. First, it will calculate the optimal amino acid substitution model and amino acid frequencies. Then, it will perform a tree search, estimating branch lengths as it goes.
```
iqtree -s ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas -bb 1000 -nt 2
```
The .iqtree file includes an ASCII graphics (text graphics) version of the tree. You can also display it by reading the .treefile (which is newick formatted) into the nw_display program:
```
nw_display ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas.treefile
```
So, let's look at it unrooted with a graphical display. nw_display can't do this, so... Let's again use our R script . Because there are so many genes, we will make the size of the text labels smaller (0.4), and set the label lengths to 15.
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas. treefile ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas.treefile.pdf 0.4 15
```

Midpoint rooting
We will use a type of rooting called midpoint - we'll hope that the root is halfway along the longest branch on the tree. 
```
gotree reroot midpoint -i ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.fas.treefile -o ~/lab05-$MY GIT/LDHD/LDHD.homologsf.al.mid.treefile
```

Now, we can look at the rooted tree at the command line:
```
nw_order -c n ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile | nw_display -
```

please make a graphic instead
```
nw_order -c n ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity :0' -s > ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile.svg -
```

convert this svg image to a pdf:
```
convert ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/LDHD/LDHD.homolog sf.al.mid.treefile.pdf
```
Branch lengths
The tree shown by default in nw_display (and most other programs) is a phylogram. This means that the lengths of each branch are proportional to the number of substitutions that have accumulated in the sequence along that branch.

If there are very short branch lengths, clades can be hard to visualize on a phylogram. Try switching the view to a cladogram, using the following command.

```
nw_order -c n ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.midCl.treefile.svg -

convert ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.midCl.treefile.pdf
```

The midpoint-rooted ASCII tree for D-2-hydroxy-acid dehydrogenases.
```
                                                /--+ D.rerio ldhd probable Dlactate dehydrogenase mitochondrial                       
                                                |                                                                                     
                                                //42-+ G.aculeatus ldhd probable Dlactate dehydrogenase mitochondrial                 
                                                \+ 69                                                                                 
                                                |\-+ S.salar ldhd LOW QUALITY PROTEIN probable Dlactate dehydrogenase                 
 /----------------------------------------------+ 100                                                                                 
 |                                              |   /---+ C.carcharias ldhd probable Dlactate dehydrogenase mitochondrial             
 |                                              |   |                                                                                 
 |                                              |   | /--+ X.laevis ldhd.L probable Dlactate dehydrogenase mitochondrial              
 |                                              \---+ 90                                                                              
 |                                                  | |  /-+ S.townsendi LDHD probable Dlactate dehydrogenase mitochondrial           
 |                                                  \-+ 93                                                                            
 |                                                    | /+/58+ C.mydas LDHD probable Dlactate dehydrogenase mitochondrial isofor      
 |                                                    | |\+ 62                                                                        
 |                                                    \-+ 100------+ G.gallus LDHD probable Dlactate dehydrogenase mitochondrial isofo
=+                                                      |                                                                             
 |                                                      |   /+ H.sapiens LDHD probable Dlactate dehydrogenase mitochondrial isof      
 |                                                      |   |                                                                         
 |                                                      \---/+100caballus LDHD probable Dlactate dehydrogenase mitochondrial          
 |                                                          \ 95                                                                      
 |                                                          \-+ F.catus LDHD probable Dlactate dehydrogenase mitochondrial isofor     
 |                                                                                                                                    
 |                                                   /+ S.salar LOC106575380 D2hydroxyglutarate dehydrogenase mitochondri             
 |                                                   |                                                                                
 \---------------------------------------------------+/100G.aculeatus d2hgdh D2hydroxyglutarate dehydrogenase mitochondrial           
                                                     ||                                                                               
                                                     \+/49+ D.rerio d2hgdh D2hydroxyglutarate dehydrogenase mitochondrial pre         
                                                      ||                                                                              
                                                      \+ 54/-----+ C.carcharias d2hgdh D2hydroxyglutarate dehydrogenase mitochondria  
                                                       |   |                                                                          
                                                       |   /-------+ X.laevis d2hgdh.S D2hydroxyglutarate dehydrogenase mitochondrial 
                                                       \---| 89                                                                       
                                                           | /--+ S.townsendi D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial   
                                                           \ 80                                                                       
                                                           |/+/93 C.mydas D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial iso   
                                                           ||\+ 85                                                                    
                                                           \+ 71-+ G.gallus D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial     
                                                            |                                                                         
                                                            |   /-+ F.catus D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial iso 
                                                            |   |                                                                     
                                                            \---+/100caballus D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial   
                                                                \+ 54                                                                 
                                                                 \-+ H.sapiens D2HGDH D2hydroxyglutarate dehydrogenase mitochondrial i
```

## Reconciliation for the D-2-hydroxy-acid dehydrogenases
```
cp ~/lab05-$MYGIT/LDHD/LDHD.homologsf.al.mid.treefile ~/l
ab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile
```
Reconcile the gene and species tree using Notung; Use the following command to perform the reconciliation:

```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile --reconcile --
speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/LDHD/
```
We can see the node names that notung used/assigned to these internal nodes with this command. But for our species tree, all internal nodes have names, so this doesn't provide any additional information.
```
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/LDHD/LDHD.homolog
s.al.mid.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```
Generate a RecPhyloXML object and view the gene-within-species tree via thirdkind
Use the following command:

```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile.rec.ntg --include.species
\argument  -g : /home/bio312-user/lab06-s1-ehtam/LDHD/LDHD.homologs.al.mid.treefile.rec.ntg
--include.species
```
To create a gene-reconciliation-within species tree reconciliation graphic, use thirdkind:
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile.rec.svg
```
Convert to a pdf for easy viewing:
```
convert  -density 150 ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile.rec.svg ~/lab06-$MYGIT/LDHD/LDHD.homologs.al.mid.treefile.rec.pdf

```
The thirdkind and notung visualizations should be in your repository, along with all intermediate files.

1. What is the cost of the reconciled tree?

Cost = 27

2. How many copies were in the common ancestor of Gnathostomes?

6 copies

3. Describe the events leading to the common ancestor of Tetrapods.

5 speciation events to 2 speciation events and 4 losses


