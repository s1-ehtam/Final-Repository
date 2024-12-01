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


