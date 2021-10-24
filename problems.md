# Problems

## Task 1

RAP-DB is a public repository for rice genome annotation data. The database was last updated on the 10th of May, 2021. Use the data available from this database to complete this task.
Proteins are made up of 20 amino acids. These amino acids occur in varying degrees of abundance. **Find the rarest amino acid in rice. Build a simple linear regression model to predict the content of the most abundant amino acid using protein length. Use the model to find the outlier protein that has the largest discrepancy between the prediction and the actual number.**

**BONUS:** Find the amino acid that yields the most robust linear model.

## Task 2

Almost all prokaryotic genes have continuous coding sequences which directly translate to protein products. Eukaryotic genes are more obscure as they have intervening sequences called introns between their coding regions.  
Use the GenBank ID U00096.3 to retrieve the complete genome of Escherichia coli strain K-12 substrain MG1655. **Find the discontinuous protein coding genes in this genome.**

Tip: All the data needed to solve this problem are present within the GenBank file and you will not need to look anywhere else.

## Task 3

CRISPR-Cas9 is a powerful and efficient tool for genome editing. To specifically downregulate a gene by knockout using a CRISPR-Cas9 system, a guide RNA sequence complementary to the target is designed. There are some special criteria that guide RNA sequences must fulfill to ensure maximum on-target and minimum off-target effects in the host.
Os03t0752300-01 is the RAP-DB ID of a two-pore potassium channel protein of special interest. **Design an appropriate guide RNA for CRISPR knockout of this gene in rice.**

Tip: There are open-source software developed for estimating on-target and genome-wide off-target scores of CRISPR gRNA sequences.

**BONUS:** Design a plasmid construct for delivery of the Cas9 system along with your designed gRNA using Agrobacteria-mediated transformation.

## Task 4 https://youtu.be/xh_wpWj0AzM?t=2052

Collect RNA-seq raw read count data from supplementary table 2 of [Formentin et al. 2018](https://doi.org/10.3389/fpls.2018.00204). **Find differentially expressed genes under their two experimental conditions at a significance level of 0.01. Create a heatmap to show the fold change in differentially expressed genes.**

**BONUS:** Test all these differentially expressed genes for homology with the two-pore potassium channel and find genes that share at least 60% overall homology.

---

# FASTA Format

## Nucleic Acids

|Code|Meaning|Meaning|
|----|-------|-------|
|A| A |Adenine|
|C| C |Cytosine|
|G| G |Guanine|
|T| T |Thymine|
|U|U |Uracil|
|(i)| i |inosine (non-standard)|
|R|A or G (I) |puRine|
|Y|C, T or U |pYrimidines|
|K|G, T or U |bases which are Ketones|
|M| A or C |bases with aMino groups|
|S| C or G |Strong interaction|
|W| A, T or U |Weak interaction|
|B| not A (i.e. C, G, T or U) |B comes after A|
|D| not C (i.e. A, G, T or U) |D comes after C|
|H| not G (i.e., A, C, T or U) |H comes after G|
|V| neither T nor U (i.e. A, C or G) |V comes after U|
|N| A C G T U |Nucleic acid|
|\-| |gap of indeterminate length |

## Amino Acids

|Code|Amino Acids|
|----|-----------|
|A |Alanine|
|B |Aspartic acid (D) or Asparagine (N)|
|C |Cysteine|
|D |Aspartic acid|
|E |Glutamic acid|
|F |Phenylalanine|
|G |Glycine|
|H |Histidine|
|I |Isoleucine|
|J |Leucine (L) or Isoleucine (I)|
|K |Lysine|
|L |Leucine|
|M |Methionine/Start codon|
|N |Asparagine|
|O |Pyrrolysine (rare)|
|P |Proline|
|Q |Glutamine|
|R |Arginine|
|S |Serine|
|T |Threonine|
|U |Selenocysteine (rare)|
|V |Valine|
|W |Tryptophan|
|Y |Tyrosine|
|Z |Glutamic acid (E) or Glutamine (Q)|
|X |any|
|\*| translation stop|
|\-| gap of indeterminate length|
