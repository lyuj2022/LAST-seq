#### retrieve transcript length

```
You can use our "Table Browser" tool (see our Table Browser User's Guide) to accomplish this by following the steps below:

1. Navigate to "Tools > Table Browser" in the top horizontal blue menu bar from our home page.

2. Set your conditions:

Clade: Mammal
Genome: Human
Assembly: Dec. 2013 (GRCh38/hg38)
Group: Genes and Gene Prediction
Track: GENCODE v22
Table: knownCanonical
Region: genome
Identifiers (names/accessions): Click "paste list" or "upload list" to attach your list of genes.
Output format: selected fields from primary and related tables

3. Click "get output" to move to the next step.

4. Under the section, "Select Fields from hg38.knownCanonical," check the following checkboxes:

chrom
chromStart
chromEnd
transcript

5. Under the hg38.kgXref fields, check the checkbox for "Gene Symbol."

6. Click the "get output" button at the top of the page.

From your output, you can now find the transcript length by calculating (chromEnd - chromStart)​.​
```

