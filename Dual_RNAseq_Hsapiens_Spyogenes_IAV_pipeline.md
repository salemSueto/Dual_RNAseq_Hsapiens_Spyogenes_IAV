#Dual RNA-seq analysis of Influenza A and Streptococcus pyogenes infection in epithelial human cells

The project has twelve RNA-seq samples. These samples are divided into four groups each with three replicates. The four groups are:

* Homo sapiens epithelial cells (Control)
* Homo sapiens epithelial cells infected by S. pyogenes AP1 (Inf_AP1)
* Homo sapiens epithelial cells infected by S. pyogenes NZ131 (Inf_NZ131)
* Homo sapiens epithelial cells infected by Influenza A virus (IAV)

##Table of Contents

1.	[Directories creation](#section1.)
2.	[Preparation](#section2.)
 	1.	[Tools](#section2.1.)
 	2.	[Sample metadata](#section2.2.)
 	3.	[Reference](#section2.3.)
 	4.	[R Package](#section2.4.)
 	5.	[STAR: genome reference index generation](#section2.5.)
 	6.	[Raw FastQ files](#section2.6.)
 	7.	[Trimmomatic adapters](#section2.7.)
 	8.	[Scotty](#section2.8.)
 	9.	[DrInsight](#section2.9.)
 	10.	[Pathogen gene length table](#section2.10.)
 	11.	[BLAST database index](#section2.11.)
3.	[Pre-processing analysis](#section3.)
4. [Mapping & Quantification analysis](#section4.)
5. [Gene Count Table Generation](#section5.)
6. [Scotty: Gene Count Table Quality Control](#section6.)
7. [Homo sapiens DEGs identification & enrichment analysis](#section7.)
8. [Drug repurposing analysis](#section8.)
	1.	[DGIdb database](#section8.1.)
	2. [L1000CDS2 tool](#section8.2.)
	3. [DrInsight tool](#section8.3.)
	4. [Drug Selection](#section8.4.)
9. [S. pyogenes TPM Analysis](#section9.)
10. [Streptococcus pyogenes strains Homology Analysis](#section10.)


---
[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtEdWFsIFJOQS1zZXFdIC0tPkJbUHJlcGFyYXRpb25dXG4gICAgQiAtLT5DW1JhdyBkYXRhXVxuICAgIEMgLS0-RFtQcmUtcHJvY2Vzc2luZ11cbiAgICBEIC0tPkVbTWFwcGluZyAmIFF1YW50aWZpY2F0aW9uXVxuICAgIEUgLS0-RltBbmFseXNpc11cbiAgICBGIC0tPkdbU2NvdHR5XVxuICAgIEYgLS0-SFtIb21vIHNhcGllbnNdXG4gICAgSCAtLT5JW0RFR3MgaWRlbnRpZmljYXRpb25dXG4gICAgSSAtLT5KW0VucmljaG1lbnRdXG4gICAgSSAtLT5LW0RydWcgcmVwdXJwb3NpbmddXG4gICAgRiAtLT5MW1BhdGhvZ2VuIFRQTV1cbiAgICBMIC0tPk1bUy4gcHlvZ2VuZXMgaG9tb2xvZ3ldXG5cbiAgICAgICAgICAgICIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0IiwidGhlbWVWYXJpYWJsZXMiOnsiYmFja2dyb3VuZCI6IndoaXRlIiwicHJpbWFyeUNvbG9yIjoiI0VDRUNGRiIsInNlY29uZGFyeUNvbG9yIjoiI2ZmZmZkZSIsInRlcnRpYXJ5Q29sb3IiOiJoc2woODAsIDEwMCUsIDk2LjI3NDUwOTgwMzklKSIsInByaW1hcnlCb3JkZXJDb2xvciI6ImhzbCgyNDAsIDYwJSwgODYuMjc0NTA5ODAzOSUpIiwic2Vjb25kYXJ5Qm9yZGVyQ29sb3IiOiJoc2woNjAsIDYwJSwgODMuNTI5NDExNzY0NyUpIiwidGVydGlhcnlCb3JkZXJDb2xvciI6ImhzbCg4MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJwcmltYXJ5VGV4dENvbG9yIjoiIzEzMTMwMCIsInNlY29uZGFyeVRleHRDb2xvciI6IiMwMDAwMjEiLCJ0ZXJ0aWFyeVRleHRDb2xvciI6InJnYig5LjUwMDAwMDAwMDEsIDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxKSIsImxpbmVDb2xvciI6IiMzMzMzMzMiLCJ0ZXh0Q29sb3IiOiIjMzMzIiwibWFpbkJrZyI6IiNFQ0VDRkYiLCJzZWNvbmRCa2ciOiIjZmZmZmRlIiwiYm9yZGVyMSI6IiM5MzcwREIiLCJib3JkZXIyIjoiI2FhYWEzMyIsImFycm93aGVhZENvbG9yIjoiIzMzMzMzMyIsImZvbnRGYW1pbHkiOiJcInRyZWJ1Y2hldCBtc1wiLCB2ZXJkYW5hLCBhcmlhbCIsImZvbnRTaXplIjoiMTZweCIsImxhYmVsQmFja2dyb3VuZCI6IiNlOGU4ZTgiLCJub2RlQmtnIjoiI0VDRUNGRiIsIm5vZGVCb3JkZXIiOiIjOTM3MERCIiwiY2x1c3RlckJrZyI6IiNmZmZmZGUiLCJjbHVzdGVyQm9yZGVyIjoiI2FhYWEzMyIsImRlZmF1bHRMaW5rQ29sb3IiOiIjMzMzMzMzIiwidGl0bGVDb2xvciI6IiMzMzMiLCJlZGdlTGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsImFjdG9yQm9yZGVyIjoiaHNsKDI1OS42MjYxNjgyMjQzLCA1OS43NzY1MzYzMTI4JSwgODcuOTAxOTYwNzg0MyUpIiwiYWN0b3JCa2ciOiIjRUNFQ0ZGIiwiYWN0b3JUZXh0Q29sb3IiOiJibGFjayIsImFjdG9yTGluZUNvbG9yIjoiZ3JleSIsInNpZ25hbENvbG9yIjoiIzMzMyIsInNpZ25hbFRleHRDb2xvciI6IiMzMzMiLCJsYWJlbEJveEJrZ0NvbG9yIjoiI0VDRUNGRiIsImxhYmVsQm94Qm9yZGVyQ29sb3IiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJsYWJlbFRleHRDb2xvciI6ImJsYWNrIiwibG9vcFRleHRDb2xvciI6ImJsYWNrIiwibm90ZUJvcmRlckNvbG9yIjoiI2FhYWEzMyIsIm5vdGVCa2dDb2xvciI6IiNmZmY1YWQiLCJub3RlVGV4dENvbG9yIjoiYmxhY2siLCJhY3RpdmF0aW9uQm9yZGVyQ29sb3IiOiIjNjY2IiwiYWN0aXZhdGlvbkJrZ0NvbG9yIjoiI2Y0ZjRmNCIsInNlcXVlbmNlTnVtYmVyQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvciI6InJnYmEoMTAyLCAxMDIsIDI1NSwgMC40OSkiLCJhbHRTZWN0aW9uQmtnQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvcjIiOiIjZmZmNDAwIiwidGFza0JvcmRlckNvbG9yIjoiIzUzNGZiYyIsInRhc2tCa2dDb2xvciI6IiM4YTkwZGQiLCJ0YXNrVGV4dExpZ2h0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0RGFya0NvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dE91dHNpZGVDb2xvciI6ImJsYWNrIiwidGFza1RleHRDbGlja2FibGVDb2xvciI6IiMwMDMxNjMiLCJhY3RpdmVUYXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwiYWN0aXZlVGFza0JrZ0NvbG9yIjoiI2JmYzdmZiIsImdyaWRDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQmtnQ29sb3IiOiJsaWdodGdyZXkiLCJkb25lVGFza0JvcmRlckNvbG9yIjoiZ3JleSIsImNyaXRCb3JkZXJDb2xvciI6IiNmZjg4ODgiLCJjcml0QmtnQ29sb3IiOiJyZWQiLCJ0b2RheUxpbmVDb2xvciI6InJlZCIsImxhYmVsQ29sb3IiOiJibGFjayIsImVycm9yQmtnQ29sb3IiOiIjNTUyMjIyIiwiZXJyb3JUZXh0Q29sb3IiOiIjNTUyMjIyIiwiY2xhc3NUZXh0IjoiIzEzMTMwMCIsImZpbGxUeXBlMCI6IiNFQ0VDRkYiLCJmaWxsVHlwZTEiOiIjZmZmZmRlIiwiZmlsbFR5cGUyIjoiaHNsKDMwNCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGUzIjoiaHNsKDEyNCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIiwiZmlsbFR5cGU0IjoiaHNsKDE3NiwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU1IjoiaHNsKC00LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTYiOiJoc2woOCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU3IjoiaHNsKDE4OCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIn19LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtEdWFsIFJOQS1zZXFdIC0tPkJbUHJlcGFyYXRpb25dXG4gICAgQiAtLT5DW1JhdyBkYXRhXVxuICAgIEMgLS0-RFtQcmUtcHJvY2Vzc2luZ11cbiAgICBEIC0tPkVbTWFwcGluZyAmIFF1YW50aWZpY2F0aW9uXVxuICAgIEUgLS0-RltBbmFseXNpc11cbiAgICBGIC0tPkdbU2NvdHR5XVxuICAgIEYgLS0-SFtIb21vIHNhcGllbnNdXG4gICAgSCAtLT5JW0RFR3MgaWRlbnRpZmljYXRpb25dXG4gICAgSSAtLT5KW0VucmljaG1lbnRdXG4gICAgSSAtLT5LW0RydWcgcmVwdXJwb3NpbmddXG4gICAgRiAtLT5MW1BhdGhvZ2VuIFRQTV1cbiAgICBMIC0tPk1bUy4gcHlvZ2VuZXMgaG9tb2xvZ3ldXG5cbiAgICAgICAgICAgICIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0IiwidGhlbWVWYXJpYWJsZXMiOnsiYmFja2dyb3VuZCI6IndoaXRlIiwicHJpbWFyeUNvbG9yIjoiI0VDRUNGRiIsInNlY29uZGFyeUNvbG9yIjoiI2ZmZmZkZSIsInRlcnRpYXJ5Q29sb3IiOiJoc2woODAsIDEwMCUsIDk2LjI3NDUwOTgwMzklKSIsInByaW1hcnlCb3JkZXJDb2xvciI6ImhzbCgyNDAsIDYwJSwgODYuMjc0NTA5ODAzOSUpIiwic2Vjb25kYXJ5Qm9yZGVyQ29sb3IiOiJoc2woNjAsIDYwJSwgODMuNTI5NDExNzY0NyUpIiwidGVydGlhcnlCb3JkZXJDb2xvciI6ImhzbCg4MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJwcmltYXJ5VGV4dENvbG9yIjoiIzEzMTMwMCIsInNlY29uZGFyeVRleHRDb2xvciI6IiMwMDAwMjEiLCJ0ZXJ0aWFyeVRleHRDb2xvciI6InJnYig5LjUwMDAwMDAwMDEsIDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxKSIsImxpbmVDb2xvciI6IiMzMzMzMzMiLCJ0ZXh0Q29sb3IiOiIjMzMzIiwibWFpbkJrZyI6IiNFQ0VDRkYiLCJzZWNvbmRCa2ciOiIjZmZmZmRlIiwiYm9yZGVyMSI6IiM5MzcwREIiLCJib3JkZXIyIjoiI2FhYWEzMyIsImFycm93aGVhZENvbG9yIjoiIzMzMzMzMyIsImZvbnRGYW1pbHkiOiJcInRyZWJ1Y2hldCBtc1wiLCB2ZXJkYW5hLCBhcmlhbCIsImZvbnRTaXplIjoiMTZweCIsImxhYmVsQmFja2dyb3VuZCI6IiNlOGU4ZTgiLCJub2RlQmtnIjoiI0VDRUNGRiIsIm5vZGVCb3JkZXIiOiIjOTM3MERCIiwiY2x1c3RlckJrZyI6IiNmZmZmZGUiLCJjbHVzdGVyQm9yZGVyIjoiI2FhYWEzMyIsImRlZmF1bHRMaW5rQ29sb3IiOiIjMzMzMzMzIiwidGl0bGVDb2xvciI6IiMzMzMiLCJlZGdlTGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsImFjdG9yQm9yZGVyIjoiaHNsKDI1OS42MjYxNjgyMjQzLCA1OS43NzY1MzYzMTI4JSwgODcuOTAxOTYwNzg0MyUpIiwiYWN0b3JCa2ciOiIjRUNFQ0ZGIiwiYWN0b3JUZXh0Q29sb3IiOiJibGFjayIsImFjdG9yTGluZUNvbG9yIjoiZ3JleSIsInNpZ25hbENvbG9yIjoiIzMzMyIsInNpZ25hbFRleHRDb2xvciI6IiMzMzMiLCJsYWJlbEJveEJrZ0NvbG9yIjoiI0VDRUNGRiIsImxhYmVsQm94Qm9yZGVyQ29sb3IiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJsYWJlbFRleHRDb2xvciI6ImJsYWNrIiwibG9vcFRleHRDb2xvciI6ImJsYWNrIiwibm90ZUJvcmRlckNvbG9yIjoiI2FhYWEzMyIsIm5vdGVCa2dDb2xvciI6IiNmZmY1YWQiLCJub3RlVGV4dENvbG9yIjoiYmxhY2siLCJhY3RpdmF0aW9uQm9yZGVyQ29sb3IiOiIjNjY2IiwiYWN0aXZhdGlvbkJrZ0NvbG9yIjoiI2Y0ZjRmNCIsInNlcXVlbmNlTnVtYmVyQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvciI6InJnYmEoMTAyLCAxMDIsIDI1NSwgMC40OSkiLCJhbHRTZWN0aW9uQmtnQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvcjIiOiIjZmZmNDAwIiwidGFza0JvcmRlckNvbG9yIjoiIzUzNGZiYyIsInRhc2tCa2dDb2xvciI6IiM4YTkwZGQiLCJ0YXNrVGV4dExpZ2h0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0RGFya0NvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dE91dHNpZGVDb2xvciI6ImJsYWNrIiwidGFza1RleHRDbGlja2FibGVDb2xvciI6IiMwMDMxNjMiLCJhY3RpdmVUYXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwiYWN0aXZlVGFza0JrZ0NvbG9yIjoiI2JmYzdmZiIsImdyaWRDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQmtnQ29sb3IiOiJsaWdodGdyZXkiLCJkb25lVGFza0JvcmRlckNvbG9yIjoiZ3JleSIsImNyaXRCb3JkZXJDb2xvciI6IiNmZjg4ODgiLCJjcml0QmtnQ29sb3IiOiJyZWQiLCJ0b2RheUxpbmVDb2xvciI6InJlZCIsImxhYmVsQ29sb3IiOiJibGFjayIsImVycm9yQmtnQ29sb3IiOiIjNTUyMjIyIiwiZXJyb3JUZXh0Q29sb3IiOiIjNTUyMjIyIiwiY2xhc3NUZXh0IjoiIzEzMTMwMCIsImZpbGxUeXBlMCI6IiNFQ0VDRkYiLCJmaWxsVHlwZTEiOiIjZmZmZmRlIiwiZmlsbFR5cGUyIjoiaHNsKDMwNCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGUzIjoiaHNsKDEyNCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIiwiZmlsbFR5cGU0IjoiaHNsKDE3NiwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU1IjoiaHNsKC00LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTYiOiJoc2woOCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU3IjoiaHNsKDE4OCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIn19LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)
---

<h2 id="section1">1.	Directories creation</h2>

```bash
$	mkdir dual_RNAseq_Hs_Sp_IAV
$	cd  dual_RNAseq_Hs_Sp_IAV
$	mkdir Tool
$	mkdir Tool/Adapters
$	mkdir Reference
$	mkdir Reference/Hsapiens
$	mkdir Reference/Hsapiens/star_index
$	mkdir Reference/Spyogenes_AP1
$	mkdir Reference/Spyogenes_AP1/star_index
$	mkdir Reference/Spyogenes_NZ131
$	mkdir Reference/Spyogenes_NZ131/star_index
$	mkdir Reference/IAV
$	mkdir Reference/IAV/star_index
$	mkdir 0_Raw_Data
$	mkdir 0_Raw_Data/QC
$	mkdir 1_Trimmed_Data
$	mkdir 1_Trimmed_Data/QC
$	mkdir 2_Hsapiens_Mapped
$	mkdir 2_Hsapiens_UnMapped
$	mkdir 2_Spyogenes_AP1_Mapped
$	mkdir 2_Spyogenes_NZ131_Mapped
$	mkdir 2_IAV_Mapped
$	mkdir 3_Quantification
$	mkdir 3_Quantification/Hsapiens
$	mkdir 3_Quantification/Spyogenes_AP1
$	mkdir 3_Quantification/Spyogenes_NZ131
$	mkdir 3_Quantification/IAV
$	mkdir 4_DEG_Analysis
$	mkdir 4_DEG_Analysis/Hsapiens
$	mkdir 5_Drug_DrInsight_Analysis/Hsapiens/
$	mkdir 6_TPM_Analysis
$	mkdir 6_TPM_Analysis/Spyogenes_AP1
$	mkdir 6_TPM_Analysis/Spyogenes_NZ131
$	mkdir 6_TPM_Analysis/IAV
$	mkdir 7_Homology_Analysis
```

<h2 id="section2">2.	Preparation </h2>

<h3 id="section2.1.">2.1.	Tools </h3>

The Conda package manager was used to install all the command-line tools, except for BLAST. Conda installation can be found [here](https://bioconda.github.io/user/install.html). In case of conflict between the user's system and Conda, the command-line tools can be manually downloaded and installed. R's tools were installed and uploaded directly from R.

| Tool         | Version | Type         | Date Accessed | Manual                                                                                                | Installation                                | Source                                                                          | 
|--------------|---------|--------------|---------------|-------------------------------------------------------------------------------------------------------|---------------------------------------------|---------------------------------------------------------------------------------| 
| FastQC       | 0.11.8  | Commad-Line  | NA            | [Manual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)                                          | conda install -c bioconda fastqc            | [Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)           | 
| MultiQC      | 1.2     | Commad-Line  | NA            | [Tutorial](https://multiqc.info/docs/)                                                                | conda install -c bioconda multiqc           | [Github](https://github.com/ewels/MultiQC)                                      | 
| Trimmomatic  | 0.39    | Commad-Line  | NA            | [Manual](http://www.usadellab.org/cms/?page=trimmomatic)                                              | conda install -c bioconda trimmomatic       | [Website](http://www.usadellab.org/cms/?page=trimmomatic)                       | 
|  STAR        | 2.7.1a  | Commad-Line  | NA            | [Manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) | conda install -c bioconda star              | [Github](https://github.com/alexdobin/STAR)                                     | 
|  HTSeq-count | 0.6     | Commad-Line  | NA            | [Manual](https://htseq.readthedocs.io/en/release_0.11.1/count.html)                                   | conda install -c bioconda htseq             | [Website](https://pypi.org/project/HTSeq/)                                      | 
|  Scotty      | NA      | Matlab       | Apr-20        | [Help](http://scotty.genetics.utah.edu/help.html)                                                     | NA                                          | [Github](https://github.com/mbusby/Scotty)                                      |                                                                         | 
|  REVIGO      | NA      | Website      | May-20        | [Tool Website](http://revigo.irb.hr/index.jsp?error=expired)                                          | NA                                          | NA                                                                              | 
|  BLAST       | 2.9.0   | Commad-Line  | NA            | [Manual](https://www.ncbi.nlm.nih.gov/books/NBK279668/)                                               | NA                                          | [Latest Version](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | 
|  Pfam        | NA      | Website      | Jun-20        | [Website] (https://www.ebi.ac.uk/Tools/hmmer/)                                                        | NA                                          | NA                                                                              | 
|  DGIdb       | NA      | Website      | May-20        | [Website](https://www.dgidb.org/)                                                                     | NA                                          | NA                                                                              |                                                                              | 
|  L1000CDS2   | NA      | Website      | May-20        | [Website](https://maayanlab.cloud/L1000CDS2/#/index)                                                  | NA                                          | NA                                                                              | 
|  Cytoscape   | 3.7.2   | Desktop Tool | NA            | [Tutorial](https://github.com/cytoscape/cytoscape-tutorials/wiki)                                     | [Page](https://cytoscape.org/download.html) | NA                                                                              | 
|  STITCH      | v11.0   | App          | NA            | NA                                                                                                    | Cytoscape App Installation                  | NA                                                                              | 

<h3 id="section2.2.">2.2.	Sample metadata </h3>

Create a metadata text file for the samples called "samples\_metadata.txt" (columns called ID and Group) in the dual\_RNAseq\_Hs\_Sp\_IAV directory.

```
ID	Group
RNA_01	Control
RNA_02	Control
RNA_03	Control
RNA_04	Inf_AP1
RNA_05	Inf_AP1
RNA_06 Inf_AP1
RNA_07	Inf_NZ131
RNA_08	Inf_NZ131
RNA_09	Inf_NZ131
RNA_10	Inf_IAV
RNA_11	Inf_IAV
RNA_12	Inf_IAV
```

<h3 id="section2.3.">2.3.	Reference </h3>

The following references were downloaded and saved in the path assigned.

| Organism                      | Genome                  | Annotation              | Proteome | GenBank         | Download                                                                                | Path                      | 
|-------------------------------|-------------------------|-------------------------|---------------|-----------------|-----------------------------------------------------------------------------------------|---------------------------| 
|  Homo sapiens                 | GRCh38                  | GRCh38.95               | Not Needed    | GRCh38          | [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)                              | Reference/Hsapiens        | 
|  Streptococcus pyogenes AP1   | ASM99376v1              | ASM99376v1.37           | ASM99376v1    | GCA_000993765.1 | [Ensembl](https://bacteria.ensembl.org/Streptococcus_pyogenes_gca_000993765/Info/Index) | Reference/Spyogenes_AP1   | 
|  Streptococcus pyogenes NZ131 | ASM1812v1               | ASM1812v1.37            | ASM99376v1    | GCA_000018125.1 | [Ensembl](https://bacteria.ensembl.org/Streptococcus_pyogenes_nz131/Info/Index)         | Reference/Spyogenes_NZ131 | 
|  Influenza A virus            | ViralMultiSegProj15622  | ViralMultiSegProj15622  | Not Needed    | NA              | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000865085.1/)                          | Reference/IAV             | 


<h3 id="section2.4.">2.4.	R Package </h3>


| Program    | Version  | Purpose                                            | Installation                       | Citation                                                                                                                                                                                                                            | 
|------------|----------|----------------------------------------------------|------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
| Biostrings | 2.56.0   | Manipulation of biological strings.                | BiocManager::install("Biostrings") | H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.56.0.                                                                                        | 
| DESeq2     | 1.28.1   | Differential expression                            | BiocManager::install("DESeq2")     | Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)                                                                                  | 
| gprofiler2 | 2_0.2.0  | Gene list functional enrichment analysis           | install.packages("gprofiler2")     | Liis Kolberg and Uku Raudvere (2020). gprofiler2: Interface to the 'g:Profiler' Toolset. R package version 0.2.0. https://CRAN.R-project.org/package=gprofiler2                                                                     | 
| gridExtra  | 2.3      | Arrange multiple grid-based plots on a page        | install.packages("gridExtra")      | Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3. https://CRAN.R-project.org/package=gridExtra                                                                                 | 
| rvest      | 0.3.6    | Download, then manipulate, HTML and XML            | install.packages("rvest")          | Hadley Wickham (2020). rvest: Easily Harvest (Scrape) Web Pages. R package version 0.3.6. https://CRAN.R-project.org/package=rvest                                                                                                  | 
| scales     | 1.1.1    | Graphical scales map data to aesthetics            | install.packages("scales")         | Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales                                                                                | 
| stringi    | 1.5.3    | Language processing tools                          | install.packages("stringi")        | Gagolewski M. and others (2020). R package stringi: Character string processing facilities. http://www.gagolewski.com/software/stringi/.                                                                                            | 
| tidytext   | 0.2.6    | Text mining tasks and plot generation.             | install.packages("tidytext")       | Silge J, Robinson D (2016). "tidytext: Text Mining and Analysis Using Tidy Data Principles in R." _JOSS_, *1*(3). doi: 10.21105/joss.00037 (URL:https://doi.org/10.21105/joss.00037), <URL: http://dx.doi.org/10.21105/joss.00037>. | 
| tidyverse  | 1.3.0    | Collection of R packages designed for data science | install.packages("tidyverse")      | Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686                                                                                                 | 
| DrInsight  | 0.1.2    | Drug Connectivity Mapping                          | install.packages("DrInsight")      | Jinyan Chan and Jinghua Gu (2020). DrInsight: Drug Repurposing: Integration and Systematic Investigation of Genomic High Throughput Data. R package version 0.1.2. https://CRAN.R-project.org/package=DrInsight                     | 

<h3 id="section2.5.">2.5.	STAR: genome reference index generation </h3>

STAR genome references indexes were created for all four organisms based on their reference genome and their annotation.

```bash
#Homo sapiens
$	star --runMode genomeGenerate --runThreadN 10 --genomeDir Reference/Hsapiens/star_index/ --genomeFastaFiles Reference/Hsapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa --genomeSAsparseD 4 --genomeSAindexNbases 10 --sjdbGTFfile Reference/Hsapiens/Homo_sapiens.GRCh38.95.gtf --sjdbOverhang 74

#S.pyogenes AP1
$	star --runMode genomeGenerate --runThreadN 7 --genomeDir Reference/Spyogenes_AP1/star_index/ --genomeFastaFiles Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.cdna.all.fa --genomeSAindexNbases 9 --sjdbGTFfile Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.37.gtf --sjdbOverhang 74

#S.pyogenes NZ131
$	star --runMode genomeGenerate --runThreadN 7 --genomeDir Reference/Spyogenes_NZ131/star_index/ --genomeFastaFiles Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.cdna.all.fa --genomeSAindexNbases 9 --sjdbGTFfile Reference/Spyogenes_AP1/Streptococcus_pyogenes_nz131.ASM1812v1.37.gtf --sjdbOverhang 74

#Influenza A virus
$	star --runMode genomeGenerate --runThreadN 7 --genomeDir Reference/IAV/star_index/ --genomeFastaFiles Reference/IAV/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna --genomeSAindexNbases 5 --limitGenomeGenerateRAM 14000000000 --limitSjdbInsertNsj 250000 --sjdbGTFfile Reference/IAV/GCF_000865085.1_ViralMultiSegProj15622_genomic.gtf
```

<h3 id="section2.6.">2.6.	Raw FastQ files </h3>

Save all the raw FastQ files in the directory called "0\_Raw\_Data". Unzip the eventual compacted raw files (gz extension) by using the most apt commad showcased below.

```bash
$	tar xvzf 0_Raw_Data/*.gz
$	gunzip 0_Raw_Data/*.gz
$	gzip -d 0_Raw_Data/*.gz
```

<h3 id="section2.7.">2.7.	Trimmomatic Adapters </h3>

[Download](https://github.com/timflutre/trimmomatic/tree/master/adapters) the adapter sequence files from the Tool/Adapters directory.

<h3 id="section2.8.">2.8.	Scotty </h3>

The Scotty files were written using an older version of Matlab.  At some point Matlab depricated one of the functions used which will cause running errors on newer versions (Matlab 2013+). It is imperative to change the function "RandStream.setDefaultStream(s)" into "RandStream.setGlobalStream(s)" in the following files: getProbabilitiesFromFit (line 23) and getProbSequencedMonteCarloByReads (line 7).

Furthermore, in order to save the tables used to create the Power Optimization and Poisson Noise plots, the following changes were made to the file "getOptimizationCharts.m" after line 80 which showcase the line: "readsSX=readsS/10^6;".  

```
csvwrite('powerOptimization.csv', pwrsCalc)
csvwrite('poissonNoise.csv', pwrBiases)
csvwrite('seqDepth.csv', readsSX)
csvwrite('nReps.csv', nReps)
```

<h3 id="section2.9.">2.9.	DrInsight </h3>

Download the drug data matrix from the Download section at the [Broad Institute](https://portals.broadinstitute.org/cmap/), or through this download [link](ftp://ftp.broadinstitute.org/pub/cmap/rankMatrix.txt.zip). Renamed it into "rankMatrix.txt" and saved it inside the "dual\_RNAseq\_Hs\_Sp\_IAV" directory.

<h3 id="section2.10.">2.10.	Pathogen gene length table </h3>

Create a text file called "geneLength.txt" in the "dual\_RNAseq\_Hs\_Sp\_IAV" directory which has the pathogen's gene and its length shown in two columns. The first column is called "ID" which has the gene ID and the second called "Length" which has the length of the gene. Run the R script called "pathogen\_geneLength\_creation.R".

```R
#!/pathogen_geneLength_creation.R
#The script takes, as inputs, the annotation files of the three pathogens and creates a text file called "geneLength.txt".
#The geneLength.txt file has the following columns: ID (geneID), Length (gene's length), Reference (name of the annotation), Type (Gene/Transcript).

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
annoPath <- c("Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.37.gtf",
             "Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.gtf",
             "Reference/IAV/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff")
##########################

geneLength_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(geneLength_df) <- c("ID", "Length", "Reference", "Type")

for (iAnnoPath in annoPath) 
{
  anno <- rtracklayer::import(iAnnoPath) 
  
  if (grepl( "gtf", iAnnoPath, fixed = TRUE) == TRUE)
  {
    gtf_df <- as.data.frame(anno) %>% filter(type=="exon" & gene_id != "") %>% mutate (Length = end - start + 1)
    gtf_df$ID <- gsub("gene:*","", gtf_df$gene_id)
    
    filter_gtf_df <- gtf_df %>% select(ID, Length)
    filter_gtf_df$Reference <- tail(str_split(iAnnoPath, pattern = "/")[[1]], n=1)
    filter_gtf_df$Type <- "Gene"
    geneLength_df <- rbind(geneLength_df, filter_gtf_df)
  }
  
  if (grepl( "gff", iAnnoPath, fixed = TRUE) == TRUE)
  {
    gff_df <- as.data.frame(anno) %>% filter(type=="gene") %>% mutate (Length = end - start + 1)
    
    filter_gff_df <- gff_df %>% select(ID, Length)
    filter_gff_df$Reference <- tail(str_split(iAnnoPath, pattern = "/")[[1]], n=1)
    filter_gff_df$Type <- "Gene"
    geneLength_df <- rbind(geneLength_df, filter_gff_df)
  }
}

write.table (geneLength_df, file="geneLength.txt", sep="\t", row.names=FALSE)
```

<h3 id="section2.11.">2.11.	BLAST database index </h3>

BLAST database indexes were created for the proteome of both S. pyogenes strains (AP1 and NZ131).

```bash
#S. pyogenes AP1
$	./Tool/ncbi-blast-2.11.0+/bin/makeblastdb -in Reference/Spyogenes_AP1/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.fa -parse_seqids -blastdb_version 5 -title "Spyogenes_NZ131" -dbtype prot

#S. pyogenes NZ131 
$	./Tool/ncbi-blast-2.11.0+/bin/makeblastdb -in Reference/Spyogenes_NZ131/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.fa -parse_seqids -blastdb_version 5 -title "Spyogenes_AP1" -dbtype prot
``` 

<h2 id="section3.">3.	Pre-processing analysis </h2>

First, the raw data quality control were accessed with FastQC and their result files were merged together with MultiQC. Then, the raw data were trimmed by the Trimmomatic tool. Last, the trimmed data quality control were accessed with FastQC and their result files were merged together with MultiQC.

```bash
#Raw Data Quality Control
$	fastqc 0_Raw_Data/*fastq -o 0_Raw_Data/QC
$	multiqc 0_Raw_Data/QC -o 0_Raw_Data/QC

#Raw Data Trimming
$	for f in 0_Raw_Data/*fastq; do trimmomatic SE $f 1_Trimmed_Data/${f%%fastq}"trim.fastq" ILLUMINACLIP:Tool/Adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50; done

#Trimmed Data Quality Control
$	fastqc 1_Trimmed_Data/*fastq -o 1_Trimmed_Data/QC
$	multiqc 1_Trimmed_Data/QC -o 1_Trimmed_Data/QC
```

<h2 id="section4.">4.	Mapping & Quantification analysis </h2>

The STAR tool was used for the mapping and the integrated HTSeq-count tool was used to measured the gene count for all three types of sequencing (un-stranded, stranded-yes, and stranded-reverse).

```bash
#Homo sapiens
$	for f in 1_Trimmed_Data/*fastq; do star --genomeDir Reference/Hsapiens/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --outFileNamePrefix 2_Hsapiens_Mapped/${f%%fastq}; done
$	mv 2_Hsapiens_Mapped/*UnMapped* 2_Hsapiens_UnMapped/
$	mv 2_Hsapiens_Mapped/*ReadsPerGene* 3_Quantification/Hsapiens

#Streptococcus pyogenes AP1
$	for f in 2_Hsapiens_UnMapped/*; do star --genomeDir Reference/Spyogenes_AP1/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_Spyogenes_AP1_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_Spyogenes_AP1_Mapped/*ReadsPerGene* 3_Quantification/Spyogenes_AP1

#Streptococcus pyogenes NZ131
$	for f in 2_Hsapiens_UnMapped/*mate1; do star --genomeDir Reference/Spyogenes_NZ131/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_Spyogenes_NZ131_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_Spyogenes_NZ131_Mapped/*ReadsPerGene* 3_Quantification/Spyogenes_NZ131

#Influenza A virus
$	for f in 2_Hsapiens_UnMapped/*mate1; do star --genomeDir Reference/IAV/star_index --readFilesIn $f --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix 2_IAV_Mapped/${f%%_Unmapped.out.mate1}; done
$	mv 2_IAV_Mapped/*ReadsPerGene* 3_Quantification/IAV
```

<h2 id="section5.">5.	Gene Count Table Generation</h2>

The following R script merges all gene count files obtained from STAR and HTSeq-count saved in the sub-directories of 3\_Quantification directory. The script outputs two different type of files (count table and summary table) for the unstranded, strand, and reverse strand. Furthermore, it also outputs barplots of the three quantification types. The output file called "htseq\_all\_sample\_strand\_rev\_count.txt" was selected for Scotty, DEGs, and enrichment analyses. Run the "merge\_geneCountTables\_from\_starHTSeq.R".


```R
#!/merge_geneCountTables_from_starHTSeq.R
#The script merges the files generated after the mapping-quantification step with STAR and HTSeq-count with the flag "--quantMode".
#Each input file has the information for one single sample under study and it is divided in two parts.
#The upper section has general information of the quantification result.
#The below section has four columns without headers.
#The first column has the gene's ID. Then, counts for unstranded reads. Then, counts for stranded reads. Finally, counts for reversed stranded reads.
#The script outputs two type of files (count table and summary table) for the unstranded, strand, and reverse strand.
#And, it also outputs one pdf file of barplots for the three quantification types.
#These output files are created for each directory these input files are found.

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gridExtra)	#Provides the ability to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
pathFiles <- c("3_Quantification/Hsapiens", "3_Quantification/Spyogenes_AP1", "3_Quantification/Spyogenes_NZ131") #List of directories where input files reside
patternFiles <- "ReadsPerGene" #Pattern for the right selection of files


###### Loop through each Directory to find the count files and merge them
for (iDr in pathFiles) 
{
  #Get filenames
  files <- list.files(path=iDr, pattern=patternFiles)
  
  # Save all files
  list_of_files <- list()
  for (i in files) 
  {list_of_files[[i]] <- read.table(paste(iDr, "/", i, sep = ""), header=FALSE)}
  
  summary_plot_df <- as.data.frame (matrix(nrow = 0, ncol = 4))
  colnames(summary_plot_df) <- c("ID", "variable", "value", "MapType")
  selectColID <- c("unstranded", "strand_yes", "strand_rev")
  
  #Merge the files & separate the count table from the summary
  for (selectCol in c(2,3,4)) #2: unstranded - #3: 1st read strand aligned with RNA (htseq-count -s yes) - #4: 1st read strand aligned with RNA (htseq-count -s reverse)
  {
    # Select the column for each files
    merge_htseq_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]])-4, ncol = length(list_of_files)))
    summary_htseq_df <- as.data.frame (matrix(data = NA, nrow = 4, ncol = length(list_of_files)))
    
    colnames(merge_htseq_df) <- names(list_of_files)
    colnames(summary_htseq_df) <- names(list_of_files)
    
    rownames(merge_htseq_df) <- list_of_files[[1]]$V1 [-(1:4)]
    rownames(summary_htseq_df) <- list_of_files[[1]]$V1 [1:4]
    
    for (i in 1:length(list_of_files)) 
    {
      sFile <- list_of_files[[i]] %>% select (contains(as.character(selectCol)))
      summary_htseq_df[,i] <- sFile[1:4, 1]
      merge_htseq_df [,i] <- tail(sFile, -4)
    }
    
    merge_htseq_colsums <- colSums(merge_htseq_df)
    merge_htseq_df <- merge_htseq_df %>% rownames_to_column(var = "ID")
    
    summary_htseq_df <- rbind (summary_htseq_df, N_mapped_Quantified = merge_htseq_colsums) %>% rownames_to_column(var = "ID")
    summary_htseq_df$MapType <- selectColID[selectCol-1]
    
    temp_summary <- reshape2::melt(summary_htseq_df,id.vars=c("ID"))
    temp_summary$MapType <- selectColID[selectCol-1]
    summary_plot_df <- rbind(summary_plot_df, temp_summary)
    
    fileName <- paste("htseq_all_sample_", selectColID[selectCol-1], sep = "")
    write.table (merge_htseq_df, paste (iDr, "/", fileName, "_count.txt", sep = ""), sep="\t", row.names=FALSE)
    write.table (summary_htseq_df,paste(iDr, "/", fileName, "_summary.txt", sep = ""), sep="\t", row.names=FALSE)
    
    ####PLOT SUMMARY INFO #####
    summary_plot_df <- summary_plot_df %>% filter(value!=MapType)
    plots_list <- list()
    
    for (i in unique(summary_plot_df$ID))
    {
      plots_list[[paste("Variable:", i)]] <- summary_plot_df %>%
        filter(ID==i) %>%
        ggplot(aes(x=variable, y=value, fill=MapType)) + 
        geom_bar(position = "dodge", stat="identity") +
        labs(title=paste("Variable:", i)) + 
        theme(axis.text.x = element_text(angle=65, hjust=1.0))
    }
    
    pdf(paste(iDr, "/", "htseq_summary_plots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(plots_list))) {grid.arrange (plots_list[[i]])}
    dev.off()
  }
}
```

<h2 id="section6.">6.	Scotty: Gene Count Table Quality Control</h2>

The Scotty tool was used for two main objectives:

1. Accertain the overall quality of the count table by using the rarefaction curves of the samples.
2. Predict the ability of the count table to identify all differential expressed genes.

Scotty analysis was carried out for the following three group comparisons:

1.  S. pyogenes AP1 infected samples against H. sapiens healthy control samples.
2.  S. pyogenes NZ131 infected samples against H. sapiens healthy control samples.
3.  Influenza A virus infected samples against H. sapiens healthy control samples.
 
The following command was used in the Matlab software as explained in the ReadMe.pdf file. The input file called "scotty.txt" has the following structure.

```
#GeneID	#Group1Sample	#Group1Sample	#Group1Sample	#Group2Sample	#Group2Sample	#Group2Sample
```

Run the following command in the MatLab software.

```Matlab
scottyEstimate('./scotty.txt', '3', '3', '1990634', '2.0', '0.05', '10', '0', '0', '0', 'Inf', '10', '10000000', '100000000', '50', '50', '87', 'outDirectory');
```


<h2 id="section7.">7.	Homo sapiens DEGs identification & enrichment analysis</h2>

The DEGs identification was carried out by DESeq2 at P-value < 0.0005. Whereas, the enrichment was carried out with g:Profiler. The analysis was conducted only for Homo sapiens species because it is the only organism in the project which posses a control group. Make sure that the ID in the samples\_metadata.txt file has the same name in the "htseq\_all\_sample\_strand\_rev\_count.txt". Run the deg\_enrich\_analysis.R script.

```R
#!/deg_enrich_analysis.R
#The script takes several inputs from the user.
#The path to the count table files saved in the vector object called "path_RawCount".
#The raw table files need to have the first column as ID followed by the samples with their gene count.
#The "path_samplesMetadata" object takes the path to the metadata about the samples.
#The samples's metadata have two columns: ID column with the sample name and Group column with the sample's group membership.
#The "groupComparisons" object vector has the list of group comparisosns to perform saved as Group1 - Group2.
#The DESeq2 input section has the value to identify DEGs. 
#The "cutOff" remove genes which row sum is less than its value. The "padj_FDR_value" is the adjusted p-value.
#The "gProfiler_supported_org" need the scientific name of the organism to perform enrichment analysis.
#The "REVIGO's Input format" has all the inputs need to perform at the REVIGO website.

#The following output files are created:
#1) DESeq2 output file (only DEGs) for each group comparison (format: Group1_vs_Group2_DESeq2_cutOff_adjPvalue.csv).
#2) One single file for the enrichment analysis for all group comparisons (format: OrganismName_gProfiler_enrich.csv).
#3) REVIGO results for each GO class (format: OrganismName_revigo_GO/class.csv).
#4) Tables obtained after the cluster generation from REVIGO's GO semantic space (format: OrganismName_revigoCluster_GO/class.csv) 
#5) Table used to generate the heatmap for the three GO class (format: OrganismName_revigoClusterHeatmap_GO/class.csv)
#6) Table used to generate the heatmap for KEGG (format: OrganismName_heatmapKEGG.csv)
#7) A PDF file with heatmap plots (format: OrganismName_heatmapPlots.pdf).

#Load Library
library(tidyverse)	#Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2)	#Differential gene expression analysis based on the negative binomial distribution
library(gprofiler2)	#Gene list functional enrichment analysis and namespace conversion
library(rvest)	#Wrappers around the 'xml2' and 'httr' packages to make it easy to download, then manipulate, HTML and XML.
library(tidytext)	#Text mining tasks, and plot generation are easier and more effective with tools already in wide use: dplyr, broom, tidyr, and ggplot2.
library(gridExtra)	#Provides the ability to arrange multiple grid-based plots on a page, and draw tables.
library(stringi)	#A multitude of language processing tools: pattern searching, concatenation, sorting, padding, wrapping, and many more.
library(scales)	#Graphical scales map data to aesthetics, and provide methods for automatically determining breaks and labels for axes and legends.

rm(list=ls(all=TRUE))

##### User's INPUT #####
path_RawCount <- c("3_Quantification/Hsapiens/htseq_all_sample_strand_rev_count.txt")   #Path to the count table files
path_samplesMetadata <- "samples_metadata.txt"                                          #Path to the metadata file (#ID #Group)
groupComparisons <- c("Inf_AP1 - Control", "Inf_NZ131 - Control", "Inf_IAV - Control")  #Specify the group comparisons as: Group1 - Group2

##### DESeq2 input
cutOff <- 0
padj_FDR_value <- 0.0005

##### Enrichment Analysis -> Supported Organims for g:Profiler tool (e.i. "Homo sapiens") -> (https://biit.cs.ut.ee/gprofiler/page/organism-list)
gProfiler_supported_org <- c("Homo sapiens") #Scientific Name 

# REVIGO's Input format
baseurl <- "http://revigo.irb.hr/"  #REVIGO website link (Do not change)
cutoff <- "0.40" #Allowed values: "0.90" "0.70" "0.50" "0.40" 
isPValue <- "yes" #Allowed values: "yes"  "no"
whatIsBetter <- "higher" #Allowed values: "higher" "lower" "absolute" "abs_log"
goSizes <- "0" #("whole UniProt")
measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"
######


###### Loop through each count table to identify the DEGs ######
for (iCT in path_RawCount) 
{
  deg_list <- list() #List of the DEGs identified
  
  # Get the count table & samples' metadata files
  rawCount_df <- read.table(iCT, header=TRUE)
  samplesMetadata_df <- read.table(path_samplesMetadata, header=TRUE)
  
  # Loop through the different group comparisons
  for (iGC in groupComparisons) 
  {
    #Check that there are two saved groups
    if (length (str_split(iGC, pattern = " - ")[[1]]) == 2)
    {
      groups <- str_split(iGC, pattern = " - ")[[1]]
      
      #Group1
      gr1ColNames <- samplesMetadata_df %>% filter(Group==groups[1])
      group1 <- grep(paste(gr1ColNames$ID, collapse="|"), colnames(rawCount_df))
      group1_df <- rawCount_df[, group1]
      
      #Group2
      gr1ColNames <- samplesMetadata_df %>% filter(Group==groups[2])
      group2 <- grep(paste(gr1ColNames$ID, collapse="|"), colnames(rawCount_df))
      group2_df <- rawCount_df[, group2]
      
      #Replicates
      replicate <- c(ncol(group1_df), ncol(group2_df))
      
      #Check that each group has more than one replicate
      if (all(replicate > 1)) 
      {
        cat("The Differential Expressed Analisis will be started!\n")
        
        #Fix the Count Table -> nowTable
        nowTable <- cbind(group1_df, group2_df)
        nowTable <- nowTable %>% mutate_if (is.character,as.numeric)
        rownames(nowTable) <- rawCount_df[,1]
        
        #Loop through each instance of cutoff value
        for (c in cutOff) 
        {
          print(paste(iGC, " * ", "CutOff: ", c, sep = ""))
          
          colSelect <- append(c(1:replicate[1]), c((replicate[1]+1):(sum(replicate))))
          colSelect <- colnames(nowTable)[colSelect]
          
          compaDF <- nowTable %>% 
            dplyr::select(colSelect) %>% 
            rownames_to_column() %>%
            mutate(Sum=rowSums(dplyr::select(., -starts_with("rowname")))) %>% 
            filter(Sum >= c) %>% 
            dplyr::select(-Sum) %>%
            column_to_rownames()
          
          # Group Factors
          sampleGroups <- append (rep(1, replicate[1]), rep(2, replicate[2]))
          myfactors_D_E <- factor(sampleGroups)             #DESeq2 factors
          
          #------- DESeq2 ------#
          matrixCountTable <- round (as.matrix(compaDF))
          coldata <- data.frame(row.names=colnames(matrixCountTable), myfactors_D_E)
          dds <- DESeqDataSetFromMatrix(countData=matrixCountTable, colData=coldata, design=~myfactors_D_E)
          res <- as.data.frame(results(DESeq(dds)))
          
          #Adjusted P-values
          for (pv in padj_FDR_value) 
          {
            #DESeq2
            deseq2_deg <- res %>% filter(padj < pv)
            #deseq2_deg$GroupComparison <- iGC
            write.csv (deseq2_deg, file = paste("4_DEG_Analysis/Hsapiens/", groups[1], "_vs_", groups[2], "_DESeq2_", c, "_", pv, ".csv", sep = ""))
            deg_list[[paste (groups[1], groups[2], "DESeq2", c, pv, sep = "_")]] <- rownames(deseq2_deg)
          }
        }
      } else {cat(paste(gc, " * ", "Replicates -> ", "The samples' replicates need to be above 1\n"))}
    } else {cat(paste(gc, "-> ", "The present comparison needs exactly two groups for the comparison 1\n", sep = ""))}
  }
}


###### g:Profiler Supported Organism Enrichment Analysis ######
enrichPlot <- list()

for (iSpecies in gProfiler_supported_org) 
{
  ##### Enrichment Analysis -> g:Profiler #####
  gProfiler_Res <- NULL
  gPr_iSpecies <- paste(substr (tolower(iSpecies), 1, 1), str_split (tolower(iSpecies), pattern = " ")[[1]][2], sep = "")
  
  tryCatch (
    {
      gProfiler_Res <- gost(deg_list, organism = gPr_iSpecies, ordered_query = FALSE,
                            multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                            measure_underrepresentation = FALSE, evcodes = FALSE,
                            user_threshold = 0.05, correction_method = "gSCS",
                            domain_scope = "annotated", custom_bg = NULL,
                            numeric_ns = "", sources = NULL)
    },
    error = function(e) 
    {
      message(e)
      return(gProfiler_Res <- NULL)
    }
  )
  
  if (is.null(gProfiler_Res) == FALSE) 
  {
    enrich_go_kegg <- gProfiler_Res$result[,c(1,3,9,10,11)] %>% filter(source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
    write.csv(enrich_go_kegg, file = paste("4_DEG_Analysis/Hsapiens/", gPr_iSpecies, "_gProfiler_enrich.csv", sep = ""))
    
    ##### REVIGO Analysis #####
    revigo_session <- html_session(baseurl)
    revigo_form <- html_form(revigo_session)[[1]] 
    
    # Loop through the GO terms
    for (go in c("GO:BP", "GO:CC", "GO:MF")) 
    {
      # Number of Cluster
      nCluster <- 40
      
      # Gene Ontology
      filterGO <- enrich_go_kegg %>% filter (source==go)
      goList <- paste(unique(filterGO$term_id), collapse="\n")
      
      filled_form <- set_values(revigo_form, 'goList'= goList, 
                                'cutoff'= cutoff, 'isPValue'= isPValue, 
                                'whatIsBetter'= whatIsBetter, 'goSizes'= goSizes, 'measure'= measure)
      
      # Get the Revigo Summary Table & Save it as a CSV file
      result_page <- submit_form(revigo_session, filled_form, submit='startRevigo')
      GORevigo <- result_page %>% follow_link("Export results to text table (CSV)")
      httr::content(GORevigo$response, as="text") %>% write_lines("revigoSummaryOnline.csv") 
      
      # Get the REVIGO summary result from the saved file
      revigoSummary <- read.csv2("revigoSummaryOnline.csv", header=TRUE, sep=",")
      revigoSummary$frequency <- gsub('%', '', revigoSummary$frequency)
      revigoSummary$frequency <- as.numeric( as.character(revigoSummary$frequency) );
      
      #Get only the GO terms represented in the "semantic space"
      uniqueGO <- revigoSummary %>% filter(plot_X != "null")
      iNum <- c(4,5,6,7,8) # Change the column into numeric: plot_X - plot_Y - plot_size - uniqueness - dispensability
      uniqueGO [ , iNum] <- apply(uniqueGO [ , iNum], 2, function(x) as.numeric(as.character(x)))
      
      # Add "HeadGO" Column which shows which GO REVIGO term has been put as the main GO for that group of terms
      revigoSummary$HeadGO <- NA
      headGO <- "NA"
      
      for (i in 1:nrow(revigoSummary)) 
      {
        if(revigoSummary$plot_X[i]!="null") {headGO <- revigoSummary$term_ID[i]}
        revigoSummary$HeadGO[i] <- headGO
      }
      write.csv(revigoSummary, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigo_%s.csv", gPr_iSpecies, go))   #Save the REVIGO result
      
      
      # GENE ONTOLOGY REVIGO CLUSTER ANALYSIS
      #Find the clusters present in the REVIGO results based on their sematic space
      clusterDF <- uniqueGO %>% remove_rownames %>% column_to_rownames(var="term_ID")
      if (nCluster >= nrow(clusterDF)) {nCluster <- nrow(clusterDF)-1}
      
      set.seed(100)
      clusters <- kmeans(clusterDF[,c(3,4)], nCluster)
      clusterDF$Cluster <- clusters$cluster
      
      # Group elements from the same cluster
      revigoCluster <- data.frame(matrix(ncol=10, nrow=nCluster))
      colnames(revigoCluster) <- c("Cluster", "PlotX", "PlotY", "RevigoGOs", "RevigoRep", "RevigoDescription", "AllGOs", "AllWords", "AllDescription", "FinalDescription")
      
      for (i in sort(unique(clusterDF$Cluster)))
      {
        filterGO_cluster <- clusterDF %>% filter(Cluster==i)
        
        # REVIGO GOs
        revigoCluster$Cluster[i] <- i
        revigoCluster$PlotX[i] <- mean(filterGO_cluster$plot_X)
        revigoCluster$PlotY[i] <- mean(filterGO_cluster$plot_Y)
        revigoCluster$RevigoGOs[i] <- paste (rownames(filterGO_cluster), collapse=" ")
        revigoCluster$RevigoRep[i] <- filterGO_cluster[which.max(filterGO_cluster$representative),][1,1] #Term with the max value of "representative of the cluster
        
        # REVIGO Cluster Description -> Concagenate the HEADs' term description of the cluster - Split - Delete generic terms (i.e. "of", "a", "an")
        wordOcc <- paste(filterGO_cluster$description, collapse=" ") %>% str_split(pattern = " ") 
        wordOcc <- wordOcc[[1]]
        indexWords <- !(wordOcc %in% stop_words$word)
        wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
        wordOcc <- as.vector(wordOcc$Var1[1:5])
        revigoCluster$RevigoDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
        
        #All GOs Columns
        filterAllGO <- revigoSummary %>% filter(HeadGO %in% rownames(filterGO_cluster))
        
        revigoCluster$AllGOs[i] <- paste (filterAllGO$term_ID, collapse=" ")
        revigoCluster$AllWords[i] <- paste(filterAllGO$description, collapse=" ")
        
        wordOcc <- paste(filterAllGO$description, collapse=" ") %>% str_split(pattern = " ") 
        wordOcc <- wordOcc[[1]]
        indexWords <- !(wordOcc %in% stop_words$word)
        wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
        wordOcc <- as.vector(wordOcc$Var1[1:5])
        revigoCluster$AllDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
        
        #Final Description of the Cluster
        revigoCluster$FinalDescription[i] <- paste(revigoCluster$RevigoRep[i], revigoCluster$AllDescription[i], nrow(filterAllGO), sep = " * ")
      }
      write.csv(revigoCluster, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigoCluster_%s.csv", gPr_iSpecies, go))   #Save the REVIGO_Cluster result
      
      
      #Get the lowest pvalue of the group comparison's GO terms and the cluster combination
      heatmapGO <- data.frame(matrix(ncol = length (unique(filterGO$query)), nrow = nrow(revigoCluster)))
      colnames(heatmapGO) <- unique(filterGO$query)
      rownames(heatmapGO) <- revigoCluster$FinalDescription
      
      for (n in 1:nrow(revigoCluster)) 
      {
        allGO <- str_split(revigoCluster$AllGOs[n], " ") [[1]] #Get all GOs 
        filterGOGroup <- filterGO %>% filter(term_id %in% allGO)  #Get all the rows which have those GOs
        
        uniqueGroup <- unique (filterGOGroup$query)
        for (k in uniqueGroup) 
        {
          pvaluesList <- filterGOGroup %>% filter(query == k)
          heatmapGO[n, grep(k, colnames(heatmapGO))] <- min(pvaluesList$p_value)
        }
      }
      write.csv(heatmapGO, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigoClusterHeatmap_%s.csv", gPr_iSpecies, go))   #Save the REVIGO_Cluster_Heatmap result
      
      #Create the Heatmap plot
      GroupName <- c()
      for (i in 1:nrow(heatmapGO)) 
      {
        groupPresence <- colnames(heatmapGO)[sapply(heatmapGO[i,], function(x) any(!is.na(x)))]
        groupPresence <- paste(groupPresence, collapse = " * ")
        GroupName <- append(GroupName, groupPresence)
      }
      
      heatmapGO$GroupName <- GroupName
      heatmapGO$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapGO))
      heatmapGO <- heatmapGO %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
      
      heatmapGO <- tibble::rownames_to_column(heatmapGO, "Description")
      heatmapGO_melt <- reshape2::melt(heatmapGO,id.vars=c("Description"))
      
      heatmapGO_melt$Description <- factor(heatmapGO_melt$Description, levels = heatmapGO$Description)
      
      plotName <- stri_replace_all_fixed (paste(iSpecies, go), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
      enrichPlot[[plotName]] <- heatmapGO_melt %>% 
        ggplot(aes(x=variable, y=factor(Description), fill=value)) +
        geom_raster(aes(fill = value)) +
        scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
        geom_text(size=3, aes(label=scientific(value, digits = 3))) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
        theme(axis.text.y = element_text(size=10, face="bold")) +
        coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
    }
    
    ##### KEGG Heatmap creation #####
    filterKEGG <- filter(enrich_go_kegg, source=="KEGG")
    heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
    colnames(heatmapKEGG) <- unique(filterKEGG$query)
    rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
    
    for (i in 1:nrow(filterKEGG)) 
    {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
    write.csv(heatmapKEGG, file = paste("4_DEG_Analysis/Hsapiens/", gPr_iSpecies, "_heatmapKEGG.csv", sep = ""))   #Save the gProfiler_Heatmap result
    
    #Create the Heatmap plot
    GroupName <- c()
    for (i in 1:nrow(heatmapKEGG)) 
    {
      groupPresence <- colnames(heatmapKEGG)[sapply(heatmapKEGG[i,], function(x) any(!is.na(x)))]
      groupPresence <- paste(groupPresence, collapse = " * ")
      GroupName <- append(GroupName, groupPresence)
    }
    
    heatmapKEGG$GroupName <- GroupName
    heatmapKEGG$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapKEGG))
    heatmapKEGG <- heatmapKEGG %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
    
    heatmapKEGG <- tibble::rownames_to_column(heatmapKEGG, "Description")
    heatmapKEGG_melt <- reshape2::melt(heatmapKEGG,id.vars=c("Description"))
    heatmapKEGG_melt$Description <- factor(heatmapKEGG_melt$Description, levels = heatmapKEGG$Description)
    heatValues <- heatmapKEGG_melt$value [!heatmapKEGG_melt$value %in% c(NA)]
    
    plotName <- stri_replace_all_fixed (paste(iSpecies, "KEGG"), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
    enrichPlot[[plotName]] <- heatmapKEGG_melt %>% 
      ggplot(aes(x=variable, y=Description, fill=value)) +
      geom_raster(aes(fill = value)) +
      scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
      geom_text(size=3, aes(label=scientific(value, digits = 3))) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
      theme(axis.text.y = element_text(size=10, face="bold")) +
      coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
    
    ##### Delete the "revigoSummaryOnline.csv" file
    if (file.exists("revigoSummaryOnline.csv")) {file.remove("revigoSummaryOnline.csv")}
    
    ##### Save Heatmap plots #####
    pdf(paste("4_DEG_Analysis/Hsapiens/", iSpecies, "_heatmapPlots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(enrichPlot))) {grid.arrange (enrichPlot[[i]])}
    dev.off()
  }
}
```

<h2 id="section8.">8.	Drug repurposing analysis</h2>

<h3 id="section8.1.">8.1.	DGIdb database </h3>

* Access the [DGIdb] (https://www.dgidb.org/) website. 
* Click on the "Search Drug-Gene Interactions" button.
* Insert the list of human differential expressed genes in the window under the "Gene" tab.
* Click on the "Find Drug-Gene Interactions" button
* Click on the "Download as TSV" button to download the results.

<h3 id="section8.2.">8.2.	L1000CDS2 tool </h3>

The tool can be access at the following [website](https://maayanlab.cloud/L1000CDS2/#/index).
Paste the list of differential expressed genes with their respectiful Log2 Fold Change value from the DESeq2 analysis separated by a comma.

```
GeneID,Log2FoldChange
```

Set the configuration window settings to:
* "reverse" which search for drugs capable of reversing the input data  
* "latest" which set the latest updated database for the drug search

Click on the "Search" button and save the table by clicking on the download button on the upper-right side of the screen.

<h3 id="section8.3.">8.3.	DrInsight tool </h3>

DrInsight tool was applied on each DEGs analysis results obtained with DESeq2 tool. Run the "drInsight\_analysis.R" script.

```R
#!/drInsight_analysis.R
#The script takes the result of DEGs identified by DESeq2 tool and find drugs which can revert the expression showcased by the input.
#The present script has the following input data:
#1) Path to the directory where DESeq2 output files reside saved in the dirPath object.
#The DESeq2 file has the following columns: GeneID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
#2) Specify a pattern in order to select the right files from the directory saved in the "patternFiles" object.
#3) Path to the rank matrix file obtained from the Broad Institute saved in the "Drug_rankMatrix" object.
#The present script has the following output files:
#1) A DrInsight result file table for each DESeq2 file found (format: drInsight_inputFileName_res.csv)
#The output file has the following columns: drug, pval, FDR, pval_pass, FDR_pass.

#Load Library
library(DrInsight)	#Drug Repurposing: Integration and Systematic Investigation of Genomic High Throughput Data
library(tidyverse)	#Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gprofiler2)	#Gene list functional enrichment analysis and namespace conversion

rm(list=ls(all=TRUE)) #Remove all object saved in the system

###### USER's Input Data
dirPath <- "4_DEG_Analysis/Hsapiens"   #Path to the directory where the DESeq2 output files are saved.
patternFiles <- "DESeq2" #Specify the pattern in order to select all the input data
Drug_rankMatrix <- "rankMatrix.txt" #Drug CMap downloaded from the Broad Institute

###### Upload the rankMatrix from the Broad Institute
cmap.ref.profiles = get.cmap.ref(cmap.data.path = Drug_rankMatrix, probe.to.genes = probe.to.genes, drug.info = drug.info)

###### UPLOAD INPUT FILES
list_of_files <- list()
files <- list.files(path=dirPath, pattern=patternFiles)

if (dirPath==".")
{
  for (i in files) 
  {list_of_files[[i]] <- read.csv2 (i, header=TRUE, sep = ",")}
} else {
  for (i in files) 
  {list_of_files[[i]] <- read.csv2 (paste(dirPath, "/", i, sep = ""), header=FALSE, sep = ",")}
}

##########
for (i in 1:length(list_of_files))
{
  nowInput <- list_of_files[[i]] [,c(1,3)]
  colnames(nowInput) <- c("ID", "score")
  
  #Get Table with gene symbol and score (Log2FC)
  idConvert <- gconvert(nowInput[,1], organism = "hsapiens", target = "ENSG", numeric_ns = "", mthreshold = Inf, filter_na = TRUE) %>% arrange(input)
  nowInput <- nowInput %>% filter(ID %in% idConvert$input) %>% arrange(ID)
  nowInput <- cbind(nowInput, idConvert) %>% select(name, score)
  colnames(nowInput) <- c("geneSymbol", "score")
  
  # DrInsight Drug Identification
  drug.ident.res <- drug.ident(query.data=nowInput, cmap.ref.profiles=cmap.ref.profiles, CEG.threshold=0.05, repurposing.unit="treatment", connectivity="negative")
  drug.pvals <- drug.ident.res$drug.pvals %>% filter(pval < 0.05)
  drug.pvals$pval_pass <- 1
  drug.pvals$FDR_pass <- ifelse(drug.pvals$FDR < 0.05, 1, 0)
  
  #Save result
  write.csv(drug.pvals, file = paste("5_Drug_DrInsight_Analysis/Hsapiens/","drInsight_", names(list_of_files)[i], "_res.csv", sep = ""))
}
```

<h3 id="section8.4.">8.4.	Drug Selection </h3>

1. Each drug repurposing platform outputs a different type of file. Therefore, the drug selection was specific for the platform:

	*  DGIdb -> drug in common between the three infections were selected for further analysis.
	*  L1000CDS2 -> all identified drugs were selected from each infection analysis.
	*  DrInsight -> drugs with a pvalue lower than 0.05 were selected.

2. The selected drugs were then filtered for drugs with known antibiotics and/or antiviral properties from the following websites:

	* [Drug](https://www.drugs.com/)
	* [DrugBank](https://go.drugbank.com/)
	* [DailyMed](https://dailymed.nlm.nih.gov/dailymed/)
	* [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
	* [MedChemExpress](https://www.medchemexpress.com/)
	* [MerckMillipore](https://www.merckmillipore.com/)
	* [Tocris ](https://www.tocris.com/)

3. Literature research was conducted for each drug to identify their bacterial and viral target.
4. Drug - Protein interaction analysis. The Stitch app was used in the Cytoscape software to find known drug-protein interaction from the list of DEGs and the list of anti-infective drugs. A higher score, in the "score" column, defines a well documented interaction; whereas, the drug-protein interaction is specified by the word "pc" under the interaction column.  

---
[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtTLiBweW9nZW5lcyBBUDEgSW5mZWN0aW9uXVxuICBCW1MuIHB5b2dlbmVzIE5aMTMxIEluZmVjdGlvbl1cbiAgQ1tJbmZsdWVuemEgQSB2aXJ1cyBJbmZlY3Rpb25dXG4gIERbSC4gc2FwaWVucyBEaWZmZXJlbnRpYWwgRXhwcmVzc2VkIEdlbmVzXVxuICBBIC0tPkRcbiAgQiAtLT5EXG4gIEMgLS0-RFxuICBFW0RHSWRiIC0-IGNvbW1vbiBkcnVnc11cbiAgRltMMTAwMENEUzIgLT4gYWxsIGRydWdzXVxuICBHW0RySW5zaWdodCAtPiBwdmFsdWUgPCAwLjA1XVxuICBEIC0tPiBFXG4gIEQgLS0-IEZcbiAgRCAtLT4gR1xuICBIW0FudGktSW5mZWN0aXZlIGRydWdzIHNlbGVjdGlvbl1cbiAgRSAtLT4gSFxuICBGIC0tPiBIXG4gIEcgLS0-IEhcbiAgSVtCYWN0ZXJpYWwgYW5kIFZpcmFsIFRhcmdldCBpZGVudGlmaWNhdGlvbl1cbiAgSltTdGl0Y2ggLT4gRHJ1Zy1Qcm90ZWluIGludGVyYWN0aW9uXVxuICBIIC0tPiBJXG4gIEggLS0-IEpcblxuXG5cbiAgICAgICAgICAgICIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0IiwidGhlbWVWYXJpYWJsZXMiOnsiYmFja2dyb3VuZCI6IndoaXRlIiwicHJpbWFyeUNvbG9yIjoiI0VDRUNGRiIsInNlY29uZGFyeUNvbG9yIjoiI2ZmZmZkZSIsInRlcnRpYXJ5Q29sb3IiOiJoc2woODAsIDEwMCUsIDk2LjI3NDUwOTgwMzklKSIsInByaW1hcnlCb3JkZXJDb2xvciI6ImhzbCgyNDAsIDYwJSwgODYuMjc0NTA5ODAzOSUpIiwic2Vjb25kYXJ5Qm9yZGVyQ29sb3IiOiJoc2woNjAsIDYwJSwgODMuNTI5NDExNzY0NyUpIiwidGVydGlhcnlCb3JkZXJDb2xvciI6ImhzbCg4MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJwcmltYXJ5VGV4dENvbG9yIjoiIzEzMTMwMCIsInNlY29uZGFyeVRleHRDb2xvciI6IiMwMDAwMjEiLCJ0ZXJ0aWFyeVRleHRDb2xvciI6InJnYig5LjUwMDAwMDAwMDEsIDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxKSIsImxpbmVDb2xvciI6IiMzMzMzMzMiLCJ0ZXh0Q29sb3IiOiIjMzMzIiwibWFpbkJrZyI6IiNFQ0VDRkYiLCJzZWNvbmRCa2ciOiIjZmZmZmRlIiwiYm9yZGVyMSI6IiM5MzcwREIiLCJib3JkZXIyIjoiI2FhYWEzMyIsImFycm93aGVhZENvbG9yIjoiIzMzMzMzMyIsImZvbnRGYW1pbHkiOiJcInRyZWJ1Y2hldCBtc1wiLCB2ZXJkYW5hLCBhcmlhbCIsImZvbnRTaXplIjoiMTZweCIsImxhYmVsQmFja2dyb3VuZCI6IiNlOGU4ZTgiLCJub2RlQmtnIjoiI0VDRUNGRiIsIm5vZGVCb3JkZXIiOiIjOTM3MERCIiwiY2x1c3RlckJrZyI6IiNmZmZmZGUiLCJjbHVzdGVyQm9yZGVyIjoiI2FhYWEzMyIsImRlZmF1bHRMaW5rQ29sb3IiOiIjMzMzMzMzIiwidGl0bGVDb2xvciI6IiMzMzMiLCJlZGdlTGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsImFjdG9yQm9yZGVyIjoiaHNsKDI1OS42MjYxNjgyMjQzLCA1OS43NzY1MzYzMTI4JSwgODcuOTAxOTYwNzg0MyUpIiwiYWN0b3JCa2ciOiIjRUNFQ0ZGIiwiYWN0b3JUZXh0Q29sb3IiOiJibGFjayIsImFjdG9yTGluZUNvbG9yIjoiZ3JleSIsInNpZ25hbENvbG9yIjoiIzMzMyIsInNpZ25hbFRleHRDb2xvciI6IiMzMzMiLCJsYWJlbEJveEJrZ0NvbG9yIjoiI0VDRUNGRiIsImxhYmVsQm94Qm9yZGVyQ29sb3IiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJsYWJlbFRleHRDb2xvciI6ImJsYWNrIiwibG9vcFRleHRDb2xvciI6ImJsYWNrIiwibm90ZUJvcmRlckNvbG9yIjoiI2FhYWEzMyIsIm5vdGVCa2dDb2xvciI6IiNmZmY1YWQiLCJub3RlVGV4dENvbG9yIjoiYmxhY2siLCJhY3RpdmF0aW9uQm9yZGVyQ29sb3IiOiIjNjY2IiwiYWN0aXZhdGlvbkJrZ0NvbG9yIjoiI2Y0ZjRmNCIsInNlcXVlbmNlTnVtYmVyQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvciI6InJnYmEoMTAyLCAxMDIsIDI1NSwgMC40OSkiLCJhbHRTZWN0aW9uQmtnQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvcjIiOiIjZmZmNDAwIiwidGFza0JvcmRlckNvbG9yIjoiIzUzNGZiYyIsInRhc2tCa2dDb2xvciI6IiM4YTkwZGQiLCJ0YXNrVGV4dExpZ2h0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0RGFya0NvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dE91dHNpZGVDb2xvciI6ImJsYWNrIiwidGFza1RleHRDbGlja2FibGVDb2xvciI6IiMwMDMxNjMiLCJhY3RpdmVUYXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwiYWN0aXZlVGFza0JrZ0NvbG9yIjoiI2JmYzdmZiIsImdyaWRDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQmtnQ29sb3IiOiJsaWdodGdyZXkiLCJkb25lVGFza0JvcmRlckNvbG9yIjoiZ3JleSIsImNyaXRCb3JkZXJDb2xvciI6IiNmZjg4ODgiLCJjcml0QmtnQ29sb3IiOiJyZWQiLCJ0b2RheUxpbmVDb2xvciI6InJlZCIsImxhYmVsQ29sb3IiOiJibGFjayIsImVycm9yQmtnQ29sb3IiOiIjNTUyMjIyIiwiZXJyb3JUZXh0Q29sb3IiOiIjNTUyMjIyIiwiY2xhc3NUZXh0IjoiIzEzMTMwMCIsImZpbGxUeXBlMCI6IiNFQ0VDRkYiLCJmaWxsVHlwZTEiOiIjZmZmZmRlIiwiZmlsbFR5cGUyIjoiaHNsKDMwNCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGUzIjoiaHNsKDEyNCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIiwiZmlsbFR5cGU0IjoiaHNsKDE3NiwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU1IjoiaHNsKC00LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTYiOiJoc2woOCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU3IjoiaHNsKDE4OCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIn19LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)](https://mermaid-js.github.io/mermaid-live-editor/#/edit/eyJjb2RlIjoiZ3JhcGggVERcbiAgQVtTLiBweW9nZW5lcyBBUDEgSW5mZWN0aW9uXVxuICBCW1MuIHB5b2dlbmVzIE5aMTMxIEluZmVjdGlvbl1cbiAgQ1tJbmZsdWVuemEgQSB2aXJ1cyBJbmZlY3Rpb25dXG4gIERbSC4gc2FwaWVucyBEaWZmZXJlbnRpYWwgRXhwcmVzc2VkIEdlbmVzXVxuICBBIC0tPkRcbiAgQiAtLT5EXG4gIEMgLS0-RFxuICBFW0RHSWRiIC0-IGNvbW1vbiBkcnVnc11cbiAgRltMMTAwMENEUzIgLT4gYWxsIGRydWdzXVxuICBHW0RySW5zaWdodCAtPiBwdmFsdWUgPCAwLjA1XVxuICBEIC0tPiBFXG4gIEQgLS0-IEZcbiAgRCAtLT4gR1xuICBIW0FudGktSW5mZWN0aXZlIGRydWdzIHNlbGVjdGlvbl1cbiAgRSAtLT4gSFxuICBGIC0tPiBIXG4gIEcgLS0-IEhcbiAgSVtCYWN0ZXJpYWwgYW5kIFZpcmFsIFRhcmdldCBpZGVudGlmaWNhdGlvbl1cbiAgSltTdGl0Y2ggLT4gRHJ1Zy1Qcm90ZWluIGludGVyYWN0aW9uXVxuICBIIC0tPiBJXG4gIEggLS0-IEpcblxuXG5cbiAgICAgICAgICAgICIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0IiwidGhlbWVWYXJpYWJsZXMiOnsiYmFja2dyb3VuZCI6IndoaXRlIiwicHJpbWFyeUNvbG9yIjoiI0VDRUNGRiIsInNlY29uZGFyeUNvbG9yIjoiI2ZmZmZkZSIsInRlcnRpYXJ5Q29sb3IiOiJoc2woODAsIDEwMCUsIDk2LjI3NDUwOTgwMzklKSIsInByaW1hcnlCb3JkZXJDb2xvciI6ImhzbCgyNDAsIDYwJSwgODYuMjc0NTA5ODAzOSUpIiwic2Vjb25kYXJ5Qm9yZGVyQ29sb3IiOiJoc2woNjAsIDYwJSwgODMuNTI5NDExNzY0NyUpIiwidGVydGlhcnlCb3JkZXJDb2xvciI6ImhzbCg4MCwgNjAlLCA4Ni4yNzQ1MDk4MDM5JSkiLCJwcmltYXJ5VGV4dENvbG9yIjoiIzEzMTMwMCIsInNlY29uZGFyeVRleHRDb2xvciI6IiMwMDAwMjEiLCJ0ZXJ0aWFyeVRleHRDb2xvciI6InJnYig5LjUwMDAwMDAwMDEsIDkuNTAwMDAwMDAwMSwgOS41MDAwMDAwMDAxKSIsImxpbmVDb2xvciI6IiMzMzMzMzMiLCJ0ZXh0Q29sb3IiOiIjMzMzIiwibWFpbkJrZyI6IiNFQ0VDRkYiLCJzZWNvbmRCa2ciOiIjZmZmZmRlIiwiYm9yZGVyMSI6IiM5MzcwREIiLCJib3JkZXIyIjoiI2FhYWEzMyIsImFycm93aGVhZENvbG9yIjoiIzMzMzMzMyIsImZvbnRGYW1pbHkiOiJcInRyZWJ1Y2hldCBtc1wiLCB2ZXJkYW5hLCBhcmlhbCIsImZvbnRTaXplIjoiMTZweCIsImxhYmVsQmFja2dyb3VuZCI6IiNlOGU4ZTgiLCJub2RlQmtnIjoiI0VDRUNGRiIsIm5vZGVCb3JkZXIiOiIjOTM3MERCIiwiY2x1c3RlckJrZyI6IiNmZmZmZGUiLCJjbHVzdGVyQm9yZGVyIjoiI2FhYWEzMyIsImRlZmF1bHRMaW5rQ29sb3IiOiIjMzMzMzMzIiwidGl0bGVDb2xvciI6IiMzMzMiLCJlZGdlTGFiZWxCYWNrZ3JvdW5kIjoiI2U4ZThlOCIsImFjdG9yQm9yZGVyIjoiaHNsKDI1OS42MjYxNjgyMjQzLCA1OS43NzY1MzYzMTI4JSwgODcuOTAxOTYwNzg0MyUpIiwiYWN0b3JCa2ciOiIjRUNFQ0ZGIiwiYWN0b3JUZXh0Q29sb3IiOiJibGFjayIsImFjdG9yTGluZUNvbG9yIjoiZ3JleSIsInNpZ25hbENvbG9yIjoiIzMzMyIsInNpZ25hbFRleHRDb2xvciI6IiMzMzMiLCJsYWJlbEJveEJrZ0NvbG9yIjoiI0VDRUNGRiIsImxhYmVsQm94Qm9yZGVyQ29sb3IiOiJoc2woMjU5LjYyNjE2ODIyNDMsIDU5Ljc3NjUzNjMxMjglLCA4Ny45MDE5NjA3ODQzJSkiLCJsYWJlbFRleHRDb2xvciI6ImJsYWNrIiwibG9vcFRleHRDb2xvciI6ImJsYWNrIiwibm90ZUJvcmRlckNvbG9yIjoiI2FhYWEzMyIsIm5vdGVCa2dDb2xvciI6IiNmZmY1YWQiLCJub3RlVGV4dENvbG9yIjoiYmxhY2siLCJhY3RpdmF0aW9uQm9yZGVyQ29sb3IiOiIjNjY2IiwiYWN0aXZhdGlvbkJrZ0NvbG9yIjoiI2Y0ZjRmNCIsInNlcXVlbmNlTnVtYmVyQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvciI6InJnYmEoMTAyLCAxMDIsIDI1NSwgMC40OSkiLCJhbHRTZWN0aW9uQmtnQ29sb3IiOiJ3aGl0ZSIsInNlY3Rpb25Ca2dDb2xvcjIiOiIjZmZmNDAwIiwidGFza0JvcmRlckNvbG9yIjoiIzUzNGZiYyIsInRhc2tCa2dDb2xvciI6IiM4YTkwZGQiLCJ0YXNrVGV4dExpZ2h0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0Q29sb3IiOiJ3aGl0ZSIsInRhc2tUZXh0RGFya0NvbG9yIjoiYmxhY2siLCJ0YXNrVGV4dE91dHNpZGVDb2xvciI6ImJsYWNrIiwidGFza1RleHRDbGlja2FibGVDb2xvciI6IiMwMDMxNjMiLCJhY3RpdmVUYXNrQm9yZGVyQ29sb3IiOiIjNTM0ZmJjIiwiYWN0aXZlVGFza0JrZ0NvbG9yIjoiI2JmYzdmZiIsImdyaWRDb2xvciI6ImxpZ2h0Z3JleSIsImRvbmVUYXNrQmtnQ29sb3IiOiJsaWdodGdyZXkiLCJkb25lVGFza0JvcmRlckNvbG9yIjoiZ3JleSIsImNyaXRCb3JkZXJDb2xvciI6IiNmZjg4ODgiLCJjcml0QmtnQ29sb3IiOiJyZWQiLCJ0b2RheUxpbmVDb2xvciI6InJlZCIsImxhYmVsQ29sb3IiOiJibGFjayIsImVycm9yQmtnQ29sb3IiOiIjNTUyMjIyIiwiZXJyb3JUZXh0Q29sb3IiOiIjNTUyMjIyIiwiY2xhc3NUZXh0IjoiIzEzMTMwMCIsImZpbGxUeXBlMCI6IiNFQ0VDRkYiLCJmaWxsVHlwZTEiOiIjZmZmZmRlIiwiZmlsbFR5cGUyIjoiaHNsKDMwNCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGUzIjoiaHNsKDEyNCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIiwiZmlsbFR5cGU0IjoiaHNsKDE3NiwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU1IjoiaHNsKC00LCAxMDAlLCA5My41Mjk0MTE3NjQ3JSkiLCJmaWxsVHlwZTYiOiJoc2woOCwgMTAwJSwgOTYuMjc0NTA5ODAzOSUpIiwiZmlsbFR5cGU3IjoiaHNsKDE4OCwgMTAwJSwgOTMuNTI5NDExNzY0NyUpIn19LCJ1cGRhdGVFZGl0b3IiOmZhbHNlfQ)
---

<h2 id="section9.">9.	S. pyogenes TPM Analysis</h2>

The infected samples do not have control samples for the specific pathogen. Therefore, TPM normalisation was conducted in order to identify the highest expressed genes during each infection. Run the "count\_to\_tpm.R" script.

```R
#!/count_to_tpm.R
#The script takes a list of count tables and calculate the TPM values 
#based on the gene's length and the mean of the fragment length sequenced.
#The script has the following input:
#1) Path for the selected count table files saved in the "pathFiles" vector object.
#The count tables need to have the first column as ID followed by the samples. 
#2) Path to a text file which shows the gene's ID and its length in two columns called: ID and Length.
#3) The mean of the Fragment Length sequenced for the selected samples.
#The script outputs one file for each count table in the analysis called "star_htseq_rev_tpm.csv" 
#and it is saved in the directory the input count table is present. 
#The output has the following columns: ID, sample1, sampleN, Length, sample1_tpm, sampleN_tpm, TPM_Average.

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list = ls()) #Remove all saved objects

########### User's Input ###########

#Path to the count table files
pathFiles <- c("3_Quantification/Spyogenes_AP1/htseq_all_sample_strand_rev_count.txt", 
               "3_Quantification/Spyogenes_NZ131/htseq_all_sample_strand_rev_count.txt",
               "3_Quantification/IAV/htseq_all_sample_strand_rev_count.txt")

#The file needs to have atleast two columns called "ID" (list of gene/transcript ID) and "Length" (length of the gene/transcript).
geneLength <- "geneLength.txt"

# Mean of the Fragment Length sequenced
mean_read_seq <- 75  
################################################

###### Get the IDs' length file
ID_length <- read.table(geneLength, header = TRUE, sep = "\t")

###### TPM calculation function (https://gist.github.com/slowkow/c6ab0348747f86e2748b)
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length (NOT IN THIS CASE)
  idx <- apply(effLen, 1, function(x) min(x) > 1) 
  
  ######
  idx_true <- ID_length$ID [which(idx == "TRUE")]   #IDs where it is true
  idx_false <- ID_length$ID [which(idx == "FALSE")] #IDs where it is false
  ######
  
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  #Prepare TPM output
  colnames(tpm) <- colnames(counts) # Copy the column names from the original matrix.
  rownames(tpm) <- idx_true
  
  idx_false_df <- matrix(data = 0, nrow = length(idx_false), ncol = ncol(counts)) 
  colnames(idx_false_df) <- colnames(counts)
  rownames(idx_false_df) <- idx_false
  
  final_tpm <- as.data.frame(rbind(tpm, idx_false_df)) %>% rownames_to_column(var = "ID") %>% arrange_at(1) %>% column_to_rownames("ID")
  return(final_tpm)
}


###### Calculate TPM values for each count table
for (iFile in pathFiles) 
{
  iFile_name <- str_split(iFile, pattern = "/")[[1]]
  iCT <- read.table (iFile, header=TRUE, sep = "") %>% arrange(ID) %>% column_to_rownames(var = "ID")
  filterInfo <- ID_length %>% filter(ID %in% rownames(iCT)) %>% arrange(ID)
  
  #Calculate TPM from the count table
  meanFragmentLength <- rep(mean_read_seq, ncol(iCT))
  iCT_tpm <- counts_to_tpm (iCT, filterInfo$Length, meanFragmentLength)
  colnames(iCT_tpm) <- paste (colnames(iCT_tpm), "_tpm", sep = "")
  
  iCT$Length <- filterInfo$Length
  iCT_final <- cbind(iCT, iCT_tpm)
  
  tpm_col_index <- grep(pattern = "tpm", colnames(iCT_final))
  iCT_final$TPM_Average <- rowSums(iCT_final[tpm_col_index])/length(tpm_col_index)
  iCT_final <- iCT_final %>% arrange(colnames(iCT_final)[ncol(iCT_final)])
  iCT_final <- iCT_final %>% arrange(desc(TPM_Average))
  
  write.csv (iCT_final, file = paste("6_TPM_Analysis", "/", iFile_name[2], "/", "star_htseq_rev_tpm.csv", sep = ""))   #Save the TPM result
}
```


<h2 id="section10.">10.	Streptococcus pyogenes strains Homology Analysis </h2>

The homology analysis was conducted for the top 25 expressed genes, based on the average value of the TPM, between the two bacterial strains. 



* Select the top 25 expressed genes from both S. pyogenes AP1 and NZ131 strains based on the highest TPM average values.

```bash
$	head -25 6_TPM_Analysis/Spyogenes_AP1/star_htseq_rev_tpm.csv | awk -F"," '{print $1}' | tr -d '"' | sed '/^[[:space:]]*$/d' > 7_Homology_Analysis/AP1_selected.txt

$	head -25 6_TPM_Analysis/Spyogenes_NZ131/star_htseq_rev_tpm.csv | awk -F"," '{print $1}' | tr -d '"' | sed '/^[[:space:]]*$/d' > 7_Homology_Analysis/NZ131_selected.txt
```

* Get the proteine sequence of the selected genes from both strains.

```bash
#Linearise the proteome fasta sequence on one line
$	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.fa > Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.linearise.fa

$	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.fa > Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.linearise.fa

#Get the protein sequence of the selected genes
$	grep -A1 -f 7_Homology_Analysis/AP1_selected.txt Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.linearise.fa > 7_Homology_Analysis/AP1_selected.fa

$	grep -A1 -f 7_Homology_Analysis/NZ131_selected.txt Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.linearise.fa > 7_Homology_Analysis/NZ131_selected.fa
```

* BLAST Homology Analysis

```bash
#S. pyogenes AP1 Homology 
$	Tool/ncbi-blast-2.11.0+/bin/blastp -query 7_Homology_Analysis/AP1_selected.fa -db Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.fa -out 7_Homology_Analysis/Sp_selected_AP1_vs_NZ131_db_blast_res.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

#S. pyogenes NZ131 Homology 
$	Tool/ncbi-blast-2.11.0+/bin/blastp -query 7_Homology_Analysis/NZ131_selected.fa -db Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.fa -out 7_Homology_Analysis/Sp_selected_NZ131_vs_AP1_db_blast_res.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
``` 

* Find the target with the lowest E-value and save the protein sequence for both queries and targets. Run the "blast\_lowest\_evalue\_fasta.R" script.

```R
#!/blast_lowest_evalue_fasta.R
#The script takes the result of the BLAST analysis and it identifies for the target with the lowest evalue.
#Then, it creates a fasta file with the sequence of the query and the target.
#The script has the following input:
#1) Path to the BLAST results (in this case for AP1 and NZ131). The input file has the following columns but it is without #headers: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, stitle.
#2) Path to the proteome fasta file (in this case for AP1 and NZ131).
#The script output files are:
#1) A text file with the same structure as the BLAST result (for AP1 and NZ131) but it only showcases the target with the #lowest evalue.
#2) The output fasta file has the protein sequence for both query and target for both AP1 and NZ131 BLAST analysis.

#Load Library
library(tidyverse)	#Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(Biostrings)	#Efficient manipulation of biological strings

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
ap1_hom_res_path <- "7_Homology_Analysis/Sp_selected_AP1_vs_NZ131_db_blast_res.txt"
ap1_proteome_path <- "Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.fa"

nz131_hom_res_path <- "7_Homology_Analysis/Sp_selected_NZ131_vs_AP1_db_blast_res.txt"
nz131_proteome_path <- "Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.fa"
##### ####### ####### ######

######Proteome
ap1_proteome_fa <- readDNAStringSet(ap1_proteome_path)
nz131_proteome_fa <- readDNAStringSet(nz131_proteome_path)
ap1_nz131_proteome_fa <- c(ap1_proteome_fa, nz131_proteome_fa)

###### S. pyogenes AP1 Homology Analysis
ap1_hom_res <- read.table(ap1_hom_res_path, header=FALSE, fill = TRUE, sep = "\t")
colnames(ap1_hom_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle") 
ap1_hom_res_final <- ap1_hom_res %>% group_by(qseqid) %>% filter(evalue==min(evalue))
write.csv(ap1_hom_res_final, file = "7_Homology_Analysis/AP1_hom_blast_res_final.csv")

###### S. pyogenes NZ131 Homology Analysis
nz131_hom_res <- read.table(nz131_hom_res_path, header=FALSE, fill = TRUE, sep = "\t")
colnames(nz131_hom_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle") 
nz131_hom_res_final <- nz131_hom_res %>% group_by(qseqid) %>% filter(evalue==min(evalue))
write.csv(nz131_hom_res_final, file = "7_Homology_Analysis/NZ131_hom_blast_res_final.csv")

###### Get the sequence for both query and target from both homology analyses
ap1_queryTarget <- c(ap1_hom_res_final$qseqid, ap1_hom_res_final$stitle)
nz131_queryTarget <- c(nz131_hom_res_final$qseqid, nz131_hom_res_final$stitle)
ap1_nz131_queryTarget <- unique (c(ap1_queryTarget, nz131_queryTarget))

selIndex <- c()
for (qT in ap1_nz131_queryTarget) {selIndex <- append(selIndex, grep(qT, names(ap1_nz131_proteome_fa)))}
ap1_nz131_queryTarget_fa <- ap1_nz131_proteome_fa[selIndex]
writeXStringSet(ap1_nz131_queryTarget_fa, "7_Homology_Analysis/ap1_nz131_queryTarget.fa", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
```

* Protein domain identification

The protein domains were identified for the protein sequences identified during the homology analysis between the two S. pyogenes strains, saved in the file 7\_Homology\_Analysis/ap1\_nz131\_queryTarget.fa. The analysis was conducted in the [HMMER website](https://www.ebi.ac.uk/Tools/hmmer/).

1. Click on the "Search" tab.
2. Click on the "Upload a file" option.
3. Click on the "Choose file" and select the "7\_Homology\_Analysis/ap1\_nz131\_queryTarget.fa" fasta file.
4. Insert your email to obtain the results.
5. Click the "Submit" button.
6. Download the results from the email received upon completion of the analysis.



