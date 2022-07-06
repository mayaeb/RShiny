# RShiny
Shiny apps for interactive data viz and exploration 


| File Name              | Description   |
| ---------------------- | ------------- |
| finregen_RNAseq  | Application to explore differential gene expression in bulk RNAseq data comparing uninjured caudal fin tissue to regenerating tissue. Allows user to select gene of interest (searchable from dropdown list), and returns expression information and a plot of normalized counts. Data points within count plot can be clicked on to show counts. Also provides a MA plot with adjustable p-value slider. |

Open the finregen_RNAseq app by running the following code:
```
runGitHub( "RShiny", "mayaeb", subdir = "/finregen_RNAseq")
```
