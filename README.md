# Mers Splicer: User Manual

## Input Files

The program receives two input files:
* An MGF file containing tandem mass spectrometry (MS/MS) data.
* A FASTA file containing protein sequences.

These files must be successfully uploaded before the remaining inputs on the first tab are enabled to the user.

![Input Mgf File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputmgffile.png)   ![Input Fasta File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputfasta.png)

**Figure 1: an example input MGF File (left) and Fasta File (Right).**

### Uploading a file

The files are uploaded by clicking the “Select MGF File” and “Upload FASTA File” buttons. Additional FASTA files can be uploaded by selecting the “Add Another FASTA” button. On click of these buttons, a window will open to select the relevant file. Once selected, the program will either confirm the file upload or inform you of an issues.

If you click the “Select MGF File” or “Upload FASTA File” buttons after having successfully uploaded a file, the program will behave exactly as before. If you successfully upload a second file, it will simply replace the first one. You can add as many additional FASTA files as you please, but currently there is not method of deleting them. If you upload an incorrect FASTA file via the “Add Another FASTA” button you will need to restart close and re-open the program.

The program will check that the files uploaded have the correct extensions, and will ensure that if trans output is selected the input FASTA contains at least two proteins. However, it is unlikely that other miscellaneous errors in the files will be picked up and such inconsistencies are likely to produce bugs. Figure 2 details the process for uploading a Fasta or MGF visually.

![Uploading Files](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/uploadfiles.png)

**Figure 2: infographic showing how to upload a Fasta or MGF file to the program.**

### Protein ID Syntax

The naming style for records in the FASTA file is important. The program returns the name of the protein that a peptide was found in, which it extracts from the protein ID in the input file. If you desire this data, please ensure your protein is named using the format:

				>sp|P61981|1433G_HUMAN 14-3-3 protein gamma……..

The program treats the text between the two | symbols as the name of the protein, and this is what will be included in the output files. For example, the above protein name will be simplified to P61981.

### Protein Sequences and Miscellaneous Characters

Proteins sequences which include characters that do not represent an amino-acid will not stop the program from running. Peptides are created as normal, and if they happen to contain one of these unwanted characters they are simply ignored.

### MGF Spectra Without Charge State

If some of the spectra in the input MGF file are missing the charge-state parameter, they are ignored. Spectra with charge states greater than +6 are similarly ignored.

Are there any other conditions in the input file which are ignored??

### Produce Intensity Plot

Next to the “Upload MGF File” button, there is a checkbox labelled “Produce Intensity Plot”. If this checkbox is ticked when the “Upload MGF File” button is clicked, the program will present the “Intensity Plot” to the user if an MGF file is successfully uploaded. This plot shows the distribution of maximum ion intensities of each spectra. That is, it plots ion intensity on the x-axis and the percentage of spectra in the MGF file which possess at least one ion with a greater intensity value on the y-axis. This plot can be used to choose the desired “Intensity Threshold” (discussed later in the section “First Tab User Inputs”). Figure 3 details how to create an intensity plot visually.

![Intensity Plot](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/mgfplot.png)

**Figure 3: infographic showing how to produce the intensity plot from an MGF File.**

### No MGF Comparison

Next to the “Select FASTA File” button, there is a checkbox labelled “No MGF Comparison”. Checking this box ensures that a peptide database will be created from the input FASTA file/s without filtering against an MGF File. 

## First Tab User Inputs

After a FASTA File and MGF File have been successfully uploaded (as opposed to uploading an MGF File, the user could also check the “No MGF Comparison” checkbox) the user inputs in the first tab are enabled. 

This section details what these inputs are and the impact they will have on the output peptide database. Note that valid entries must be given to these inputs before the next tab becomes available to the user.


![First tab](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/firsttab.png)

**Figure 4: image highlighting the location of the first tab user inputs being discussed (left). The three text boxes must receive valid inputs before the ‘Next Tab’ button is enabled (right).**

### PPM

The first text box on the first tab is labelled “PPM (0.1-1000)”. This text box takes a number between 0.1 and 1000 inclusive and assigns it as the PPM constant. This value is used when potential peptides are being filtered via comparison to the spectra in the MGF file, and has units of parts per million. It details how close a the monoisotopic mass of a potential peptide must be to the precursor ion mass of an MGF spectra for the two to be considered a match. 

### Intensity Threshold

The second text box on the first tab is labelled “Intensity Threshold”. This text box takes any non-negative number (including zero) and assigns it as the Intensity Threshold constant.Any spectra in the MGF file with a max ion intensity below this value are ignored in the filtering process. This is useful if the user wishes to only use high quality spectra to filter their database. The Intensity Plot detailed in the “Input Files” section is useful for selecting a desirable Intensity Threshold; the plot clearly details what proportion of the spectra will be ignored at any given Intensity Threshold. 

### Apply b/y Ion Comparison

The checkbox labelled “Apply b/y Ion Comparison” is located at the bottom of the first tab. If not checked, the filtering process only checks that the monoisotopic mass of a peptide is within the given PPM of the precursor mass of one of the spectra. If checked, a spectra and peptide pair that match on precursor mass undergo further comparison before the peptide can pass filtering. The additional test checks that a certain percentage of the b/y ions of the potential peptide have a matching ion in the spectra. The flowchart detailing the filtering process is shown in figure X. The b/y Ion comparison requires an additional two parameters, “b/y Ion Accuracy” and “Minimum b/y Ion Matches”, which will be detailed later in this section. Thus, checking the “Apply b/y Ion Comparison” checkbox also enables the text boxes which takes these inputs. 

### Minimum b/y Ion Matches and b/y Ion Accuracy

The third and fourth text boxes on the first tab are labelled “Minimum b/y Ion Matches(%)” and “b/y Ion Accuracy (0.01-0.2)”. The first of these text boxes takes any value between 0 and 100 inclusive and has units of percentage. The second takes a decimal between 0.01 and 0.2 and has units Dalton. After a peptide has matched with a spectra based on precursor mass, the user may selected to check that this spectra matches sufficiently with the peptide’s b/y ions. A b/y ion is considered matched if it’s mass and the mass of an ion in the spectra are within the prescribed “b/y Ion Accuracy”. The peptide ultimately passes the filtering stage if the percentage of b/y ions which are considered matched exceeds the value input to ““Minimum b/y Ion Matches”.

## Second Tab User Inputs

The second tab of the interface becomes available once a valid input to the first tab is received. The user navigates to the second tab by clicking the “Next Tab” button or by selecting the “Input Parameters” label in the top right corner of the window. This tab largely receives inputs which alter the nature of the peptides created and includes factors such as peptide length, modifications and splicing type. It also receives information regarding the which output files the uses wishes to produce, and which precursor charge states in the MGF file are to be considered. The process for inputting the information on the second tab is summarised in Figure 6.

![Second tab](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/secondtab.png)

**Figure 6: an infographic detailing how to input the information required on the second tab.**

### Minimum and Maximum Length

The first two dropdown boxes on the second tab are titled “Minimum Length” and “Maximum Length”. These boxes contain values from 2 to 25 and set the minimum and maximum length a peptide can be to be included in the final database. These values are inclusive; a minimum length of 2 will allow peptides of length 2 into the output. Thus if minimum and maximum are set to be equal, all output peptides will have the same length.

### Modifications

The default layout of the second tab includes three modification dropdown boxes. This can be increased to six by clicking the “Add Additional Mods” button. The program offers 32 default modifications which are detailed in the following format:
* Modification Name (ATS).

The characters in brackets represent which amino acids are modified. Some modifications which only apply to one amino acid do not have the brackets included. Table 1 summarises each of the 32 default modifications.

**Table 1: summary of the default PTMs available for selection by the user.**

| Modification Name  | Modification Name | Mass Change |
| ------------- | ------------- | ------------- |
| 4-hydroxynonenal  | HNE  | 156.11504  |
| Acetylation  | K  | 42.010567  |
| Beta-methylthiolation  | C  | 45.98772  |
| Carbamidomethylation  | C  | 57.021465  |
| Carboxylation | E  | 43.98983  |
| Carboxymethyl  | C  | 58.005478  |
| Citrullination  | R  | 0.984016  |
| Deamidation  | NQ  | 0.984016  |
| Demethylation | KR  | 28.0313  |
| Dioxidation  | M  | 31.989828  |
| FAD  | CHY  | 783.1415  |
| Farnesylation  | C  | 204.1878  |
| Geranyl-geranyl  | C  | 272.2504  |
| Guanidination  | K  | 42.021797  |
| HexNAcylation  | N  | 203.07938  |
| Hexose  | NSY  | 162.0528  |
| Lipoyl  | K  | 188.03296  |
| Acetylation  | K  | 42.010567  |
| Methylation  | KR  | 14.01565  |
| Methylation  | TSCHDE | 14.01565  |
| Oxidation  | HW  | 15.994915  |
| Oxidation  | M  | 15.994915  |
| Palmitoylation  | CSTK | 238.22966  |
| Phosphopantetheine  | S | 340.0858  |
| Phosphorylation  | HDCR  | 79.96633  |
| Phosphorylation  | STY  | 79.96633  |
| Propionamide  | C  | 71.03712  |
| Pyridoxal phosphate  | K  | 229.014  |
| S-pyridylethylation  | C | 105.057846  |
| Sulfation  | YST  | 31.989828  |
| Sulphone | M | 156.11504  |
| Ubiquitin  | TSCK | 114.04293  |
| Ubiquitination | K | 383.2281  |
