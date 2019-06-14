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

The ‘Max Mods Per Pep’ dropdown box is included under these, allowing the user to select the maximum number of PTMs that can be applied to a single peptide.

The user can concurrently select to add their own custom modification. The bottom entry in a modification dropdown box is called ‘Custom Modification’. Selecting this entry opens a pop-up window which allows the user to enter the amino acids to be modified, the associated mass change and a name for the modification. When the “Create Modification” button is clicked, a message appears either instructing the user to check the validity of their input, or instructing them that the custom modification has been successfully created. If successful, the modification will have been added to the modification dropdown box and selected. This process is detailed visually in Figure 7.

![Custom Mod 1](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/custommod1.png)
![Custom Mod 2](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/custommod2.png)

**Figure 7: infographic explaining how to create and select a new, custom modification.**

### Splice Type

The second tab also provides options for selecting the type of splicing to be performed. Separate checkboxes for linear, trans and cis are included; for each splice type selected a unique output will be created. Note that the cis splicing output will not include any peptides created via linear splicing. That is, if a peptide can be created by both linear and cis splicing of the input proteins, it will appear only in the linear output. The trans output will similarly have both cis and linear peptides removed.

The splice type selected will also limit the size of the input files the program can manage. There is no recommended maximum input size for linear splicing. For example, the full human proteome, containing approximately 20,000 protein sequences, can take between four and ten hours to run on a 16 core computer depending on the input parameters and mgf size. An input of this size is not recommended when completing cis splicing as it requires far more computational power. It is recommended that no more than 200 protein sequences be present in the input file if the output is to be completed in under 24 hours. The computation strain generated by trans splicing is even more extreme, and a maximum of 2,000 amino acids is recommended in the input Fasta file when running trans. If trans is selected and the input file contains more than 2,000 amino-acids, a message box is displayed warning the user to reconsider their input.

![Confirmation Box](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/confirmation.png)

**Figure 8: The above message box appears on screen after the ‘Generate Output’ button is clicked if trans has been selected and more than 2,000 amino acids are present in the input Fasta file.**

### Maximum Distance

The drop-down box labelled “Maximum Distance” is only relevant to cis splicing and thus is only enabled if cis splicing is selected. It contains values from 2 to 25 and a ‘None’ option. Cis spliced peptides are created when two cleavages from the same peptide are recombined together. The “Maximum Distance” sets a limit to how far apart these two cleavages can be in the origin peptide for them to be combined and included in the cis spliced output. If at their closest point, the two cleavages are not within the “Maximum Distance” the cis peptide will not be included in the output. For example, when applying a maximum distance of 3 to the sequence ASTQRALL,

* AST and ALL can be combined to create a cis peptide, but
* AST and LL cannot as they are 4 amino acids apart.

If the ‘None’ option is selected there will be no limitation of the distance between paired cleavages in the cis splicing process.

### Cis Overlap Off

This checkbox is again relevant only to cis splicing, and is thus only enabled when the cis checkbox is selected. Overlap splicing is where two cleavages which share amino-acids from their parent proteins are combined to form a new peptide. This would require the presence of two instances of the parent protein, as each amino-acid can only be used once, This type of splicing is thus by default included in the trans output, and can additionally be included in the cis output if the ‘Cis Overlap Off’ checkbox is unchecked. 

### Charge States

At the bottom of the second tab there are five checkboxes labelled +1, +2, +3, +4, +5 which are collectively titled ‘Charge States’. Checking a charge state ensures that all spectra in the MGF of that charge state will be considered when filtering. If a charge state is left unchecked, all corresponding spectra in the MGF are ignored.


## Output File Options

There are numerous default output files created when the program is run, and additional files which can be output if desired. The checkboxes which are currently unexplained relate to these files; they are discussed below and shown in figure x.

### Filtered Peptide Database Fasta Output

When the program is run, the filtered database created will be written to a Fasta file. This occurs every time the program is run and cannot be turned on or off. If multiple splice types are selected when running the program, there will be separate Filtered Peptide Database Fasta file for each splice type. 

The linear and cis filtered peptide database Fasta have the same format. Each peptide has an ID which details the record number within the output file, and the name of each protein from the input Fasta file that the peptide was found in. Under each ID the peptide sequence is listed. This is shown in figure X.

![Cis Linear Output](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/cisoutput.png)

**Figure x: Annotated snapshot of the cis and linear filtered peptide database output Fasta File.**

The trans splicing filtered peptide database Fasta file has slightly different ID convention. Expressing which proteins a trans spliced peptide is generated from is more complicated than cis and linear, as trans splicing requires to protein sequences to occur. Trans spliced peptides therefore include in their IDs a more detailed description of their origin. All possible combinations of protein sequences which can produce a given trans spliced peptide are reported in the ID. If the cleavage used to create a trans peptide is greater than 5 amino acids length, the location of this cleavage within the input protein is reported in brackets after the protein name. Overlap peptides are detailed with “Overlap/” followed by the peptide name and the range within the protein that the overlap spliced peptide derived from. This is summarised in figure x.

![Trans Output](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/transoutput.png)

**Figure x: annotated screenshot showing how to interpret the peptide ID in the trans spliced filtered peptide database Fasta.**

### To be explained

* Info file
* All user selectable options.
* Multiple outputs when memory is too large.

### Confirming the Output

Detail file naming and post output confirm program behaviour. Also detail the behaviour of the program while completing the output.

Debugging: instructions on how to check console etc.

