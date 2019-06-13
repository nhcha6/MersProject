# Mers Splicer: User Manual

## Input Files

The program receives two input files:
* An MGF file containing tandem mass spectrometry (MS/MS) data.
* A FASTA file containing protein sequences.

These files must be successfully uploaded before the remaining inputs on the first tab are enabled to the user.

![Input Mgf File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputmgffile.png)   ![Input Fasta File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputfasta.png)

** Figure 1: an example input MGF File (left) and Fasta File (Right). **

## Uploading a file

The files are uploaded by clicking the “Select MGF File” and “Upload FASTA File” buttons. Additional FASTA files can be uploaded by selecting the “Add Another FASTA” button. On click of these buttons, a window will open to select the relevant file. Once selected, the program will either confirm the file upload or inform you of an issues.

If you click the “Select MGF File” or “Upload FASTA File” buttons after having successfully uploaded a file, the program will behave exactly as before. If you successfully upload a second file, it will simply replace the first one. You can add as many additional FASTA files as you please, but currently there is not method of deleting them. If you upload an incorrect FASTA file via the “Add Another FASTA” button you will need to restart close and re-open the program.

The program will check that the files uploaded have the correct extensions, and will ensure that if trans output is selected the input FASTA contains at least two proteins. However, it is unlikely that other miscellaneous errors in the files will be picked up and such inconsistencies are likely to produce bugs. Figure 2 details the process for uploading a Fasta or MGF visually.

![Uploading Files](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/uploadfiles.png)