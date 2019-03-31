# PMD-FDR-for-paper
This is the code and datasets used in the paper to be submitted - 
Current title of paper: Challenges in Peptide-Spectrum Matching: a Robust and Reproducible Statistical Framework for Removing Low-Accuracy, High-Scoring Hits 

## Installation and running
- Pull this project to a local directory
- UnZip each of the zip files found in the Data directory back into the Data directory (IMPORTANT)
- Load the project in RStudio
- Ensure that the stringr package is installed
- Run all of the code in main.R

## Expected effects of running code
- Software will load six files from the Data directory (2 data sets have two files and 2 data sets have 1 file)
- Each data set will be processed
- After each step of processing, a "verify" step occurs; this was a rudimentary debugging step, where information about the new component was sent to standard output and to the Plots tab. If everything went correctly then this code will run (although noisily).  If an error occurred the program will crash.  This was by design to avoid creating invalid data files. Ungraceful but effective.
- Note that under the current configuration, the program will reprocess the data every time that main.r is run.  It is possible to reconfigure the code to save intermediary files, skipping the loading/processing steps.  These have been commented out.
- Separate plot windows are created for each of figures 3 - 8 and supplementary figures A-D. Figures 1 and 2 were created by hand in Adobe Illustrator and PowerPoint, respectively. The data for Figure 2 and Table 4 are sent to standard output
- The figures were generated in such a way that they scale appropriately to how they will be displayed in the paper - either they will have a width of 7 inches or 3.25 inches. The windows themselves are approximately the same size but the line widths, text, etc. is set so that they all have the same size when shrunk down to the paper's size. For example, while Figures 3 and 4 have very different fonts when displayed as windows, they are both 8 point when displayed in the paper. (This is all managed through inheritance from the Plot_Image RefClass object.)

## Libraries required:
- stringr

## Code was run and tested using...
- R version 3.5.1 (2018-07-02)
- RStudio version 1.1.463

## Caveats
While we hope that others will use this work, there are a few caveats:

- This is (will be) part of a published work and should be considered copyrighted by Rhapsody Data LLC.
- The only goal of this project is to satisfy reproducibility requirements for the figures in the paper.
- A newer, refactored version is in the works that will allow others to more easily use and maintain code.
- An offshoot of this newer software will be a module for Galaxy-P.
- While some effort has been made to simplify the code, it was the result of scientific exploration; I have deliberately left some code that has been commented out that was part of the original exploration.  The primary purpose of this is to show how to add functionality to the existing structure.
