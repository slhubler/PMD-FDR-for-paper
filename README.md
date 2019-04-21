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
- Note that under the current configuration, the program will reprocess the data every time that main.r is run. Previous versions of this software saved intermediate files to speed-up later analyses but I decided to remove that functionality.
- Separate plot windows are created for each of figures 3 - 8 and supplementary figures A-D. Figures 1 and 2 were created by hand in Adobe Illustrator and PowerPoint, respectively. The data for Figure 2 and Table 4 are sent to standard output
- The figures were generated in such a way that they scale appropriately to how they will be displayed in the paper - either they will have a width of 7 inches or 3.25 inches. The windows themselves are approximately the same size but the line widths, text, etc. is set so that they all have the same size when shrunk down to the paper's size. For example, while Figures 3 and 4 have very different fonts when displayed as windows, they are both 8 point when displayed in the paper. (This is all managed through inheritance from the Plot_Image RefClass object.)

## Libraries required:
- stringr

## Code was run and tested using...
- R version 3.5.1 (2018-07-02)
- RStudio version 1.1.463

## Test code
The script found in main_test.R refers to the objects defined in test.R and to my reversion of RUnit, which is found in a different project. I left these scripts in the project in order to illustrate component use and to create a template for future test development.

## Under the hood
Here I describe the basic structure of my code.  Note that I have used a convention that makes working with this code much simpler in RStudio - I created sections using repeated # signs.  This allows me to grow and shrink sections.  For example, every class has a structure that looks like

````
####
# <class name>
####
````

Also, while R assumes that all variables are global, Ref-Classes allow some amount of encapsulation.  In order to hide "private" functions, I often define functions inside the function that uses it.  While this often appears to be no gain in terms of readability, I have found that it becomes very readable when I collapse the inside function during viewing or debugging.

Below I describe some of the primary features of this current release, from a programmer's point of view.

### Data_Processor
Data_Processor is _not_ an abstract class (although it certainly could be in a later version).  It's purpose is to allow a single entry point for the entire class dependancy hierarchy that is used for this project. The goal of this hierarchy is to perform the primary calculations from end-to-end; namely, producing an individual False Discovery Rate (FDR) for each Peptide-Spectral Match (PSM) in the input files. The relationship between all of the various objects is explicitly described and built in the `initialize()` function. Note that this is the only function within Data_Processor; all of the functionality arises from the structure of Data_Object and its implementations. 

### Data_Object
I nominally created an abstract class (which is pretty tricky in R, to be honest). The purpose of this class is to abstract dependency interactions between subclasses. The Data_Processor uses the `append_child()` directly to make a parent aware of child.  It also uses `append_parent()` (indirectly) to complete the bidirection association, through the use of functions named (for example) `set_info()` and `set_raw_data()`.

The rules for Data_Object are described in test.R under the definition of the Test_Data_Object class (I used Test-Driven Development for this piece). The key features are that you can use `load_data()` to force an object to load or `ensure()` to load data if and only if it is needed; this is a Just In Time feature. Furthermore, `load_data()` calls `ensure()` on all parents.  Thus, you need only call ensure on the object of interest, not on all of the objects on which it depends. Also, calling `ensure()` multiple times in a row only performs a single load so it doesn't hurt you. Finally, `load_data()` calls `m_load_data()`, which is implemented as an abstract method - it raises an error if the programmer has not provided a version for the subclass.

### Data_Object_Info
The top of the dependency hiearchy stored in the Data_Processor is the Data_Object_Info.  As its name implies, it is an implementation of the Data_Object.  It is mostly a data class, with a few functions thrown in for the purpose of "computing" file paths. I also consider this an abstract class, although it can be used explicitly.  The implementations in this project (Data_Object_Info_737_two_step, Data_Object_Info_737_combined, Data_Object_Pyrococcus_tr, Data_Object_Mouse_Mutations) represent the locations and desciptions of input files.

### Plot_Image
Base R has some issues when being used for plots in publication.  Rather than spend a _lot_ of time fighting with cross-platform uses of Adobe Illustrator, I decided to implement a class that ensures that an image has the correct proportions relative to itself.  This class is, once again, an abstract class that puts an image in a window of arbitrary size.  The user selects the size of the window and the subclass has been set so that it can change the size of all lines and text in proportion to the size of the window.  I have also added a standard header for my own documentation purposes but that is not used in the published images.

In terms of publication, I created the largest possible window on my machine, then copied it to the Word document.  Since I had told Plot_Image what size I wanted the final product, all of the fonts will be the correct size (at least 6 pt) and line widths will be correct (at least 1 pt, I believe). This will be true regardless of the size of the original window.  However, larger initial windows have reduced pixelation effects. The "scale" feature refers to the relationship between the initial window and the final image size; generally, larger scales (e.g. 4) better than smaller (e.g. 1 or 0.5) for reducing pixelation.

### Legend_Object
In order to standardize the image to particular font/line size, I needed an appropriate way to scale the legend (if it existed). This object manages that detail.

### Plots_for_Paper
This object creates windows for each of the plots in the published paper.  It also prints some results needed for tables and figures made by hand.

## Caveats
While we hope that others will use this work, there are a few caveats:

- This is (will be) part of a published work and should be considered copyrighted by Rhapsody Data LLC.
- The only goal of this particular project is to satisfy reproducibility requirements for the figures in the paper.
- A newer version is in the works that will will allow an implementation within Galaxy-P.

## Future Uses
All that being said, I have found Plot_Image and Legend_Object immensely useful for side-stepping a gnarly publication issue, saving me hundreds of hours in reformatting/restructuring images.  I am also happy with the current implementation of the Data_Object; it also appears to be a simple way to manage a commplex problem - organizing an ever-changing workflow.

If anyone else finds these useful, feel free to use or copy them.
