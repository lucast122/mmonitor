## MMonitor quick Tutorial with test data

This shall serve as a small example on how to use MMonitor with some test data. We will be doing taxonomic analysis of 6 16S rRNA samples using a small coverage, followed by an evaluation using the MMonitor dashboard including some basic statistics and metadata correlations.

1. Go to the website www.mmonitor.de and register a new user. Remember the username and password for later.
2. Clone this repository:
   `git clone https://github.com/lucast122/MMonitor`
3. There should be a filled called 'sample_list.csv' in the MMonitor/desktop/mmonitor_tutorial subfolder. Change the paths in the 'sample folder' column so they match the actual path of the files on your    computer. For example in the row describing sample1 you should enter the full path to the folder sample 1, do the same for the other samples.   
4. Create a new conda environment, install the required packages and change PYTHONPATH:  
   `conda create --name mmonitor python=3.11`  
   `conda activate mmonitor`  
   `pip install -r MMonitor/desktop/requirements.txt`  
   `conda install -c bioconda minimap2`  
   `export PYTHONPATH=$PYTHONPATH:MMonitor/desktop/`  
   `export PYTHONPATH=$PYTHONPATH:MMonitor/desktop/src/`  
   
6. Start the GUI using `python MMonitor/desktop/src/mmonitor/__main__.py` If you get an error telling you that some module was not you probably didn't export the PYTHONPATH correctly. The command uses a relative path to the folder that you downloaded, so maybe change it to the full absolute path if the relative path doesn't work.
7. After starting the app the main window should pop up. Click user authentication and change the username and password to the name and password you chose on the website. Save the config.
8. Click 'Process sequencing data' and select 'Quick taxonomy 16s nanopore' then click 'continue'.
9. The 'Sample Data Input' window should open. Now click 'Add multiple samples from CSV'.
10. In the file selection window select the sample_list.csv that you changes earlier.
11. You should receive the message 'Processing 6 samples. This may take a while.' if everything worked. If you receive an error message check the file paths in the csv again and make sure they are   
    correct.
12. Click 'submit' and check the console. Analysis should start now. If you receive an error message make sure that minimap2 is installed an on your sytem PATH. Wait for the analysis to finish, you  
    should see a popup window saying 'Analysis complete'.
13. Go back to the website, login with your credentials and click on 'Dashboard'.
14. Check out the taxonomy of the samples in the taxonomy dashboard. Try out different plot types and different taxonomic ranks by using the dropdown selections.
15. Go to diversity and check the alpha diversity of the samples.
16. Go to QC and check basic statistics of the samples.
17. Go to Correlations. Upload the 'metadata.csv' from the folder 'mmonitor_tutorial'. If the metadata parsed correctly you should get a notification telling you that 6 metadata entries were parsed.
18. Refresh the page and click 'Download Taxonomy-metadata Correlations' to download a CSV file containing the correlation between the taxonomic profiles and the metadata. Check out the resulting   
    correlations to see if any taxa correlate with the metadata.
19. If at any stage something didn't work please report the error using the 'Submit feedback' function. Please also report general feedback and other issues you had using the software.




## App Requirements

- Install `R` and packages `("jpeg", "png", "RColorBrewer", "lattice", "latticeExtra")`
- Add `Rscript` to path!
- Download the latest version of [centrifuge](https://github.com/DaehwanKimLab/centrifuge), install and make sure executable is on system path
- Either download and [build new index](https://ccb.jhu.edu/software/centrifuge/manual.shtml) or download pre-built index for centrifuge


## Develop Project

- Install Python>=3.9 (must include `tkinter`)
- `pip install -r requirements.txt`
- `cd src && python -m mmonitor`


## R Horizon Plot Troubleshooting

`latticeExtra` and its dependencies needs to be installed in the library that `Rscript` uses.
Find out what libraries `Rscript` uses with `Rscript -e ".libPaths()"` and make sure they all have the necessary packages installed with
`install.packages(c("jpeg", "png", "RColorBrewer", "lattice", "latticeExtra"), lib = "<library path>")`.
You might need to start your R IDE with admin privileges depending on the library location. 


## Build App


1. Download code on machine with desired operating system 
2. `pip install pyinstaller`[1](https://pyinstaller.readthedocs.io/en/stable/spec-files.html) [2](https://pyinstaller.readthedocs.io/en/stable/usage.html#supporting-multiple-operating-systems)

3. `python build.py`

Tested on Windows 10 and Ubuntu 20.

Executable is most likely not backwards compatible.

## How to use

1. Create or choose Project
  A project represents the database MMonitor uses. The user always has to create or choose a previously created database in order to use the tool. Creating a project         will create a called <name>.sqlite3 on the dard drive, where <name> is the name the user chooses.
2. Run analysis pipeline
  Click the button "Run analysis pipeline" to open up the pipeline selection window. Please select the analysis you want to perform. For taxonomic analysis select either     "Quick taxonomy nanopore" or "Quick taxonomy 16s nanopore" depending on your data input. Use the first if you have whole genome sequencing reads and, if you have 16s       reads choose the second. The functional analysis takes way longer than taxonomic because the reads will first need to get assembled into contigs and then mapped against    a database for binning. Especially binning may take longer depending on your system.
4. Start Monitoring
   If you choose a project that has been previously created and you already ran an analysis for it then you can also skip step 2 and start monitoring directly. However,       clicking the button "Start monitoring" with an empty project will result in an error message popping up. If you ran an analysis in step 2 and received a notification       that the analysis has finished you can then go ahead and click "Start monitoring". A browser window with your default browser will popup along with the visualization       application.
5. Analyze data in Dashboard
   After step 3 you can now look at the results of the analysis in the browser dashboard. Click on one of the links to change between the different apps. Following is a       quick summary of each apps' functionality.
Optional: Add metadata from csv
   If you have metadata you can add it as a csv file. In the analysis dashapp you can then see correlation values for your metadata and taxonomy.

Taxonomy: View the absolute and relative abundances of your samples after running a taxonomic analysis.
Horizon plot: Create a horizon plot to quickly see changes in abundance between your samples.
Correlations: If you previously added metadata you can here view the different correlation values between your metagenome abundances and the metadata.
KEGG: View the KEGG metabolic maps created by the analysis pipeline if all functional analysis steps were selected.
Gene browser: Check out your assembled metagenomes using the gene browser nNot tested yet, may not work on every system).

Be sure to press 'Quit' to terminate the app, as it also terminates the dash server.

