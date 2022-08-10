# PaCent
PaCent algorithm benefits from the pairwise pattern of stable-isotope unlabelled and labelled samples in data-independent acquisition MS maps for trace proteomics studies. PaCent is developed by [Guomics Lab](http://www.guomics.com).
PaCent composes of two parts. 

## Spectra Complement
R environment (3.6.0 or higher) should be installed on your computer beforehand. 
You can run the script with RStudio, or with command line as follows:
<br>Rscript.exe C:\Users\lenovo\Desktop\spectra_complement.R <input original library> <output consensus library>
<br>Example: C:\Users\R\R-3.6.0\bin\Rscript.exe C:\Users\spectra_complement.R C:\Users\tlib.tsv C:\Users\ttlib.tsv
<br>Please note that Spectra Complement currently only supports SpectraST5.0 generated libraries. For other libraries, manual modification of the dataform should be conducted.

## Pairing
Python environment (3.0 or higher) should be installed on your computer beforehand.
You can run the script with Spyder, etc., or with command line as follows:
<br>python pairing.py -rt_shift_tolerance <value> -input <file> -output <file>
<br>Example: python pairing.py -rt_shift_tolerance 10 -input openswath.tsv -output pacent.tsv.
<br>Please note that Pairing currently only supports tab separated value files from OpenSWATH2.4. For other versions of OpenSWATH, manual modification of the file headers should be conducted. 
  
