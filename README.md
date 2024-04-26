# +ckc #
This should be initialized as gitmodule with name `+ckc` to use as a MATLAB package.

## Installation ##  
Navigate to the project repository root folder in a git terminal. From there, you can add this submodule using the following command:  
```(bash)
git submodule add git@github.com:Neuro-Mechatronics-Interfaces/matlab_package__ckc.git +ckc
```

## Pre-Processing Template ##  
Any `ckc.pre_process_<pipeline>.m` function is the first step in the processing pipeline for using the DEMUSE tool.  
This procedure puts the `.poly5` (or `.mat`) recording file into the correct (`.mat`) format, which can then be loaded by the experiment/application-specific `ckc.reader` function.  
Create your experiment/pipeline-specific pre-processing by copying `ckc.template__pre_process.m` into the `+ckc` folder, then renaming it as-appropriate and modifying the actual processing as appropriate based on the files that should go into the resultant `.mat` file.  

## Reader Template ##  
Any `ckc.reader_<pipeline>.m` is a CKC file-reader. It should be used with any files generated from the corresponding `ckc.pre_process_<pipeline>.m`.  
To create an experiment/pipeline-specific reader, copy `ckc.template__reader.m` to a new file in `+ckc` and name it to match your new pre-processing function.  
