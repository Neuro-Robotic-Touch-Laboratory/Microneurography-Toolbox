[![View GetWritableFolder on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/87634-getwritablefolder)

If you need to write data to a file for later use, you could store everything in the preferences. However, doing so would require you to load your data to a variable and use getpref/setpref. Especially on Octave this can be fairly slow.

This function will return the path to a folder for which Matlab/Octave has write access. This path will be dependent on the system only, not on the version of Matlab/Octave.

Licence: CC by-nc-sa 4.0