# fourgp_speclib

This python package defines classes representing arrays of spectra, and libraries to keep them in.

Spectrum libraries are a bit like having a directory full of data files on disk, each containing a spectrum. However, they also include a database which can store arbitrary metadata about each spectrum -- for example, stellar paramaters and abundances. It is possible to search a spectrum library based on metadata constraints.

Various implementations of the SpectrumLibrary class are provided, storing the metadata in different flavours of SQL database. SQLite is probably the simplest and creates portable libraries that you can transfer to a different machine with all metadata intact. MySQL is a faster database engine, and probably a better option for data which doesn't need to move around.

# Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden
