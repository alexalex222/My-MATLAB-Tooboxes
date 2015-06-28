Memory Bounded Variational Dirichlet Process Mixture Model (MB-VDP)
Implementation of the algorithm presented in: 
R. Gomes, M. Welling, and P. Perona. Incremental Learning of Nonparametric Bayesian Mixture Models.  CVPR 2008.

This implementation is mostly Matlab with some core routines written in C and called as MEX functions.  This code is provided for 
research purposes, and is not intended for heavy duty applications.  (It ought to be re-implented in another language for that.)  
An undocumented trick is used to pass matrices by reference to the MEX helper functions, in order to avoid the unnecessary memory 
copying that Matlab tends to do. 

To install, just run install.m in Matlab. Its purpose is to invoke the MEX compiler.
The main function is called MB_VDP.  See the m-file for information on its use, and also test_script.m to see how it is invoked.

The algorithm is based on the variational approximation introduced in Kurihara et al, NIPS 2006. 
Parts of the basic Variational updates are based on the code from Kenichi Kurihara's MATLAB implementation available
on his web page.

USAGE:  You are free to use and modify this code as you see fit.  If you use it in your research or as a basis for another
algorithm or implementation, then please reference our CVPR 2008 paper in your work.  Address any questions or problems to
Ryan Gomes and gomes@vision.caltech.edu.

Copyright 2008, Ryan Gomes and The California Institute of Technology.
DISCLAIMER: No guarantees of the correct operation of this code are expressed or implied.  It has not been tested to the 
level necessary for critical applications.  The authors assume no responsibilty for the consequences of its use.