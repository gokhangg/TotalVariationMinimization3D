The source codes are provided for applications of 3D and 2D total variation (TV) minimizations which are widely used for image denoising and restoration. In the application, Chambolle's dual approach [1] for the minimization of the cost function was deployed.

An example application is provided here and based on ITK. The codes are written as an ITK filter so you can also inherit the codes for your specific ITK applications. However, the performance of the codes are not guaranteed, since they have not been optimized for ITK applications and the codes may not comply with some ITK requirements.

There are two binaries of the example code which are compiled for Windows and Linux and you can directly run.

Parameters of the executables:
-in: input image
-out: output image to be written
-l: lambda (regularization weight) of the TV cost function
-it: number of iterations in the optimization
-slc: if the argument is "true" a 3D image is processed slice by slice (2D-wise), otherwise it will be processed as a whole 3D image.


[1] Chambolle, A. Journal of Mathematical Imaging and Vision (2004) 20: 89. https://doi.org/10.1023/B:JMIV.0000011325.36760.1e"# TotalVariationMinimization3D" 
