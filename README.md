The source codes here are provided for 3D and 2D total variation (TV) minimizations which are widely used for image denoising and restoration. In the application, Chambolle's dual approach [1] for the minimization of the cost function was deployed.

An example application is also provided here and based on ITK. The codes are written as an ITK filter so you can also inherit them for your specific ITK applications. However, the performance is not guaranteed, since the codes have not been optimized for ITK specifications and may not comply with some ITK requirements.

There are two binaries of the example code which are compiled for Windows and Linux OSs and you can directly run them. There is also an image in the same directory on which you can see performance of the executables.

Parameters of the executables:

-in: input image

-out: output image to be written

-l: lambda (regularization weight) of the TV cost function

-it: number of iterations in the optimization

-slc: if the argument is "true" a 3D image is processed slice by slice (2D-wise), otherwise it will be processed as a whole 3D image.

-iso: if the argument is "true" a input image is processed in an isotropic fashion, otherwise slice thickness will be used to get weights of directional derivatives in the nabla operators.



[1] Chambolle, A. Journal of Mathematical Imaging and Vision (2004) 20: 89. https://doi.org/10.1023/B:JMIV.0000011325.36760.1e"# TotalVariationMinimization3D" 
