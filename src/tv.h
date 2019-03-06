/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#ifndef tv_h
#define tv_h

#include <new>
#include <math.h>
#include <iostream>
#include "itkImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"
#include "itkSize.h"
#include "itkMath.h"
#include "itkCastImageFilter.h"

using namespace std;

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class TotalVariationMinimization: public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    using Superclass = ImageToImageFilter< TInputImage, TOutputImage >;
    using Self = TotalVariationMinimization;
    using Pointer = SmartPointer< Self >;
    using ConstPointer = SmartPointer< const Self >;
    typedef unsigned int uint32;
private:
    uint32 Dm = TInputImage::ImageDimension;
    uint32 TotalVox = 1;
    uint32 It = 5;
    uint32 Sz[TInputImage::ImageDimension];
    float *Im;
    float Sp[TInputImage::ImageDimension];
    float Sc[TInputImage::ImageDimension];
    float *P[TInputImage::ImageDimension];
    float To = 0.2f;
    float Lm = 0.0f;
    bool Allocated = false;
    bool SliceBySlice = true;
    bool IsIsotropic = false;
    bool Verbose = false;
    void TVmin2D(void);
    void Engine3D(void);
    void EngineSliceBySlice(void);
public:
    itkNewMacro (Self);itkTypeMacro(DiscreteGaussianImageFilter, ImageToImageFilter);
    void Load3D(void);
    void SetTo(float T)
    {
        To = T;
    }
    void SetIt(unsigned int I)
    {
        It = I;
    }
    void SetLambda(float Lam)
    {
        Lm = Lam;
    }
    void SlcBySlc(bool Slc)
    {
        SliceBySlice = Slc;
    }
    void Isotropic(bool iso)
    {
        IsIsotropic = iso;
    }
    void SetVerbose(bool Verb)
    {
        this->Verbose = Verb;
    }
    void
    PrintSelf(std::ostream & os, Indent indent) const override;
protected:
    TotalVariationMinimization()
    {
    }
    ~TotalVariationMinimization();
    void
    GenerateData() override;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "tv.hxx"
#endif
#endif
