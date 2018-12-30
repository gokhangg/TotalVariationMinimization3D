/*
 * temp.cpp
 *
 *  Created on: 2018 M12 29
 *      Author: gogo
 */

/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#include "tv.h"

#ifndef tv_hxx
#define tv_hxx

#define EPSILON 0.0000001

namespace itk
{

template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::GenerateData()
{
    uint32 Dim = TInputImage::ImageDimension;
    float *pP[TInputImage::ImageDimension];
    float *pIm;
    float r = 0.0;
    typename TInputImage::Pointer InIm = const_cast<TInputImage *>(this->GetInput());

    /*Getting size of the image for fast processing*/
    for (int i = 0; i < Dim; i++)
    {
        Sz[i] = InIm->GetBufferedRegion().GetSize()[i];
        Sp[i] = InIm->GetSpacing()[i];
    }

    /*If the image is 3D and desired to be processed slice by slice in X-Y plane*/
    if (SliceBySlice)
    {
        Dm = 2;
    }
    /*If the image is 3D and desired to be processed as 3D*/
    else
    {
        Dm = Dim;
    }
    /*If image slice thickness differs in each direction, get scaling weights for
     * gradient and divergence operators*/
    if (!IsIsotropic)
    {
        for (int i = 0; i < Dm; i++)
        {
            r += 1 / ((Sp[i] + EPSILON) * (Sp[i] + EPSILON));
        }
        r = sqrt(Dm / r);
        for (int i = 0; i < Dm; i++)
        {
            Sc[i] = r / (Sp[i] + EPSILON);
        }
    }
    /*Assume the same slice thickness in each direction*/
    else
    {
        for (int i = 0; i < Dm; i++)
        {
            Sc[i] = 1;
        }
    }
    /*Getting size of the memory size to be allocated for fast processing*/
    for (int i = 0; i < Dm; i++)
    {
        TotalVox *= Sz[i] + 2;
    }
    /*mem allocation*/
    for (int i = 0; i < Dm; i++)
    {
        P[i] = new float[TotalVox];
        pP[i] = P[i];
    }
    Im = new float[TotalVox];
    Allocated = true;
    /*Initialization of the allocated memory region*/
    pIm = Im;
    for (int i = (TotalVox); i >= 0; i--)
    {
        for (int j = 0; j < Dm; j++)
            *(pP[j]++) = 0;
        *(pIm++) = 0;
    }
    this->GetOutput()->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
    this->GetOutput()->Allocate();
    /**Update and get result*/
    if (SliceBySlice || (Dm == 2))
    {
        this->EngineSliceBySlice();

    }
    else if (Dm == 3)
    {
        this->Engine3D();
    }
}
/*If image is to be processes in 3D fashion, loading the entire volume*/
template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::Load3D(void)
{
    uint32 Dim = TInputImage::ImageDimension;
    float *pIm;
    uint32 Sz2[TInputImage::ImageDimension];
    uint32 SzPl;
    ImageRegionConstIterator < TInputImage > InIt(const_cast<TInputImage *>(this->GetInput()), const_cast<TInputImage*>(this->GetInput())->GetRequestedRegion());

    /*For the image edges, zero padding is used*/
    for (int i = 0; i < Dm; i++)
    {
        Sz2[i] = Sz[i] + 2;
    }
    SzPl = Sz2[0] * Sz2[1];

    /*Loading raw image from ITK image*/
    InIt.GoToBegin();
    if (Dm == 3)
    {
        pIm = Im + SzPl + Sz2[0] + 1;
        for (int k = 0; k < Sz[2]; k++)
        {
            for (int j = 0; j < Sz[1]; j++)
            {
                for (int i = 0; i < Sz[0]; i++)
                {
                    *pIm++ = static_cast<float>(InIt.Get());
                    ++InIt;
                }
                /*padding*/
                pIm += 2;
            }
            /*padding*/
            pIm += 2 * Sz2[0];
        }
    }
}
//Reads from In image, applies 2D TV minimization and saves the result to Out back.
//In order to reduce memory use, all gradient and divergence operators applied together in the for loops.
template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::TVmin2D(void)
{
    uint32 Dim = TInputImage::ImageDimension;
    uint32 Sz2[TInputImage::ImageDimension];
    float Lam = Lm + EPSILON;
    float rat;
    float *pP[TInputImage::ImageDimension];
    float *pIm;
    float diff[TInputImage::ImageDimension];
    typename TOutputImage::Pointer OutIm = this->GetOutput();
    /*For the image boundary, padding is used*/
    Sz2[0] = Sz[0] + 2;
    Sz2[1] = Sz[1] + 2;
    /*Minimization of TV*/
    for (int it = 0; it < It; it++)
    {
        pP[0] = P[0] + 1 + Sz2[0];
        pP[1] = P[1] + 1 + Sz2[0];
        pIm = Im + Sz2[0] + 1;
        for (int j = 0; j < Sz[1]; j++)
        {
            for (int i = 0; i < Sz[0]; i++)
            {
                diff[0] = Sc[0] * ((*pP[0]) - *(pP[0] - 1)) + Sc[1] * (*(pP[1]) - *(pP[1] - Sz2[0])) - *(pIm) / Lam;
                diff[1] = diff[0];
                diff[0] = Sc[0] * (Sc[0] * (*(pP[0] + 1) - *(pP[0])) + 
                			Sc[1] * (*(pP[1] + 1) - *(pP[1] - Sz2[0] + 1)) - *(pIm + 1) / Lam - diff[0]);
                diff[1] = Sc[1] * (Sc[0] * (*(pP[0] + Sz2[0]) - *(pP[0] - 1 + Sz2[0])) + 
                			Sc[1] * (*(pP[1] + Sz2[0]) - *(pP[1])) - *(pIm + Sz2[0]) / Lam - diff[1]);
                rat = (1 + To * sqrt(diff[0] * diff[0] + diff[1] * diff[1]));
                *pP[0]++ = (*pP[0] + To * diff[0]) / rat;
                *pP[1]++ = (*pP[1] + To * diff[1]) / rat;
                pIm++;
            }
            pIm += 2;
            pP[0] += 2;
            pP[1] += 2;
        }
    }
    /*For the image boundary, padding is used*/
    pP[0] = P[0] + 1 + Sz2[0];
    pP[1] = P[1] + 1 + Sz2[0];
    pIm = Im + 1 + Sz2[0];
    /*Getting processed image*/
    for (int j = 0; j < Sz[1]; j++)
    {
        for (int i = 0; i < Sz[0]; i++)
        {
            *pIm++ = (*pIm - Lam * (Sc[0] * (*pP[0] - *(pP[0] - 1)) + Sc[1] * (*pP[1] - *(pP[1] - Sz2[0]))));
            pP[0]++;
            pP[1]++;
        }
        pIm += 2;
        pP[0] += 2;
        pP[1] += 2;
    }
}
//Reads from In image, applies 3D TV minimization and saves the result to Out back.
//In order to reduce memory use, all gradient and divergence operators applied together in the for loops.
template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::Engine3D(void)
{
    uint32 Dim = TInputImage::ImageDimension;
    uint32 Sz2[TInputImage::ImageDimension];
    uint32 SzPl;
    float rat, Lam = Lm + EPSILON;
    typename TOutputImage::Pointer OutIm = this->GetOutput();
    float *pP[TInputImage::ImageDimension];
    float *pIm;
    float diff[TInputImage::ImageDimension];
    ImageRegionIterator < TOutputImage > OutIt(OutIm, OutIm->GetRequestedRegion());
    /*load the image to the allocated mem region*/
    this->Load3D();
    /*For the image boundary, padding is used*/
    Sz2[0] = Sz[0] + 2;
    Sz2[1] = Sz[1] + 2;
    Sz2[2] = Sz[2] + 2;
    SzPl = Sz2[0] * Sz2[1];
    /*Minimization of TV*/
    for (int it = 0; it < It; it++)
    {
        pP[0] = P[0] + 1 + Sz2[0] + SzPl;
        pP[1] = P[1] + 1 + Sz2[0] + SzPl;
        pP[2] = P[2] + 1 + Sz2[0] + SzPl;
        pIm = Im + SzPl + Sz2[0] + 1;
        for (int k = 0; k < Sz[2]; k++)
        {
            for (int j = 0; j < Sz[1]; j++)
            {
                for (int i = 0; i < Sz[0]; i++)
                {
                    diff[0] = Sc[0] * ((*pP[0]) - *(pP[0] - 1)) + Sc[1] * (*(pP[1]) - *(pP[1] - Sz2[0])) + Sc[2] * (*(pP[2]) - 
                    			*(pP[2] - SzPl)) - *(pIm) / Lam;
                    diff[1] = diff[0];
                    diff[2] = diff[0];
                    diff[0] = Sc[0] * (Sc[0] * (*(pP[0] + 1) - *(pP[0])) + Sc[1] * (*(pP[1] + 1) - *(pP[1] - Sz2[0] + 1)) + 
                    			Sc[2] * (*(pP[2] + 1) - *(pP[2] - SzPl + 1)) - *(pIm + 1) / Lam - diff[0]);
                    diff[1] = Sc[1] * (Sc[0] * (*(pP[0] + Sz2[0]) - *(pP[0] - 1 + Sz2[0])) + Sc[1] * (*(pP[1] + Sz2[0]) - *(pP[1])) + 
                    			Sc[2] * (*(pP[2] + Sz2[0]) - *(pP[2] - SzPl + Sz2[0])) - *(pIm + Sz2[0]) / Lam - diff[1]);
                    diff[2] = Sc[2] * (Sc[0] * (*(pP[0] + SzPl) - *(pP[0] - 1 + SzPl)) + Sc[1] * (*(pP[1] + SzPl) - *(pP[1] + SzPl - Sz2[0])) + 
                    			Sc[2] * (*(pP[2] + SzPl) - *(pP[2])) - *(pIm + SzPl) / Lam - diff[2]);
                    rat = (1 + To * sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]));
                    *pP[0]++ = (*pP[0] + To * diff[0]) / rat;
                    *pP[1]++ = (*pP[1] + To * diff[1]) / rat;
                    *pP[2]++ = (*pP[2] + To * diff[2]) / rat;
                    pIm++;
                }
                pIm += 2;
                pP[0] += 2;
                pP[1] += 2;
                pP[2] += 2;
            }
            pIm += 2 * Sz2[0];
            pP[0] += 2 * Sz2[0];
            pP[1] += 2 * Sz2[0];
            pP[2] += 2 * Sz2[0];
        }
    }
    /*For the image boundary, padding is used*/
    pP[0] = P[0] + 1 + Sz2[0] + SzPl;
    pP[1] = P[1] + 1 + Sz2[0] + SzPl;
    pP[2] = P[2] + 1 + Sz2[0] + SzPl;
    pIm = Im + 1 + Sz2[0] + SzPl;
    /*Initialization of output image filter iterator*/
    OutIt.GoToBegin();
    /*Getting processed image*/
    for (int k = 0; k < Sz[2]; k++)
    {
        for (int j = 0; j < Sz[1]; j++)
        {
            for (int i = 0; i < Sz[0]; i++)
            {
                *pIm = (*pIm - Lam * (Sc[0] * (*pP[0] - *(pP[0] - 1)) + Sc[1] * (*pP[1] - *(pP[1] - Sz2[0])) + Sc[2] * (*pP[2] - *(pP[2] - SzPl))));
                for (int l = 0; l < Dim; l++)
                    pP[l]++;
                OutIt.Set(static_cast<typename TOutputImage::PixelType>(*pIm++));
                ++OutIt;
            }
            pIm += 2;
            pP[0] += 2;
            pP[1] += 2;
            pP[2] += 2;
        }
        pIm += 2 * Sz2[0];
        pP[0] += 2 * Sz2[0];
        pP[1] += 2 * Sz2[0];
        pP[2] += 2 * Sz2[0];
    }
}
//Reads from In image, checks if the image is 2D or 3D and applies 2D TV minimization plane by plane and saves the result to Out back.
template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::EngineSliceBySlice(void)
{
    uint32 Dim = TInputImage::ImageDimension;
    uint32 Sz2[TInputImage::ImageDimension];
    uint32 SzPl;
    float *pIm;
    typename TOutputImage::Pointer OutIm = this->GetOutput();
    ImageRegionConstIterator < TInputImage > InIt(const_cast<TInputImage*>(this->GetInput()), const_cast<TInputImage*>(this->GetInput())->GetRequestedRegion());
    ImageRegionIterator < TOutputImage > OutIt(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    /*If the image is 2D then third dim size is 1*/
    if (Dim == 2)
    {
        Sz[2] = 1;
    }
    for (int i = 0; i < Dim; i++)
    {
        Sz2[i] = Sz[i] + 2;
    }
    OutIt.GoToBegin();
    InIt.GoToBegin();

    for (int k = 0; k < Sz[2]; k++)
    {
        /*Load k-th plane to the allocated mem region*/
        pIm = Im + Sz2[0] + 1;
        for (int j = 0; j < Sz[1]; j++)
        {
            for (int i = 0; i < Sz[0]; i++)
            {
                *pIm++ = static_cast<float>(InIt.Get());
                ++InIt;
            }
            pIm += 2;
        }
        /*minimize the image in the k-th plane in the allocated mem region*/
        TVmin2D();
        /*write back the k-th plane in the allocated mem region to the out image filter mem region*/
        pIm = Im + Sz2[0] + 1;
        for (int j = 0; j < Sz[1]; j++)
        {
            for (int i = 0; i < Sz[0]; i++)
            {
                OutIt.Set(static_cast<typename TOutputImage::PixelType>(*pIm++));
                ++OutIt;
            }
            pIm += 2;
        }
    }
}
template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "Lambda: " << Lm << std::endl;
    os << indent << "Iteration Num: " << It << std::endl;
    os << indent << "To: " << It << std::endl;
    os << indent << "Iteration Num: " << It << std::endl;
}
template<typename TInputImage, typename TOutputImage>
TotalVariationMinimization<TInputImage, TOutputImage>::~TotalVariationMinimization()
{
    if (Allocated)
    {
        for (int i = 0; i < Dm; i++)
            delete[] P[i];
        delete[] Im;
        Allocated = false;
        TotalVox = 1;
    }
}

}

#endif
