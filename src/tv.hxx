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

constexpr float EPSILON = 0.0000001f;



namespace itk
{

template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::GenerateData()
{
    uint32 dim = TInputImage::ImageDimension;

    typename TInputImage::Pointer InIm = const_cast<TInputImage *>(this->GetInput());

    /*Getting size of the image for fast processing*/
	std::vector<size_t> size(dim);
	std::vector<float> spacing(dim);
	const auto sz = InIm->GetBufferedRegion().GetSize();
	const auto sp = InIm->GetSpacing(); 
	std::cout<<sizeof(sp)/sizeof(sp[0]);
    for (int i = 0u; i < dim; i++)
    {
        size[i] = sz[i];
        spacing[i] = sp[i];
    }

    /*If the image is 3D and desired to be processed slice by slice in X-Y plane*/
    /*if (SliceBySlice)
    {
       size.resize(2);
	   spacing.resize(2);
    }*/

	this->GetOutput()->SetBufferedRegion(this->GetOutput()->GetRequestedRegion());
    this->GetOutput()->Allocate();

    /*If image slice thickness differs in each direction, get scaling weights for
     * gradient and divergence operators*/
	std::vector<float> scaling(std::size(spacing), sqrtf(1.f / 3.f));
    if (!IsIsotropic)
    {
		auto out = run<false>();
    }
	else
	{
		auto out = run<true>();
	}

    /**Update and get result*/
}

template<typename TInputImage, typename TOutputImage>
template<typename T>
vector<float> TotalVariationMinimization<TInputImage, TOutputImage>::computeScaling(const T spacing)
{
	float r = 0.f;
	const int dim_ = std::size(spacing);
	std::vector<float> scaling(dim_);
    for (int i = 0; i < dim_; i++)
		r += 1 / ((spacing[i] + EPSILON) * (spacing[i] + EPSILON));
    r = sqrt(static_cast<float>(dim_) / r);
    for (int i = 0; i < dim_; i++)
		scaling[i] = r / (spacing[i] + EPSILON);
	return scaling;
}

template<typename TInputImage, typename TOutputImage>
template<bool IsIso>
TVimage<IsIso> TotalVariationMinimization<TInputImage, TOutputImage>::run()
{
	typename TInputImage::Pointer InIm = const_cast<TInputImage *>(this->GetInput());
	uint32 dim = TInputImage::ImageDimension;

    /*Getting size of the image for fast processing*/
	std::vector<size_t> size(dim);
	std::vector<float> spacing(dim);
	const auto sz = InIm->GetBufferedRegion().GetSize();
	const auto sp = InIm->GetSpacing(); 
	std::cout<<sizeof(sp)/sizeof(sp[0]);
    for (int i = 0u; i < dim; i++)
    {
        size[i] = sz[i];
        spacing[i] = sp[i];
    }

	TVimage<IsIso> im(5,5);
	auto outIm = engine(im);
	return outIm;
}

template<typename TInputImage, typename TOutputImage>
template<typename T>
T TotalVariationMinimization<TInputImage, TOutputImage>::engine(T& in)
{
	const float lambda = 50.f;
	auto out = in;
	const float to = 0.15f;
	{
		auto vP = in.getGradient();
		const unsigned int dm = std::size(vP);
		for (int it = 0u; it < It; it++)
		{
			std::cout<<"It: "<<it<<std::endl;
			auto midP = in.getDivergence(vP) - (in / lambda);
			auto psi = midP.getGradient();
			auto r = (psi[0] * psi[0]);
			for (auto ind = 1u; ind < dm; ++ind)
				r += (psi[ind] * psi[ind]);
			r.transform(std::sqrtf);
			r = (r * to) + 1;
			for (auto ind = 0u; ind < dm; ++ind)
			{
				vP[ind] = (vP[ind] + psi[ind] * to) / r;
			}
		}
		out += (in.getDivergence(vP) * lambda);
	}
	return out;
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
}

}

#endif
