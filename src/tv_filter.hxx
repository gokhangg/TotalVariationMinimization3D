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

#include "tv_filter.h"

#ifndef tv_hxx
#define tv_hxx

constexpr float EPSILON = 0.0000001f;

template<typename T1, typename T2>
auto mv(T1 pIn, const size_t sz, T2* pOut)
{
	T1 const _end_ = pIn + sz;
	for (; pIn < _end_;)
	{
		*pOut++ = static_cast<T2>(*pIn++);
	}
};

namespace itk
{

template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::GenerateData()
{
    

	/*Is image to be processed isotropically?*/
	if (!m_isotropic)
	{
		run<true>();
	}
	else
	{
		run<false>();
	}
    /**Update and get result*/
}

template<typename TInputImage, typename TOutputImage>
auto TotalVariationMinimization<TInputImage, TOutputImage>::computeScaling(const std::vector<float> spacing) const -> std::vector<float>
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
void TotalVariationMinimization<TInputImage, TOutputImage>::run() 
{
	typename TInputImage::Pointer InIm = const_cast<TInputImage *>(this->GetInput());
	unsigned int dim = TInputImage::ImageDimension;

	/*Getting size of the image for fast processing*/
	std::vector<size_t> size_(dim);
	std::vector<float> spacing(dim);
	const auto sz = InIm->GetBufferedRegion().GetSize();
	const auto sp = InIm->GetSpacing();
	for (auto i = 0u; i < dim; i++)
	{
		size_[i] = sz[i];
		spacing[i] = sp[i];
	}

	unsigned int cnt = 1;
	if (m_sliceBySlice && dim > 2)
	{
		for (auto ind = 2u; ind < dim; ++ind)
			cnt *= size_[ind];
		spacing.resize(2);
		size_.resize(2);
	}
	
	auto pIn = this->GetInput()->GetBufferPointer();
	
	auto out = this->GetOutput();
	out->SetBufferedRegion(out->GetRequestedRegion());
	out->Allocate();
	auto pOut = out->GetBufferPointer();

	/*If image slice thickness differs in each direction, get scaling weights for
	* gradient and divergence operators*/
	auto scaling = computeScaling(spacing);

	TVimage<IsIso> im(size_);
	im.setScaling(scaling);
	
	for (auto ind = 0u; ind < cnt; ++ind)
	{
		const auto sz_ = std::size(im);
		mv(pIn + ind * sz_, sz_, std::data(im));
		engine(im, im);
		mv(std::data(im), sz_, pOut + ind * sz_);
	}
}

template<typename TInputImage, typename TOutputImage>
template<typename T>
void TotalVariationMinimization<TInputImage, TOutputImage>::engine(T& in, T& out) 
{
	const float lambda = EPSILON + m_lm;
	const float to = EPSILON + m_to;
	
	out = in / lambda;
	auto vP = in.getGradient();
	auto midP = TVimage<>::getDivergence(vP) - out;
	auto psi = midP.getGradient();
	auto r = (psi[0] * psi[0]);

	const auto dm = std::size(vP);
	const auto nIt = m_it - 1;
	for (auto it = 0u; it < nIt; it++)
	{
		midP = TVimage<>::getDivergence(vP) - out;
		psi = midP.getGradient();
		for (auto ind = 0u; ind < dm; ++ind)
			r += (psi[ind] * psi[ind]);
		r.transform(std::sqrtf);
		r = (r * to) + 1;
		for (auto ind = 0u; ind < dm; ++ind)
		{
			vP[ind] = (vP[ind] + psi[ind] * to);
			vP[ind] /= r;
		}
	}
	out -= out.getDivergence(vP);
	auto oper = [&](const float in)
	{
		return lambda * in;
	};
	out.transform(oper);

}

template<typename TInputImage, typename TOutputImage>
void TotalVariationMinimization<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "Lambda: " << m_lm << std::endl;
    os << indent << "Iteration Num: " << m_it << std::endl;
    os << indent << "To: " << m_it << std::endl;
    os << indent << "Iteration Num: " << m_it << std::endl;
}
template<typename TInputImage, typename TOutputImage>
TotalVariationMinimization<TInputImage, TOutputImage>::~TotalVariationMinimization()
{
}

}

#endif
