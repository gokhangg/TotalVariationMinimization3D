/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#ifndef __TV_FILTER_H__
#define __TV_FILTER_H__

#include "tv_image.h"

#include "itkImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"
#include "itkSize.h"
#include "itkMath.h"
#include "itkCastImageFilter.h"


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
    
    itkNewMacro (Self);
	itkTypeMacro(DiscreteGaussianImageFilter, ImageToImageFilter);

	/*
	@brief: Sets "to" in the algorithm.
	@param: to To value.
	@return:
	*/
    void SetTo(const float to) noexcept
    {
        m_to = to;
    }

	/*
	@brief: Sets iteration number.
	@param: it Iteration number in the optimization.
	@return:
	*/
    void SetIt(const unsigned int i) noexcept
    {
        m_it = i;
    }

	/*
	@brief: Sets "lambda" weight in the algorithm.
	@param: lam Lambda value.
	@return:
	*/
    void SetLambda(const float lam) noexcept
    {
        m_lm = lam;
    }

	/*
	@brief: Sets if the algorithm will be applied in the isotropic way or not.
	@param: iso Isotropic value.
	@return:
	*/
    void SetIsotropic(const bool iso) noexcept
    {
        m_isotropic = iso;
    }

	/*
	@brief: Sets if the algorithm will be applied to the image slice by slice.
	@param: slc Slice by slice value.
	@return:
	*/
	void SetSliceBySlice(const bool slc) noexcept
	{
		m_sliceBySlice = slc;
	}

    void PrintSelf(std::ostream & os, Indent indent) const override;

protected:
    TotalVariationMinimization()
    {
    }

    ~TotalVariationMinimization();
    void GenerateData() override;

private:
	unsigned int m_it = 10;
	float m_to = 0.15f;
	float m_lm = 0.0f;
	bool m_isotropic = false;
	bool m_sliceBySlice = false;

	template<bool IsIso>
	void run();

	/*
	@brief: Core of the algorithm.
	@param: in Input image.
	@param: out Output image to be written.
	@return:
	*/
	template<typename T>
	void engine(T& in, T& out) const;

	/*
	@brief: Computes scaling of the image dimensions.
	@param: in Pixel size vector.
	@return: Scaling vector.
	*/
	vector<float> computeScaling(const vector<float> in) const;
};

}


#ifndef ITK_MANUAL_INSTANTIATION
#include "tv_filter.hxx"
#endif
#endif
