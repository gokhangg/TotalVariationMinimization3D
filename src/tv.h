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
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>


template<bool IsIsotropic = true>
class TVimage
{
public:
	using ThisType = typename TVimage<IsIsotropic>;
    explicit TVimage()
    : m_stride{std::vector<size_t>(1,1)}
	, m_dim{1}
    {
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    explicit TVimage(T& size__)
    {
        setSize(size__);
        allocateMem();
    }

    template<typename... Args>
    explicit TVimage(const Args... T)
    : TVimage()
    {
        setSize(T...);
        allocateMem();
    }

    void findStride()
    {
        size_t st = 1;
        for (auto item: m_size)
        {
            st *= item;
            m_stride.emplace_back(st);
        }
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    void setSize(const T& size__)
    {
		m_dim = std::size(size__);
        m_size = size__;
		m_stride.resize(1);
		m_scale.resize(m_dim, 1.f);
		findStride();
    }
    
    template<typename... Args>
    void setSize(Args... args)
    {
		constexpr unsigned int dim = sizeof...(args);
		std::vector<size_t> size(dim);
		fillIndex(std::data(size), args...);
		setSize(size);
    }
    
    void allocateMem()
    {
        const unsigned int strd = getStride()[m_dim];
		m_cont = std::vector<float>(strd);
    }

	template<bool IsIso>
	void getDiff(const float* inP, const int size, const int m_stride, float* outIt, const float) const;

	template<>
    void getDiff<true>(const float* inP, const int size, const int m_stride, float* outIt, const float) const
    {
		auto p2 = inP + m_stride;
		auto pEnd = inP + size - m_stride;
        for (; inP < pEnd;)
            *outIt++ = *p2++ - *inP++;
    }

    template<>
    void getDiff<false>(const float* inP, const int size, const int m_stride, float* outIt, const float sc) const
    {
        auto p2 = inP + m_stride;
		auto pEnd = inP + size - m_stride;
        for (; inP < pEnd;)
            *outIt++ = sc * (*p2++ - *inP++);
    }
	
	auto begin()
	{
		return std::begin(m_cont);
	}

	auto end()
	{
		return std::end(m_cont);
	}

	auto begin() const
	{
		return std::begin(m_cont);
	}

	auto end() const
	{
		return std::end(m_cont);
	}

	auto data()
	{
		return std::data(m_cont);
	}

	auto data() const
	{
		return std::data(m_cont);
	}

	auto size() const
	{
		return std::size(m_cont);
	}

	auto getDerivative(const unsigned int axis, const bool forward) const -> ThisType
	{
		ThisType out;
        if (axis >= m_dim)
            return out;
		out = *this;
		const unsigned int strd = m_stride[axis];
		auto pBegin = std::data(m_cont);
        const auto inItEnd = std::prev(std::end(m_cont), strd);
        auto pOut = std::data(out) + (forward?0:strd);
		const auto fullSize = m_stride[std::size(m_stride) - 1];
		auto scl = m_scale[axis];
        getDiff<IsIsotropic>(pBegin, fullSize, strd, pOut, scl);

        const unsigned int upStrd = m_stride[axis + 1];
		auto ploc = std::data(out) + (forward?(upStrd - strd):0);
		const auto end_ = std::data(out) + fullSize;
        for (; ploc < end_;)
        {
            const auto endP = ploc + strd;
            for (auto ptr = ploc; ptr < endP;)
                *ptr++ = 0;
			ploc += upStrd;
        }
		return std::move(out);
	}

	auto getGradient() const -> std::vector<ThisType>
	{
		std::vector<ThisType> out(m_dim);
		for (auto ind = 0u; ind < m_dim; ++ind)
			out[ind] = getDerivative(ind, true);	
		return out;
	}

	static auto getDivergence(const std::vector<ThisType>& in) -> ThisType
	{
		ThisType out = in[0].getDerivative(0, false);
		const auto dim = out.getDim();
		for (auto ind = 1u; ind < dim; ++ind)
			out += in[ind].getDerivative(ind, false);	
		return out;
	}

	template<typename OperationT>
	ThisType  subOperatorWithSecondImage(const ThisType& inIm, OperationT operation)
	{
		ThisType out;
		if (getSize() != inIm.getSize())
			return out;
		out = *this;
		const unsigned int strd = std::size(out);
		auto ptrIn1 = std::data(*this);
		auto ptrIn2 = std::data(inIm);
		auto ptrOut = std::data(out);
		const auto endPtr = ptrIn1 + strd;
		operation(ptrIn1, ptrIn2, ptrOut, endPtr);
		return out;
	}

	template<typename T, typename OperationT>
	ThisType  subOperatorWithBasicVars(const T in, OperationT operation)
	{
		ThisType out = *this;
		const unsigned int strd = std::size(out);
		const float var = static_cast<float>(in);
		auto ptrIn = std::data(out);
		const auto endPtr = ptrIn + strd;
		operation(ptrIn, var, endPtr);
		return out;
	}

	ThisType operator+(const ThisType& inIm)
	{
		auto oper = [](float* ptrIn1, const float* ptrIn2, float* ptrOut, const float* endPtr)
		{
			for (; ptrIn1 < endPtr;)
				*ptrOut++ = *ptrIn1++ + *ptrIn2++;
		};
		return subOperatorWithSecondImage(inIm, oper);
	}

	ThisType operator-(const ThisType& inIm)
	{
		auto oper = [](float* ptrIn1, const float* ptrIn2, float* ptrOut, const float* endPtr)
		{
			for (; ptrIn1 < endPtr;)
				*ptrOut++ = *ptrIn1++ - *ptrIn2++;
		};
		return subOperatorWithSecondImage(inIm, oper);
	}

	ThisType operator*(const ThisType& inIm)
	{
		auto oper = [](float* ptrIn1, const float* ptrIn2, float* ptrOut, const float* endPtr)
		{
			for (; ptrIn1 < endPtr;)
				*ptrOut++ = *ptrIn1++ * *ptrIn2++;
		};
		return subOperatorWithSecondImage(inIm, oper);
	}

	ThisType operator/(const ThisType& inIm)
	{
		auto oper = [](float* ptrIn1, const float* ptrIn2, float* ptrOut, const float* endPtr)
		{
			for (; ptrIn1 < endPtr;)
				*ptrOut++ = *ptrIn1++ / *ptrIn2++;
		};
		return subOperatorWithSecondImage(inIm, oper);
	}

	template<typename T>
	ThisType operator+(const T in)
	{
		auto oper = [](float* ptrIn, const float var, const float* endPtr)
		{
			for (; ptrIn < endPtr;)
				*ptrIn++ = *ptrIn + var;
		};
		return subOperatorWithBasicVars(in, oper);
	}

	template<typename T>
	ThisType operator-(const T in)
	{
		auto oper = [](float* ptrIn, const float var, const float* endPtr)
		{
			for (; ptrIn < endPtr;)
				*ptrIn++ = *ptrIn - var;
		};
		return subOperatorWithBasicVars(in, oper);
	}

	template<typename T>
	ThisType operator*(const T in)
	{
		auto oper = [](float* ptrIn, const float var, const float* endPtr)
		{
			for (; ptrIn < endPtr;)
				*ptrIn++ = *ptrIn * var;
		};
		return subOperatorWithBasicVars(in, oper);
	}

	template<typename T>
	ThisType operator/(const T in)
	{
		auto oper = [](float* ptrIn, const float var, const float* endPtr)
		{
			for (; ptrIn < endPtr;)
				*ptrIn++ = *ptrIn / var;
		};
		return subOperatorWithBasicVars(in, oper);
	}

	template<typename T>
	ThisType operator=(const T in)
	{
		const float var = static_cast<float>(in);
		const unsigned int strd = std::size(m_cont);
		auto ptrIn = std::data(m_cont);
		const auto endPtr = ptrIn + strd;
		for (; ptrIn < endPtr;)
				*ptrIn++ = var;
	}

	const ThisType& operator+=(const ThisType& im)
	{
		if (getSize() != im.getSize())
			return *this;
		const unsigned int strd = std::size(m_cont);
		auto ptr = std::data(*this);
		auto ptr1 = std::data(im);
		const auto endPtr = ptr + strd;
		for (; ptr < endPtr;)
			*ptr++ = *ptr + *ptr1++;
		return *this;
	}

	template<typename... Args>
	float getVal(const Args... args)
	{
		constexpr unsigned int size = sizeof...(args);
		std::vector<size_t> res(size);
		fillIndex(std::data(res), args...);
		unsigned int index{0};
		for (auto ind = 0u; ind < std::size(m_size); ++ind)
			index += m_stride[i] * res[ind];
		return m_cont[index];
	}

    auto getSize() const -> std::vector<size_t>
    {
        return m_size;
    }

	auto getDim() const
	{
		return m_dim;
	}

	auto getStride() const -> std::vector<size_t>
    {
        return m_stride;
    }

	template<typename OperationT>
	void transform(const OperationT operation)
	{
		auto ptr = std::data(m_cont);
		const auto endPtr = ptr + std::size(m_cont);
		for (;ptr < endPtr; )
		*ptr++ = operation(*ptr);
	}

    private:
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value>::type>
	void fillIndex(size_t*& ptr, const T dim)
	{
	   *ptr++ = static_cast<size_t>(dim);
	}

	template<typename T, typename... ArgT>
	void fillIndex(size_t* ptr, const T dim, const ArgT... args)
	{
		fillIndex(ptr, dim);
		fillIndex(ptr, args...);
	}
    
    std::vector<float> m_cont;
    std::vector<size_t> m_size;
    std::vector<size_t> m_stride;
    std::vector<float> m_scale;
	unsigned int m_dim;
};

template<typename T, bool Iso>
static TVimage<Iso> operator + (T l, const TVimage<Iso>& r)
{
  return r + l;
}

template<typename T, bool Iso>
static TVimage<Iso> operator-(T l, TVimage<Iso>& r)
{
  return (r - l) * -1.f;
}

template<typename T, bool Iso>
static TVimage<Iso> operator*(T l, const TVimage<Iso>& r)
{
  return r * l;
}

template<typename T, bool Iso>
static TVimage<Iso> operator/(T l, const TVimage<Iso>& r)
{
  TVimage<Iso> out = r;
  out = 1;
  return out / r;
}




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

	template<typename T>
	vector<float> computeScaling(const T in);

	template<bool IsIso>
	TVimage<IsIso> run();
	
	template<typename T>
	T engine(T& in);

    void SetTo(float T)
    {
        To = T;
    }
    void SetIt(uint32 I)
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
