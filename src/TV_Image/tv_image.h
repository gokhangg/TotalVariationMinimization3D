/*
 * Project: 3D Total Variation minimization
 * Author: Gokhan Gunay, ghngunay@gmail.com
 * Copyright: (C) 2018 by Gokhan Gunay
 * License: GNU GPL v3 (see License.txt)
 */

#ifndef __TV_IMAGE__
#define __TV_IMAGE__

#include <vector>
#include <iterator>

template<bool IsIsotropic = true>
class TVimage
{
public:
	enum class DiffDir
	{
		FORWARD,
		BACKWARD
	};

	using ThisType = TVimage<IsIsotropic>;
	explicit TVimage()
		: m_stride{ std::vector<size_t>(1,1) }
		, m_dim{ 1 }
	{
	}

	explicit TVimage(const std::vector<size_t> size__, const float initialVal = 0.f)
		: TVimage()
	{
		setSize(size__);
		allocateMem(initialVal);
	}

	template<typename... Args>
	explicit TVimage(const Args... T)
		: TVimage()
	{
		setSize(T...);
		allocateMem();
	}

	/*
	@brief: Function to determine size of the image to be created.
	@param: size_ A vector containing size of the image.
	@return:
	*/
	void setSize(const std::vector<size_t> size__)
	{
		m_dim = std::size(size__);
		m_size = size__;
		m_stride.resize(1);
		m_scale.resize(m_dim, 1.f);
		fillStride();
	}

	/*
	@brief: Variadic template function to determine size of the image to be created.
	@param: args Arguments of the image size.
	@return:
	*/
	template<typename... Args>
	void setSize(Args... args)
	{
		m_dim = sizeof...(args);
		std::vector<size_t> size(m_dim);
                fillIndex(std::data(size), m_dim, args...);
		setSize(size);
	}

	/*
	@brief: Used to fill member stride vector.
	@param:
	@return:
	*/
	void fillStride()
	{
		size_t st = 1;
		for (auto item : m_size)
		{
			st *= item;
			m_stride.emplace_back(st);
		}
	}

	void allocateMem(const float initialVal = 0.f)
	{
		const unsigned int sz_ = getStride()[m_dim];
		m_cont = std::vector<float>(sz_, initialVal);
	}

	/*
		@brief: Used to get differentiation of data.
		@param: inP Pointer to the data container. The container must be 1d vector-like.
		@param: size Size of the data container.
		@param: stride Stride of the data to be processed. determines start and end points.
		@param: pOunt Pointer to where the result will be written.
		@param: sc Scaling factor for the differentiation.
		@return:
		*/
	void getDiff(const float* inP, const int size, const int stride, float* pOut, const float sc) const
	{
		auto p2 = inP + stride;
		auto const pEnd = inP + size - stride;
		if constexpr (IsIsotropic)
		{
			for (; inP < pEnd;)
				*pOut++ = *p2++ - *inP++;
			return;
		}
		else
		{
			for (; inP < pEnd;)
				*pOut++ = sc * (*p2++ - *inP++);
		}
	}

	/*
	@brief: Used to get derivative of data w.r.t. given axis.
	@param: axis Axis with respect of which the differentiation will be taken.
	@param: forward Determines either forward or backward differentiation.
	@return: Derived image.
	*/
	auto getDerivative(const unsigned int axis, const DiffDir diffDir) const -> ThisType
	{
		ThisType out;
		if (axis >= m_dim)
			return out;
		out = *this;
		const unsigned int strd = m_stride[axis];
		const auto pBegin = std::data(m_cont);
		const auto pOut = std::data(out) + (DiffDir::FORWARD == diffDir ? 0 : strd);
		const auto fullSize = m_stride[m_dim];
		const auto scl = m_scale[axis];
		getDiff(pBegin, fullSize, strd, pOut, scl);

		const unsigned int upStrd = m_stride[axis + 1];
		auto ploc = std::data(out) + (DiffDir::FORWARD == diffDir ? (upStrd - strd) : 0);
		const auto end_ = std::data(out) + fullSize;
		for (; ploc < end_;)
		{
			const auto endP = ploc + strd;
			for (auto ptr = ploc; ptr < endP;)
				*ptr++ = 0;
			ploc += upStrd;
		}
		return out;
	}

	/*
	@brief: Gets gradient of the image.
	@return: Gradient of the image.
	*/
	auto getGradient() const -> std::vector<ThisType>
	{
		std::vector<ThisType> out(m_dim);
		const auto dim = m_dim;
		for (auto ind = 0u; ind < dim; ++ind)
		{
			out[ind] = getDerivative(ind, ThisType::DiffDir::FORWARD);
		}
		return out;
	}

	/*
	@brief: Gets divergence of the vector images.
	This method is static.
	@param: in Input vector image.
	@return: Gradient of the vector image.
	*/
	template<typename T>
	static T getDivergence(const std::vector<T>& in)
	{
		auto out = in[0].getDerivative(0, T::DiffDir::BACKWARD);
		const auto dim = std::size(in);
		for (auto ind = 1u; ind < dim; ++ind)
		{
			out += in[ind].getDerivative(ind, T::DiffDir::BACKWARD);
		}
		return out;
	}

	/*
	@brief: Applies provided operation to two images (current and an external image).
	@param: inIm Input image.
	@param: operation Operation to be applied.
	@return: A new image obtained.
	*/
	template<typename OperationT>
	ThisType  subOperatorWithSecondImage(const ThisType& inIm, OperationT operation)
	{
		ThisType out;
		if (getSize() != inIm.getSize())
			return out;
		out = ThisType(getSize());
		const unsigned int strd = std::size(out);
		auto ptrIn1 = std::data(*this);
		auto ptrIn2 = std::data(inIm);
		auto ptrOut = std::data(out);
		const auto endPtr = ptrIn1 + strd;
		operation(ptrIn1, ptrIn2, ptrOut, endPtr);
		return out;
	}

	/*
	@brief: Applies provided operation to the image and a basic variable.
	@param: in Basic integral variable.
	@param: operation Operation to be applied.
	@return: A new image obtained.
	*/
	template<typename T, typename OperationT>
	ThisType  subOperatorWithBasicVars(const T in, OperationT operation) const
	{
		ThisType out = *this;
		const unsigned int strd = std::size(out);
		const float var = static_cast<float>(in);
		auto ptrIn = std::data(out);
		const auto endPtr = ptrIn + strd;
		operation(ptrIn, var, endPtr);
		return out;
	}

	/*
	@brief: Applies provided operation onto the image.
	@param: operation Operation to be applied.
	@return:
	*/
	template<typename OperationT>
	void transform(const OperationT operation)
	{
		auto ptr = std::data(m_cont);
		const auto endPtr = ptr + std::size(m_cont);
		for (; ptr < endPtr; )
			*ptr++ = operation(*ptr);
	}

	/*
	@brief: Sets scaling of the image dimensions.
	@param: scale Scale vector.
	@return:
	*/
	void setScaling(const std::vector<float> scale)
	{
		m_scale = scale;
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
	const ThisType& operator=(const T in)
	{
		const float var = static_cast<float>(in);
		const unsigned int strd = std::size(m_cont);
		auto ptrIn = std::data(m_cont);
		const auto endPtr = ptrIn + strd;
		for (; ptrIn < endPtr;)
			*ptrIn++ = var;
	}

	template<typename OperationT>
	void selfOperator(const ThisType& in, const OperationT oper)
	{
		const unsigned int strd = std::size(m_cont);
		auto ptrSelf = std::data(*this);
		auto ptrIn = std::data(in);
		const auto endPtr = ptrSelf + strd;
		oper(ptrSelf, ptrIn, endPtr);
	}

	const ThisType& operator+=(const ThisType& in)
	{
		auto oper = [](float* ptrSelf, const float* ptrIn, const float* const endPtr)
		{
			for (; ptrSelf < endPtr;)
				*ptrSelf++ += *ptrIn++;
		};

		if (getSize() != in.getSize())
			return *this;
		selfOperator(in, oper);
		return *this;
	}

	const ThisType& operator-=(const ThisType& in)
	{
		auto oper = [](float* ptrSelf, const float* ptrIn, const float* const endPtr)
		{
			for (; ptrSelf < endPtr;)
				*ptrSelf++ -= *ptrIn++;
		};

		if (getSize() != in.getSize())
			return *this;
		selfOperator(in, oper);
		return *this;
	}

	const ThisType& operator*=(const ThisType& in)
	{
		auto oper = [](float* ptrSelf, const float* ptrIn, const float* const endPtr)
		{
			for (; ptrSelf < endPtr;)
				*ptrSelf++ *= *ptrIn++;
		};

		if (getSize() != in.getSize())
			return *this;
		selfOperator(in, oper);
		return *this;
	}

	const ThisType& operator/=(const ThisType& in)
	{
		auto oper = [](float* ptrSelf, const float* ptrIn, const float* const endPtr)
		{
			for (; ptrSelf < endPtr;)
				*ptrSelf++ /= *ptrIn++;
		};

		if (getSize() != in.getSize())
			return *this;
		selfOperator(in, oper);
		return *this;
	}


	float& operator[](const std::size_t in)
	{
		return m_cont[in];
	}

	auto getSize() const noexcept
	{
		return m_size;
	}

	auto getDim() const noexcept
	{
		return m_dim;
	}

	auto getStride() const noexcept
	{
		return m_stride;
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

private:
	template<typename T, typename std::enable_if_t<std::is_integral<T>::value, int> = 0>
        void fillIndex(size_t*& ptr, const T dim)
	{
		*ptr++ = static_cast<size_t>(dim);
	}

	template<typename T, typename... ArgT>
        void fillIndex(size_t* ptr, const T dim, const ArgT... args)
	{
		fillIndex(ptr, dim);
                fillIndex(ptr, dim, args...);
	}

	std::vector<float> m_cont;
	std::vector<size_t> m_size;
	std::vector<size_t> m_stride;
	std::vector<float> m_scale;
	unsigned int m_dim;
};

template<bool Iso>
static TVimage<Iso> operator + (const float l, TVimage<Iso>& r)
{
	return r + l;
}

template<bool Iso>
static TVimage<Iso> operator-(const float l, TVimage<Iso>& r)
{
	TVimage<Iso> out(r.getSize, l);
	out -= r;
	return out;
}

template<bool Iso>
static TVimage<Iso> operator*(const float l, TVimage<Iso>& r)
{
	return r * l;
}

template<bool Iso>
static TVimage<Iso> operator/(const float l, TVimage<Iso>& r)
{
	TVimage<Iso> out(r.getSize, l);
	out /= r;
	return out;
}

#endif
