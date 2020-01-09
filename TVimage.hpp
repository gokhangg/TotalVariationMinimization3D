#include <vector>
#include <iostream>
#include <numeric>
#include <iterator>

class TVimage
{
public:
    
    explicit TVimage()
    : stride.emplace_back(1)
    {
    
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    explicit TVimage(T size__)
    : TVimage()
    , sz{size__}
    {
        findStride();
        allocateMem();
    }

    template<typename... Args>
    explicit TVimage(const Args... T)
    : TVimage()
    {
        setSize(T...);
        findStride();
        allocateMem();
    }

    void findStride()
    {
        size_t st = 1;
        for (auto item: sz)
        {
            st *= item;
            stride.emplace_back(st);
        }
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    void setSize(T size__)
    {
        sz = size__;
    }
    
    template<typename... Args>
    void setSize(Args... args)
    {
        sz.resize(0);
        sc.resize(0);
        stride.resize(1);
        sSize(args...);
        findStride();
    }
    
    void allocateMem()
    {
        auto multi = std::accumulate(std::begin(sz), std::end(sz), 1, std::multiplies<size_t>());
        std::cout<<multi<<std::endl;
        cont.resize(std::size(sz));
        for (auto &item: cont)
        	item.resize(multi, 0);
    }

    template<bool IsIsotropic, typename = typename std::enable_if<IsIsotropic>::type>
    forwardDerivative(const int axis) -> std::vector<float>
    {
        std::vector
    }

    template<bool IsIsotropic, typename = typename std::enable_if<!IsIsotropic>::type>
    forwardDerivative() -> std::vector<float>
    {
    }

    template<bool IsIsotropic, typename = typename std::enable_if<IsIsotropic>::type>
    backwardDerivative() -> std::vector<float>
    {
    }

    template<bool IsIsotropic, typename = typename std::enable_if<!IsIsotropic>::type>
    backwardDerivative() -> std::vector<float>
    {
    }

    
    std::vector<size_t> getSize()
    {
        return sz;
    }
    
    private:
    
    template<typename T, typename = typename std::enable_if<std::is_integral<T>::value>::type>
    void sSize(const T dim)
    {
        sz.emplace_back(static_cast<size_t>(dim));
        sc.emplace_back(1.f);
    }
    
    template<typename TT, typename... Args>
    void sSize(const TT dim, const Args... T)
    {
        sSize(T...);
        sSize(dim);
    }
    
    std::vector<std::vector<float>> cont;
    std::vector<size_t> sz;
    std::vector<size_t> stride;
    std::std::vector<float> sc;
};



    

    template<bool IsIsotropic>
    forwardDerivative(const int axis, const bool fo) -> TVimage
    {
        TVimage out;
        if (axis >= std::size(sz))
            return out;
        out.setSize(this->getSize();
        out.allocateMem();

        const unsigned int strd = stride[axis]
        const auto it1End = std::next(std::begin(cont), stride[axis+1] - strd);
        auto itOut = std::begin(out);
        getDiff<IsIsotropic>(std::begin(cont), it1End, itOut, sc[axis]);

        const unsigned int val1 = stride[axis + 1];
        const unsigned int val2 = stride[axis];
        for (auto itloc = std::next(std::begin(out), val1 - val2); itloc < std::end(out); std::advance(itLoc, val))
        {
            const auto endIt = std::next(itloc, val2);
            for (auto it = itloc; it < endIt; ++it)
                *it++ = 0;
        }
		return out;
    }

    template<bool IsIsotropic>
    backwardDerivative() -> TVimage
    {
		TVimage out;
        if (axis >= std::size(sz))
            return out;
        out.setSize(this->getSize();
        out.allocateMem();

        const unsigned int strd = stride[axis]
        const auto it1End = std::next(std::begin(cont), stride[axis+1] - strd);
        auto itOut = std::next(std::begin(out), strd);
        getDiff<IsIsotropic>(std::begin(cont), it1End, itOut, sc[axis]);

        const unsigned int val1 = stride[axis + 1];
        const unsigned int val2 = stride[axis];
        for (auto itloc = std::begin(out); itloc < std::end(out); std::advance(itLoc, val))
        {
            const auto endIt = std::next(itloc, val2);
            for (auto it = itloc; it < endIt; ++it)
                *it++ = 0;
        }
		return out;
    }


class TVimage
{
public:
    
    explicit TVimage()
    : stride.emplace_back(1)
    {
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    explicit TVimage(T size__)
    : TVimage()
    , sz{size__}
    {
        findStride();
        allocateMem();
    }

    template<typename... Args>
    explicit TVimage(const Args... T)
    : TVimage()
    {
        setSize(T...);
        findStride();
        allocateMem();
    }

    void findStride()
    {
        size_t st = 1;
        for (auto item: sz)
        {
            st *= item;
            stride.emplace_back(st);
        }
    }
    
    template<typename T, typename = typename std::enable_if<std::is_same<T, std::vector<size_t> >::value>::type>
    void setSize(T size__)
    {
        sz = size__;
		stride.resize(1);
		findStride();
    }
    
    template<typename... Args>
    void setSize(Args... args)
    {
		constexpr unsigned int num = sizeof...(args);
        sz.resize(num);
        sc.resize(num, 1.f);
        stride.resize(1);
		fillIndex(std::data(sz), args...);
    }
    
    void allocateMem()
    {
        cont.resize(stride[std::size(sz)], 0);
    }


	template<typname T, bool IsIsotropic, typename = typename std::enable_if<IsIsotropic>::type>
    inline void getDiff(T& inIt, const int stride, const T& inItEnd, T& outIt, const float)
    {
		auto it1 = inIt;
		auto it2 = std::next(inIt, stride);
        for (; it1 < inItEnd;)
            *outIt++ = *it2++ - *it1++;
    }

    template<typname T, bool IsIsotropic, typename = typename std::enable_if<!IsIsotropic>::type>
    inline void getDiff(T& inIt, const int stride, const T& inItEnd, T& outIt, const float sc)
    {
        auto it1 = inIt;
		auto it2 = std::next(inIt, stride);
        for (; it1 < inItEnd;)
            *outIt++ = sc * (*it2++ - *it1++);
    }

	auto getDerivative(const bool forward) -> TVimage
	{
		TVimage out;
        if (axis >= std::size(sz))
            return out;
        out.setSize(this->getSize();
        out.allocateMem();

		const unsigned int strd = stride[axis]
		const auto inItBegin = std::begin(cont)
        const auto inItEnd = std::prev(std::end(cont), strd);
        auto itOut = std::next(std::begin(out), forward?0:strd);
        getDiff<IsIsotropic>(inItBegin, inItEnd, strd, itOut, sc[axis]);

        const unsigned int upStrd = stride[axis + 1];
		auto itloc = std::next(std::begin(out), forward?(upStrd - strd):0);
        for (; itloc < std::end(out); std::advance(itLoc, upStrd))
        {
            const auto endIt = std::next(itloc, strd);
            for (auto it = itloc; it < endIt;)
                *it++ = 0;
        }
		return out;
	}

	template<typename... Args>
	float getVal(const Args... args)
	{
		constexpr unsigned int size = sizeof...(args);
		std::vector<size_t> res(size);
		fillIndex(std::data(res), args...)
		unsigned int index{0};
		for (auto ind = 0u; ind < std::size(sz); ++ind)
			index += stride[i] * res[ind];
		return cont[index];
	}
    
    std::vector<size_t> getSize() const -> std::vector<size_t>
    {
        return sz;
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
    
    std::vector<float> cont;
    std::vector<size_t> sz;
    std::vector<size_t> stride;
    std::std::vector<float> sc;
};