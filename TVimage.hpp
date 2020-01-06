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







template<class T, class TT = std::iterator<std::input_iterator_tag, typename T::value_type, int>>
class CustomIterator: public TT
{
    CustomIterator(const T it, const int strd)
    : std::iterator<std::input_iterator_tag, T::value_type, int>{it}
    , stride{strd}
    {
    }
    CustomIterator& operator++()
    {
        *this = std::next(*this, stride);
        return *this;
    }

    CustomIterator operator++(int)
    {
        CustomIterator res(*this);
        *this = std::next(*this, stride);
        return res;
    }

private:
    const int stride;
};

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
        cont.resize(stride[std::size(sz)], 0);
    }

    template<typname T, bool IsIsotropic, typename = typename std::enable_if<IsIsotropic>::type>
    inline void forwardDiff(T& inIt, const T& inEndIt, T& outIt, const float)
    {
        for (; inIt < inEndIt;)
            *outIt++ = *inIt++ - *inIt1;
    }

    template<typname T, bool IsIsotropic, typename = typename std::enable_if<!IsIsotropic>::type>
    inline void forwardDiff(T& inIt, const T& inEndIt, T& outIt, const float sc)
    {
        for (; inIt < inEndIt;)
            *outIt++ = sc * (*inIt++ - *inIt1);
    }

    template<bool IsIsotropic>
    forwardDerivative(const int axis) -> TVimage
    {
        TVimage out;
        if (axis >= std::size(sz))
            return out;
        out.setSize(this->getSize();
        out.allocateMem();

        const unsigned int strd = stride[axis]
        const auto it1End = std::next(std::begin(cont), stride[axis+1] - strd);
        auto itOut = std::begin(out);
        forwardDiff<IsIsotropic>(std::begin(cont), it1End, itOut, sc[axis]);

        const unsigned int val1 = stride[axis + 1];
        const unsigned int val2 = stride[axis];
        for (auto itloc = std::next(std::begin(out), val1 - val2); itloc < std::end(out); std::advance(itLoc, val))
        {
            const auto endIt = std::next(itloc, val2);
            for (auto it = itloc; it < endIt; ++it)
                *it++ = 0;
        }
    }

    template<bool IsIsotropic>
    backwardDerivative() -> std::vector<float>
    {
    }

    
    std::vector<size_t> getSize() const -> std::vector<size_t>
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
    
    std::vector<float> cont;
    std::vector<size_t> sz;
    std::vector<size_t> stride;
    std::std::vector<float> sc;
};