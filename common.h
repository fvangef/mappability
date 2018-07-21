// reduce space consumption.
// requires genome not to have more than ~ 4 gigabases
// multi-sequence fasta file must contain less than ~64k sequences of at most ~ 4 gigabases in total

namespace seqan {

template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};

};

struct SearchParams
{
    unsigned length;
    unsigned overlap;
    unsigned threads;
    static constexpr bool outputProgress = false;
    // bool indels;
    // unsigned threshold = 0;
};

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

template <typename TChar, typename TConfig>
void resetLimits(seqan::String<TChar, TConfig> const &, auto const &, unsigned const)
{ }

template <typename TString, typename TConfig>
void resetLimits(seqan::StringSet<TString, TConfig> const & text, auto & c, unsigned const length)
{
    auto const & limits = seqan::stringSetLimits(text);
    for (unsigned i = 1; i < seqan::length(limits) - 1; ++i)
    {
        for (unsigned j = 1; j < length; ++j)
        {
            c[limits[i] - j] = 0;
        }
    }
}

template <bool outputProgress>
inline void printProgress(uint64_t &, uint64_t const, uint64_t const);

template <>
inline void printProgress<false>(uint64_t &, uint64_t const, uint64_t const)
{ }

template <>
inline void printProgress<true>(uint64_t & progress_count, uint64_t const progress_step, uint64_t const progress_max)
{
    #pragma omp atomic
    ++progress_count; // TODO: necessary?
    if (omp_get_thread_num() == 0 && (progress_count & progress_step) == 0)
    {
        float progress = static_cast<float>(progress_count)/progress_max;
        std::cout << "\rProgress: " << (truncf(progress*10000)/100) << "%" << std::flush;
    }
}

template <bool outputProgress>
inline void initProgress(uint64_t &, uint64_t &, uint64_t &, uint64_t const, uint64_t const);

template <>
inline void initProgress<false>(uint64_t & progress_count, uint64_t & progress_step, uint64_t & progress_max, uint64_t const step_size, uint64_t const max_i)
{ }

template <>
inline void initProgress<true>(uint64_t & progress_count, uint64_t & progress_step, uint64_t & progress_max, uint64_t const step_size, uint64_t const max_i)
{
    progress_count = 0;
    progress_max = (max_i + step_size - 1) / step_size; // = ceil(max_i / step_size)
    progress_step = 2;
    uint64_t progress_step_tmp = std::max(progress_max/10000, static_cast<uint64_t>(1));
    while (progress_step < progress_step_tmp)
        progress_step <<= 1;
    --progress_step;
    // std::cout << "progress_max: " << progress_max << '\n';
    // std::cout << "progress_step: " << progress_step << '\n';
    // std::cout << "progress_step_tmp: " << progress_step_tmp << '\n';
    // exit(2);
}

using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;

template <typename TText>
using TIndex = seqan::Index<TText, TIndexConfig>;

// template <typename TString, typename TPos, typename TLength>
// auto extractNeedle(TString const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     return infix(text, pos, pos + len);
// }
//
// template <typename TString, typename TOwner, typename TPos, typename TLength>
// auto extractNeedle(seqan::StringSet<TString, TOwner> const & text, TPos const pos, TLength const len)
// {
//     using namespace seqan;
//     auto const & limits = stringSetLimits(text);
//     Pair<unsigned, unsigned> res;
//     posLocalize(res, pos, limits);
//     if (res.i2 + len > length(text[res.i1]))
//     {
//         TString r(suffix(text[res.i1], res.i2));
//         r += prefix(text[res.i1 + 2], length(text[res.i1]) - res.i2);
//         return r;
//     }
//     else
//     {
//         return infix(text[res.i1], res.i2, res.i2 + len);
//     }
// }

// template <typename T>
// struct is_int_vector
// {
//     static constexpr bool value = false;
// };
//
// template <uint8_t width_t>
// struct is_int_vector<sdsl::int_vector<width_t> >
// {
//     static constexpr bool value = true;
// };
//
// template <typename T>
// constexpr bool is_int_vector<T>::value;
//
// template <uint8_t width_t>
// constexpr bool is_int_vector<sdsl::int_vector<width_t> >::value;
