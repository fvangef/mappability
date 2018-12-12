using namespace seqan;

template <typename TIter, typename value_type, typename TVector, typename TText>
inline void extendExact(TIter it, value_type * hits, TVector & hits_occ, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();

    if (b - a + 1 == length)
    {
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        for (auto occ:getOccurrences(it))
            hits_occ.at(a-ab).push_back(occ);
        return;
    }
    //if (b + 1 <= bb)
    //{
        auto it2 = it;
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            bool success = true;
            for (uint64_t i = b + 1; i <= b_new && success; ++i)
            {
                success = goDown(it2, text[i], Rev());
            }
            if (success)
                extendExact(it2, hits, hits_occ, text, length, a, b_new, ab, bb);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, (int64_t)0);
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if(!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, hits_occ, text, length, a_new, b, ab, bb);
    }
}

// forward
template <typename TIter, typename value_type, typename TVector, typename TText>
inline void extend(TIter it, value_type * hits, TVector & hits_occ, unsigned errors_left, TText const & text,
            unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb // entire interval
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <typename TIter, typename value_type, typename TVector, typename TText>
inline void approxSearch(TIter it, value_type * hits, TVector & hits_occ, unsigned errors_left, TText const & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const &
)
{
    if (b == b_new)
    {
        extend(it, hits, hits_occ, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch(it, hits, hits_occ, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
            } while (goRight(it, Rev()));
        }
    }
    else
    {
        for (uint64_t i = b + 1; i <= b_new; ++i)
        {
            if (!goDown(it, text[i], Rev()))
                return;
        }
        extendExact(it, hits, hits_occ, text, length, a, b_new, ab, bb);
    }
}
template <typename TIter, typename value_type, typename TVector, typename TText>
inline void approxSearch(TIter it, value_type * hits, TVector & hits_occ, unsigned errors_left, TText const & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/
)
{
    if (a == a_new)
    {
        extend(it, hits, hits_occ, errors_left, text, length, a, b, ab, bb);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch(it, hits, hits_occ, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
            } while (goRight(it, Fwd()));
        }
    }
    else
    {
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if (!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact(it, hits, hits_occ, text, length, a_new, b, ab, bb);
    }
}

template <typename TIter, typename value_type, typename TVector, typename TText>
inline void extend(TIter it, value_type * hits, TVector & hits_occ, unsigned errors_left, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
)
{
    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();
    if (errors_left == 0)
    {
        extendExact(it, hits, hits_occ, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        for (auto occ:getOccurrences(it))
            hits_occ.at(a-ab).push_back(occ);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch(it, hits, hits_occ, errors_left, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev()
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, (int64_t)0);
        approxSearch(it, hits, hits_occ, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, typename TIndex, typename TText, typename TContainer>
inline void runAlgo2(TIndex & index, TText const & text, TContainer & c, SearchParams const & params)
{
    typedef typename TContainer::value_type value_type;

    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, params.overlap);

    uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string

    const uint64_t max_i = textLength - params.length + 1;
    const uint64_t step_size = params.length - params.overlap + 1;
   // #pragma omp parallel for schedule(dynamic, std::max(1ul, max_i/(step_size*params.threads*50))) num_threads(params.threads)
    for (uint64_t i = 0; i < max_i; i += step_size)
    {
        value_type hits[params.length - params.overlap + 1] = {};
        auto delegate = [&hits, i, textLength, &params, &text](auto it, auto const & /*read*/, unsigned const errors_spent) {
            uint64_t const bb = std::min(textLength - 1, i + params.length - 1 + params.length - params.overlap);
            extend(it, hits, errors - errors_spent, text, params.length,
                i + params.length - params.overlap, i + params.length - 1, // searched interval
                i, bb // entire interval
            );
        };

        auto const & needle = infix(text, i + params.length - params.overlap, i + params.length);
        Iter<TIndex, VSTree<TopDown<> > > it(index);
        _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
        uint64_t max_pos = std::min(i + params.length - params.overlap, textLength - params.length);
        for (uint64_t j = i; j <= max_pos; ++j)
            c[j] = hits[j - i];
    }

    resetLimits(indexText(index), c, params.length);
}
