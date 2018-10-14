using namespace seqan;

template <unsigned max_errors, typename TIter, typename value_type, typename TText>
inline void extendExact3(TIter it, value_type * hits, TIter * it_zero_errors, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
    , std::vector<std::vector<TIter>> & one_error_iters, int * one_error_iters_count
)
{
    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();

    if (b - a + 1 == length)
    {
        if (max_errors == 0)
            it_zero_errors[a-ab] = it; // TODO: könnte man für max_errors > 0 rausoptimieren
        else if (max_errors == 1)
        {
            one_error_iters[a-ab].push_back(it);
            one_error_iters_count[a-ab] = std::min((uint64_t) countOccurrences(it) + one_error_iters_count[a-ab], max_val);
        }
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
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
                extendExact3<max_errors>(it2, hits, it_zero_errors, text, length, a, b_new, ab, bb, one_error_iters, one_error_iters_count);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, 0l);
        for (int64_t i = a - 1; i >= a_new; --i)
        {
            if(!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a_new, b, ab, bb, one_error_iters, one_error_iters_count);
    }
}

// forward
template <unsigned max_errors, typename TIter, typename value_type, typename TText>
inline void extend3(TIter it, value_type * hits, TIter * it_zero_errors, unsigned errors_left, TText const & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb // entire interval
            , std::vector<std::vector<TIter>> & one_error_iters, int * one_error_iters_count
);

// TODO: remove text everywhere: auto & text = indexText(index(it));
template <unsigned max_errors, typename TIter, typename value_type, typename TText>
inline void approxSearch3(TIter it, value_type * hits, TIter * it_zero_errors, unsigned errors_left, TText const & text, unsigned const length,
            uint64_t a, uint64_t b, // searched interval
            uint64_t ab, uint64_t bb, // entire interval
            uint64_t b_new,
            Rev const &
            , std::vector<std::vector<TIter>> & one_error_iters, int * one_error_iters_count
)
{
    if (b == b_new)
    {
        extend3<max_errors>(it, hits, it_zero_errors, errors_left, text, length, a, b, ab, bb, one_error_iters, one_error_iters_count);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left - delta, text, length, a, b + 1, ab, bb, b_new, Rev(), one_error_iters, one_error_iters_count);
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
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a, b_new, ab, bb, one_error_iters, one_error_iters_count);
    }
}
template <unsigned max_errors, typename TIter, typename value_type, typename TText>
inline void approxSearch3(TIter it, value_type * hits, TIter * it_zero_errors, unsigned errors_left, TText const & text, unsigned const length,
                  uint64_t a, uint64_t b, // searched interval
                  uint64_t ab, uint64_t bb, // entire interval
                  int64_t a_new,
                  Fwd const & /*tag*/
                  , std::vector<std::vector<TIter>> & one_error_iters, int * one_error_iters_count
)
{
    if (a == a_new)
    {
        extend3<max_errors>(it, hits, it_zero_errors, errors_left, text, length, a, b, ab, bb, one_error_iters, one_error_iters_count);
        return;
    }
    if (errors_left > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left - delta, text, length, a - 1, b, ab, bb, a_new, Fwd(), one_error_iters, one_error_iters_count);
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
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a_new, b, ab, bb, one_error_iters, one_error_iters_count);
    }
}

template <unsigned max_errors, typename TIter, typename value_type, typename TText>
inline void extend3(TIter it, value_type * hits, TIter * it_zero_errors, unsigned errors_left, TText const & text, unsigned const length,
    uint64_t a, uint64_t b, // searched interval
    uint64_t ab, uint64_t bb // entire interval
    , std::vector<std::vector<TIter>> & one_error_iters, int * one_error_iters_count
)
{
    constexpr uint64_t max_val = std::numeric_limits<value_type>::max();

    if (errors_left == 0)
    {
        extendExact3<max_errors>(it, hits, it_zero_errors, text, length, a, b, ab, bb, one_error_iters, one_error_iters_count);
        return;
    }
    if (b - a + 1 == length)
    {
        if (max_errors == errors_left)
            it_zero_errors[a-ab] = it;
        else if (max_errors == 1)
        {
            one_error_iters[a-ab].push_back(it);
            one_error_iters_count[a-ab] = std::min((uint64_t) countOccurrences(it) + one_error_iters_count[a-ab], max_val);
        }
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev(), one_error_iters, one_error_iters_count
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        int64_t a_new = alm + std::max((int64_t) ((a - alm) - 1) >> 1, 0l);
        approxSearch3<max_errors>(it, hits, it_zero_errors, errors_left, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd(), one_error_iters, one_error_iters_count
        );
    }
}

// template <unsigned errors, typename TIndex, typename TText, typename TContainer>
// inline void runAlgo3(TIndex & index, TText const & text, TContainer & c, SearchParams const & params)
// {
//     typedef typename TContainer::value_type value_type;
//     typedef Iter<TIndex, VSTree<TopDown<> > > TIter;
//
//     auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
//     _optimalSearchSchemeComputeFixedBlocklength(scheme, params.overlap);
//
//     auto const & limits = stringSetLimits(indexText(index));
//
//     uint64_t const textLength = length(text); // lengthSum() forwards to length() for a single string
//
//     const uint64_t max_i = textLength - params.length + 1;
//     const uint64_t step_size = params.length - params.overlap + 1;
//
//     uint64_t progress_count;
//     uint64_t progress_max;
//     uint64_t progress_step;
//     initProgress<SearchParams::outputProgress>(progress_count, progress_step, progress_max, step_size, max_i);
//
//     #pragma omp parallel for schedule(dynamic, std::max(1ul, max_i/(step_size*params.threads*50))) num_threads(params.threads)
//     for (uint64_t i = 0; i < max_i; i += step_size)
//     {
//         uint64_t max_pos = std::min(i + params.length - params.overlap, textLength - params.length) + 1;
//
//         if (std::any_of(c.begin() + i, c.begin() + max_pos, [](auto value){ return value == 0; }))
//         {
//             TIter it_zero_errors[params.length - params.overlap + 1];
//             value_type hits[params.length - params.overlap + 1] = {};
//
//             auto delegate = [&hits, &it_zero_errors, i, textLength, params, &text](auto it, auto const & /*read*/, unsigned const errors_spent) {
//                 uint64_t const bb = std::min(textLength - 1, i + params.length - 1 + params.length - params.overlap);
//                 if (errors_spent == 0)
//                 {
//                     extend3<errors>(it, hits, it_zero_errors, errors - errors_spent, text, params.length,
//                         i + params.length - params.overlap, i + params.length - 1, // searched interval
//                         i, bb // entire interval
//                     );
//                 }
//                 else
//                 {
//                     extend(it, hits, errors - errors_spent, text, params.length,
//                         i + params.length - params.overlap, i + params.length - 1, // searched interval
//                         i, bb // entire interval
//                     );
//                 }
//             };
//
//             auto const & needle = infix(text, i + params.length - params.overlap, i + params.length);
//             TIter it(index);
//             _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
//             for (uint64_t j = i; j < max_pos; ++j)
//             {
//                 if (countOccurrences(it_zero_errors[j - i]) > 1) // guaranteed to exist, since there has to be at least one match!
//                 {
//                     for (auto const & occ : getOccurrences(it_zero_errors[j-i], Fwd()))
//                     {
//                         auto const occ_pos = posGlobalize(occ, limits);
//                         c[occ_pos] = hits[j - i];
//                     }
//                 }
//                 else
//                 {
//                     c[j] = hits[j - i];
//                 }
//             }
//         }
//         printProgress<SearchParams::outputProgress>(progress_count, progress_step, progress_max);
//     }
//
//     resetLimits(indexText(index), c, params.length);
// }
