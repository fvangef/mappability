#include <vector>
#include <cstdint>
#include <limits>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

struct Options
{
    unsigned errors;
    bool mmap;
    bool indels;
    bool high;
    CharString indexPath;
    CharString outputPath;
    CharString alphabet;
};

namespace seqan {
    template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TSize/*, typename TChar*/>
    inline bool
    _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it,
                   typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const &vDesc,
                   Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> &_range,
                   TSize &smaller,
                   Dna5 c) {
        typedef Index<TText, FMIndex<TOccSpec, TIndexSpec/*FMIndexConfig<TOccSpec, TLengthSum> */> > TIndex;
        typedef typename Fibre<TIndex, FibreLF>::Type TLF;

        if (c == Dna5('N'))
            return false;

        TIndex const &index = container(it);
        TLF const &lf = indexLF(index);

        _range = range(index, vDesc);

        TSize _smaller;
        _range.i1 = lf(_range.i1, c, _smaller);
        _range.i2 = lf(_range.i2, c, smaller);
        smaller -= _smaller;

        return _range.i1 < _range.i2;
    }
};

#include "common.h"
#include "algo2.hpp"
#include "algo3.hpp"
#include "algo4.hpp"

string get_output_path(Options const & opt, SearchParams const & searchParams)
{
    string output_path = toCString(opt.outputPath);
    output_path += "_" + to_string(opt.errors) + "_" + to_string(searchParams.length);
    output_path += ".gmapp" + string(opt.high ? "16" : "8");
    return output_path;
}

template <typename T>
inline void save(vector<T> const & c, string const & output_path)
{
    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename TDistance, typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams)
{
    vector<value_type> c(length(text) - searchParams.length + 1, 0);

    switch (opt.errors)
    {
        case 0:  runAlgo4<0>(index, text, c, searchParams);
                 break;
        case 1:  runAlgo4<1>(index, text, c, searchParams);
                 break;
        case 2:  runAlgo4<2>(index, text, c, searchParams);
                 break;
        case 3:  runAlgo4<3>(index, text, c, searchParams);
                 break;
        case 4:  runAlgo4<4>(index, text, c, searchParams);
                 break;
        default: cerr << "E = " << opt.errors << " not yet supported.\n";
                 exit(1);
    }

    if (SearchParams::outputProgress)
        std::cout << '\r';
    std::cout << "Progress: 100.00%\n" << std::flush;
    cout.flush();

    string output_path = get_output_path(opt, searchParams);
    save(c, output_path);
}

template <typename TChar, typename TAllocConfig, typename TDistance, typename value_type>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    typedef String<TChar, TAllocConfig> TString;

    typedef StringSet<TString, Owner<ConcatDirect<> > > TStringSet;
    TIndex<TStringSet> index;
    open(index, toCString(opt.indexPath), OPEN_RDONLY);
    auto const & text = indexText(index);
    run<TDistance, value_type>(index, text.concat, opt, searchParams);
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.high) {
        run<TChar, TAllocConfig, TDistance, uint16_t>(opt, searchParams);
    }
    else
        run<TChar, TAllocConfig, TDistance, uint8_t>(opt, searchParams);
}

template <typename TChar, typename TAllocConfig>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.indels) {
        run<TChar, TAllocConfig, EditDistance>(opt, searchParams);
    }
    else
        run<TChar, TAllocConfig, HammingDistance>(opt, searchParams);
}

template <typename TChar>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.mmap)
        run<TChar, MMap<> >(opt, searchParams);
    else
        run<TChar, Alloc<> >(opt, searchParams);
}

int main(int argc, char *argv[])
{
    // Argument parser
    ArgumentParser parser("Mappabilty");
    addDescription(parser,
        "App for calculating the mappability values. Only supports Dna4/Dna5 so far.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    // mention that file name will be prefix?
    addOption(parser, ArgParseOption("O", "output", "Path to output directory (error number, length and overlap will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("hi", "high", "Stores the mappability vector in 16 bit unsigned integers instead of 8 bit (max. value 65535 instead of 255)"));

    addOption(parser, ArgParseOption("o", "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region)", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "overlap");

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    Options opt;
    SearchParams searchParams;
    getOptionValue(opt.errors, parser, "errors");
    getOptionValue(opt.indexPath, parser, "index");
    getOptionValue(opt.outputPath, parser, "output");
    opt.mmap = isSet(parser, "mmap");
    opt.indels = isSet(parser, "indels");
    opt.high = isSet(parser, "high");

    getOptionValue(searchParams.length, parser, "length");
    getOptionValue(searchParams.threads, parser, "threads");
    getOptionValue(searchParams.overlap, parser, "overlap");

    if (searchParams.overlap > searchParams.length - 1)
    {
        cerr << "ERROR: overlap cannot be larger than K - 1.\n";
        exit(1);
    }

    if (!(searchParams.length - searchParams.overlap >= opt.errors + 2))
    {
        cerr << "ERROR: overlap should be at least K - E - 2. (K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts).\n";
        exit(1);
    }

    // searchParams.overlap - length of common overlap
    searchParams.overlap = searchParams.length - searchParams.overlap;

    if (opt.indels)
    {
        cerr << "ERROR: Indels are not supported yet.\n";
        exit(1);
    }

    CharString _indexPath = opt.indexPath;
    _indexPath += ".alphabet";
    open(opt.alphabet, toCString(_indexPath));

    if (opt.alphabet == "dna4")
    {
        run<Dna>(opt, searchParams);
    }
    else
    {
        cerr << "NOTE: k-mers with N bases are not considered and get a mappability value of 0.\n";
        run<Dna5>(opt, searchParams);
    }
}
