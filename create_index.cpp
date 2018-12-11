#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.h"

using namespace seqan;

template <typename TString, typename TStringSetConfig>
void buildIndex(StringSet<TString, TStringSetConfig> const & chromosomes, CharString const & indexPath)
{
    // TODO(cpockrandt): avoid copying
    StringSet<TString, Owner<ConcatDirect<> > > chromosomesConcat(chromosomes);
    TIndex<StringSet<TString, Owner<ConcatDirect<> > > > index(chromosomesConcat);
    // TODO: create fwd and rev separately and save them right afterwards to reduce memory peak
    indexCreate(index, FibreSALF());
    // do not store 
    clear(getFibre(getFibre(getFibre(index.rev, FibreSA()), FibreSparseString()), FibreValues()));
    clear(getFibre(getFibre(getFibre(index.rev, FibreSA()), FibreSparseString()), FibreIndicators()));
    save(index, toCString(indexPath));
}

int main(int argc, char *argv[])
{
    // Argument Parser
    // TODO: allow more than 4 gigabases
    ArgumentParser parser("Index Creation");
    addDescription(parser, "App for creating an index. Only supports Dna (with and without N's). At most 4 gigabases in total allowed.");

    addOption(parser, ArgParseOption("G", "genome", "Path to the genome", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "genome", "fa fasta fastq");
	setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "index");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, genomePath;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(genomePath, parser, "genome");

    // Read fasta input file
    StringSet<CharString> ids;
    StringSet<Dna5String> chromosomes;
    SeqFileIn seqFileIn(toCString(genomePath));
    readRecords(ids, chromosomes, seqFileIn);
    std::cout << "Number of sequences: " << length(chromosomes) << '\n' << std::flush;
    std::cout << "Sampling rate: " << TMyFastConfig::SAMPLING << '\n';

    // check whether it can be converted to Dna4
    unsigned i = 0;
    for (; i < length(chromosomes); ++i)
    {
        unsigned j = 0;
        for (; j < length(chromosomes[i]); ++j)
        {
            if (chromosomes[i][j] == 'N')
                break;
        }
        if (j < length(chromosomes[i]) && chromosomes[i][j] == 'N')
            break;
    }
    bool canConvert = i == length(chromosomes);
    std::cout << "Index will be constructed using Dna" << (5 - canConvert) << " alphabet.\n" << std::flush;

    if (canConvert)
    {
        StringSet<DnaString> chromosomes4(chromosomes);
        clear(chromosomes);
        buildIndex(chromosomes4, indexPath);
    }
    else
    {
        buildIndex(chromosomes, indexPath);
    }
    CharString _indexPath = indexPath;
    _indexPath += ".alphabet";
    CharString alphabet = canConvert ? "dna4" : "dna5";
    save(alphabet, toCString(_indexPath));

    // _indexPath = indexPath;
    // _indexPath += ".ids";
    // save(ids, toCString(_indexPath));
    std::cout << "Index created successfully.\n";

    return 0;
}
