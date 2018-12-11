#include <vector>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.h"

using namespace std;
using namespace seqan;

template <typename T>
inline void save(vector<T> const & c, string const & output_path)
{
    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename T, typename TText>
inline void insertForward(vector<T> v, TText const & chromosomes, uint64_t const K)
{
    uint64_t total_length = 0;
    for (unsigned chr = 0; chr < length(chromosomes); ++chr)
        total_length += length(chromosomes[chr]);

    unsigned i_global = 0;
    for (unsigned chr = 0; chr < length(chromosomes); ++chr)
    {
        unsigned i = 0;
        while (i < length(chromosomes[chr]))
        {
            // appendValue(chromosomes2[chr], chromosomes[chr][i]);
            ++i;
            ++i_global;

            // keep the first N, throw away all consecutive Ns
            if (chromosomes[chr][i - 1] == Dna5('N'))
            {
                unsigned j = 0;
                while (i < length(chromosomes[chr]) && chromosomes[chr][i] == Dna5('N'))
                {
                    ++j;
                    ++i;

                    // cChromosomes.insert(cChromosomes.begin() + i_global, Dna5('N'));
                    if (v.size() >= i_global)
                        v.insert(v.begin() + i_global, 0);

                    ++i_global;
                }
            }
        }
    }

    cout << "Append 0 for missing trailing Ns: " << (total_length - K + 1 - v.size()) << "x \n";
    while (v.size() < total_length - K + 1)
        v.push_back(0);
}

template <typename T, typename TText>
inline void insertForward2(vector<T> v_old, vector<T> & v2, TText const & chromosomes, uint64_t const K)
{
    uint64_t total_length = 0;
    for (unsigned chr = 0; chr < length(chromosomes); ++chr)
        total_length += length(chromosomes[chr]);

    v2.reserve(total_length - K + 1);

    uint64_t progress_count;
    uint64_t progress_max;
    uint64_t progress_step;
    initProgress<true>(progress_count, progress_step, progress_max, 1, v_old.size());

    unsigned i_global = 0;
    for (unsigned chr = 0; chr < length(chromosomes); ++chr)
    {
        unsigned i = 0;
        while (i < length(chromosomes[chr]) && i_global < v_old.size())
        {
            v2.push_back(v_old[i_global]);
            // appendValue(chromosomes2[chr], chromosomes[chr][i]);
            ++i;
            ++i_global;

            // keep the first N, throw away all consecutive Ns
            if (chromosomes[chr][i - 1] == Dna5('N'))
            {
                unsigned j = 0;
                while (i < length(chromosomes[chr]) && chromosomes[chr][i] == Dna5('N') && i_global < v_old.size())
                {
                    // cChromosomes.insert(cChromosomes.begin() + i_global, Dna5('N'));
                    v2.push_back(0);

                    ++j;
                    ++i;
                }
            }
            printProgress<true>(progress_count, progress_step, progress_max);
        }
    }

    cout << "Append 0 for missing trailing Ns: " << (total_length - K + 1 - v2.size()) << "x \n";
    while (v2.size() < total_length - K + 1)
        v2.push_back(0);
}

template <typename value_t, typename TText>
void dump(CharString const & inputPath, CharString const & outputPath, uint64_t const K, TText const & chromosomes, TText & compressedChromosomes)
{
    vector<value_t> v;
    ifstream file(toCString(inputPath), std::ios::binary);
    if (!file.eof() && !file.fail())
    {
        file.seekg(0, std::ios_base::end);
        std::streampos fileSize = file.tellg();
        v.resize(fileSize / sizeof(value_t));
        file.seekg(0, std::ios_base::beg);
        file.read((char*) & v[0], fileSize);
        file.close();

        cout << "Load successful\n";

        vector<value_t> v2;
        insertForward2(v, v2, chromosomes, K);

        std::cout << "Size: " << v.size() << " -> " << v2.size() << '\n';
        std::string out = toCString(outputPath);
        save(v2, out);
        return;
    }

    cout << "Something went wrong ...\n";
}

int main(int argc, char *argv[])
{
    // Argument Parser
    // TODO: allow more than 4 gigabases
    ArgumentParser parser("Mappability Dumper");
    addDescription(parser, "Transforms the output of the mappability program into human readable format. Please be aware that the file size will increase significantly and most likely multiply.");

    addOption(parser, ArgParseOption("I", "input", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "input", "gmapp8 gmapp16");
	setRequired(parser, "input");

    addOption(parser, ArgParseOption("G", "genome", "Path to the ORIGINAL genome file", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "genome", "fa fasta");
	setRequired(parser, "genome");

    // addOption(parser, ArgParseOption("C", "compressedgenome", "Path to the COMPRESSED genome file", ArgParseArgument::INPUT_FILE, "IN"));
	// setValidValues(parser, "compressedgenome", "fa fasta");
	// setRequired(parser, "compressedgenome");

    addOption(parser, ArgParseOption("O", "output", "Path to mappability output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");


    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    uint64_t K;
    CharString genomePath, compressedGenomePath, inputPath, outputPath;
    getOptionValue(genomePath, parser, "genome");
    // getOptionValue(compressedGenomePath, parser, "compressedgenome");
    getOptionValue(inputPath, parser, "input");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(K, parser, "length");

    StringSet<Dna5String> chromosomes, compressedChromosomes;
    {
        StringSet<CharString> ids;
        SeqFileIn seqFileIn(toCString(genomePath));
        readRecords(ids, chromosomes, seqFileIn);
        close(seqFileIn);
    }
    // {
    //     StringSet<CharString> ids;
    //     SeqFileIn seqFileIn(toCString(compressedGenomePath));
    //     readRecords(ids, compressedChromosomes, seqFileIn);
    //     close(seqFileIn);
    // }

    string _inputPath = toCString(inputPath);

    if (_inputPath.substr(_inputPath.find_last_of(".") + 1) == "gmapp8")
        dump<uint8_t>(inputPath, outputPath, K, chromosomes/*, compressedChromosomes*/);
    else
        dump<uint16_t>(inputPath, outputPath, K, chromosomes/*, compressedChromosomes*/);

    return 0;
}
