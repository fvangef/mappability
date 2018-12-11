#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace std;
using namespace seqan;

int main(int argc, char *argv[])
{
    // Argument Parser
    // TODO: allow more than 4 gigabases
    ArgumentParser parser("Mappability Dumper");
    addDescription(parser, "Transforms the output of the mappability program into human readable format. Please be aware that the file size will increase significantly and most likely multiply.");

    addOption(parser, ArgParseOption("I", "input", "Path to the mappability file", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "input", "fa fasta");
	setRequired(parser, "input");

    addOption(parser, ArgParseOption("O", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString inputPath, outputPath;
    getOptionValue(inputPath, parser, "input");
    getOptionValue(outputPath, parser, "output");

    // Read fasta input file
    StringSet<CharString> ids;
    StringSet<Dna5String> chromosomes;
    SeqFileIn seqFileIn(toCString(inputPath));
    readRecords(ids, chromosomes, seqFileIn);
    close(seqFileIn);

    StringSet<Dna5String> chromosomes2;
    resize(chromosomes2, length(chromosomes));

    for (unsigned chr = 0; chr < length(chromosomes); ++chr)
    {
        cout << chr << ';';
        unsigned i = 0;
        while (i < length(chromosomes[chr]))
        {
            appendValue(chromosomes2[chr], chromosomes[chr][i]);
            ++i;

            // keep the first N, throw away all consecutive Ns
            if (chromosomes[chr][i - 1] == Dna5('N'))
            {
                unsigned j = 0;
                while (i < length(chromosomes[chr]) && chromosomes[chr][i] == Dna5('N'))
                {
                    ++j;
                    ++i;
                }
                cout << j << ';';
            }
        }
        cout << '\n';
    }

    SeqFileOut seqFileOut(toCString(outputPath));
    writeRecords(seqFileOut, ids, chromosomes2);
    close(seqFileOut);

    return 0;
}
