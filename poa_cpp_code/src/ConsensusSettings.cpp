// Author: Armin Töpfer

#include <pacbio/ccs/ConsensusSettings.h>

namespace PacBio {
namespace CCS {
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption MaxLength{
    "max_length",
    { "maxLength" },
    "maxLength",
    "Maximum length of subreads to use for generating SMC.",
    CLI::Option::IntType(21000)
};
const PlainOption MinLength{
    "min_length",
    { "minLength" },
    "minLength",
    "Minimum length of subreads to use for generating SMC.",
    CLI::Option::IntType(10)
};
const PlainOption MinPasses{
    "min_passes",
    { "minPasses" },
    "minPasses",
    "Minimum number of subreads required to generating SMC.",
    CLI::Option::IntType(3)
};
const PlainOption MinPredictedAccuracy{
    "min_predicted_accuracy",
    { "minPredictedAccuracy" },
    "minPredictedAccuracy",
    "Minimum predicted accuracy in [0, 1].",
    CLI::Option::FloatType(0.9)
};
const PlainOption MinIdentity{
    "min_identity",
    { "minIdentity" },
    "minIdentity",
    "Minimum identity to the POA to use a subread. 0 disables this filter.",
    CLI::Option::FloatType(0.82)
};
const PlainOption MinZScore{
    "min_zscore",
    { "minZScore" },
    "minZScore",
    "Minimum z-score to use a subread. NaN disables this filter.",
    CLI::Option::FloatType(-3.4)
};
const PlainOption MaxDropFraction{
    "max_drop_fraction",
    { "maxDropFraction" },
    "maxDropFraction",
    "Maximum fraction of subreads that can be dropped before giving up.",
    CLI::Option::FloatType(0.34)
};
const PlainOption NoPolish{
    "no_polish",
    { "noPolish" },
    "noPolish",
    "Only output the initial template derived from the POA (faster, less accurate).",
    CLI::Option::BoolType(false)
};
const PlainOption Polish{
    "polish",
    { "polish" },
    "polish",
    "Emit high-accuracy SMC sequences polished using the Arrow algorithm",
    CLI::Option::BoolType(true)
};
const PlainOption PolishRepeats{
    "polish_repeats",
    { "polishRepeats" },
    "polishRepeats",
    "Polish repeats of 2 to N bases of 3 or more elements.",
    CLI::Option::IntType(0)
};
const PlainOption MinReadScore{
    "min_read_score",
    { "minReadScore" },
    "minReadScore",
    "Minimum read score of input subreads.",
    CLI::Option::FloatType(0.75)
};
const PlainOption MinSnr{
    "min_snr",
    { "minSnr" },
    "minSnr",
    "Minimum SNR of input subreads.",
    CLI::Option::FloatType(3.75)
    // See https://github.com/PacificBiosciences/pbccs/issues/86 for a more
    // detailed discussion of this default.)
};
const PlainOption ByStrand{
    "by_strand",
    { "byStrand" },
    "byStrand",
    "Generate a consensus for each strand.",
    CLI::Option::BoolType()
};
const PlainOption ForceOutput{
    "force",
    { "force" },
    "force",
    "Overwrite OUTPUT file if present.",
    CLI::Option::BoolType()
};
const PlainOption Zmws{
    "zmws",
    { "zmws" },
    "channels",
    "Generate SMC for the provided comma-separated channels ranges only. Default = all",
    CLI::Option::StringType(""),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};
const PlainOption ReportFile{
    "report_file",
    { "reportFile" },
    "reportFile",
    "Where to write the results report.",
    CLI::Option::StringType("smc_report.txt")
};
const PlainOption NumThreads{
    "num_threads",
    { "numThreads" },
    "numThreads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption LogFile{
    "log_file",
    { "logFile" },
    "logFile",
    "Log to a file, instead of STDERR.",
    CLI::Option::StringType("")
};
const PlainOption RichQVs{
    "rich_qvs",
    { "richQVs" },
    "richQVs",
    "Emit dq, iq, and sq \"rich\" quality tracks.",
    CLI::Option::BoolType()
};
const PlainOption ModelPath{
    "model_path",
    { "modelPath" },
    "modelPath",
    "Path to a model file or directory containing model files.",
    CLI::Option::StringType("")
};
const PlainOption ModelSpec{
    "model_spec",
    { "modelSpec" },
    "modelSpec",
    "Name of chemistry or model to use, overriding default selection.",
    CLI::Option::StringType("")
};
const PlainOption ZmwTimings{
    "zmw_timings",
    { "zmwTimings" },
    "channelTimings",
    "Measure individual ZMW wall clock timings.",
    CLI::Option::BoolType(),
    JSON::Json(nullptr),
    CLI::OptionFlags::HIDE_FROM_HELP
};

const PlainOption PoaErrSpanVersion{
    "poa_err_span_v",
    {"poaErrSpanV"},
    "poaErrSpanV",
    "err calculation version. CHOICES[v1/v2]",
    CLI::Option::StringType("v1")
};

const PlainOption PoaIdentityVersion{
    "poa_identity_v",
    {"poaIdentityV"},
    "poaIdentityV",
    "err identity calculation version. CHOICES[v1/v2]",
    CLI::Option::StringType("v1")
};

const PlainOption PassesCntVersion{
    "passes_cnt_version",
    {"passesCntV"},
    "passesCntV",
    "passes cnt version. CHOICES[v1/v2]",
    CLI::Option::StringType("v1")
};

const PlainOption ReComputePoaIdentityVersion{
    "re-compute poa identity",
    {"reComputePoaIdentityV"},
    "reComputePoaIdentityV",
    "re-compute poa identity. CHOICES[v1/v2]",
    CLI::Option::StringType("v1")
};

const PlainOption BamThreads{
    "bam_threads",
    { "bamThreads" },
    "bamThreads",
    "num threads that bam used for zip/unzip",
    CLI::Option::IntType(4)
};

const PlainOption BamCompressionLevel{
    "bam_compression_level",
    { "bamCompressionLevel" },
    "bamCompressionLevel",
    "bamCompressionLevel -1 as default, 1 for fast compression, 9 for best compression",
    CLI::Option::IntType(-1)
};

const PlainOption MatchScore{
    "match_score",
    { "matchScore" },
    "matchScore",
    "matchScore for alignment 3 as default",
    CLI::Option::IntType(3)
};

const PlainOption MismatchScore{
    "mismatch_score",
    { "mismatchScore" },
    "mismatchScore",
    "mismatchScore for alignment -5 as default",
    CLI::Option::IntType(-5)
};

const PlainOption InsertionScore{
    "insertion_score",
    { "insertionScore" },
    "insertionScore",
    "insertionScore for alignment -4 as default",
    CLI::Option::IntType(-4)
};

const PlainOption DeletionScore{
    "deletion_score",
    { "deletionScore" },
    "deletionScore",
    "deletionScore for alignment -4 as default",
    CLI::Option::IntType(-4)
};

// clang-format on
}  // namespace OptionNames

ConsensusSettings::ConsensusSettings(const PacBio::CLI::Results& options)
    : ByStrand{options[OptionNames::ByStrand]}
    , ForceOutput{options[OptionNames::ForceOutput]}
    , LogFile{options[OptionNames::LogFile].get<decltype(LogFile)>()}
    , LogLevel{options.LogLevel()}
    , MaxDropFraction{options[OptionNames::MaxDropFraction]}
    , MaxLength{options[OptionNames::MaxLength]}
    , MinLength{options[OptionNames::MinLength]}
    , MinPasses{options[OptionNames::MinPasses]}
    , MinPredictedAccuracy{options[OptionNames::MinPredictedAccuracy]}
    , MinReadScore{options[OptionNames::MinReadScore]}
    , MinSNR{options[OptionNames::MinSnr]}
    , MinIdentity{options[OptionNames::MinIdentity]}
    , MinZScore{options[OptionNames::MinZScore] == nullptr
                    ? NAN
                    : static_cast<float>(options[OptionNames::MinZScore])}
    , ModelPath{options[OptionNames::ModelPath].get<decltype(ModelPath)>()}
    , ModelSpec{options[OptionNames::ModelSpec].get<decltype(ModelSpec)>()}
    , PolishRepeats{options[OptionNames::PolishRepeats]}
    , ReportFile{options[OptionNames::ReportFile].get<decltype(ReportFile)>()}
    , RichQVs{options[OptionNames::RichQVs]}
    , WlSpec{options[OptionNames::Zmws].get<decltype(WlSpec)>()}
    , ZmwTimings{options[OptionNames::ZmwTimings]}
    , PoaErrSpanVersion{options[OptionNames::PoaErrSpanVersion].get<decltype(PoaErrSpanVersion)>()}
    , PoaIdentityVersion{options[OptionNames::PoaIdentityVersion].get<decltype(PoaIdentityVersion)>()}
    , PassesCntVersion{options[OptionNames::PassesCntVersion].get<decltype(PassesCntVersion)>()}
    , ReComputePoaIdentityVersion{options[OptionNames::ReComputePoaIdentityVersion].get<decltype(ReComputePoaIdentityVersion)>()}
    , BamThreads{options[OptionNames::BamThreads]}
    , BamCompressionLevel{options[OptionNames::BamCompressionLevel]}
    , MatchScore{options[OptionNames::MatchScore]}
    , MismatchScore{options[OptionNames::MismatchScore]}
    , InsertionScore{options[OptionNames::InsertionScore]}
    , DeletionScore{options[OptionNames::DeletionScore]}
{
    // N.B. If the user somehow specifies both polish and noPolish, noPolish wins.
    // Unfortunately there's no sensible way to check for this condition and error out.
    // This could be improved upon in the pbcopper API, perhaps.
    NoPolish = options[OptionNames::NoPolish] || !options[OptionNames::Polish];

    // N.B. This is the trick to resolved nthreads from either our
    // option or the "nproc" which has meaning in tool contracts.
    // Derek says he may streamline the API in the future.
    int requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NThreads = ThreadCount(requestedNThreads);
}

size_t ConsensusSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

PacBio::CLI::Interface ConsensusSettings::CreateCLI(const std::string& description,
                                                    const std::string& version)
{
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"smc", description, version};

    i.AlternativeToolContractName("gssmc");

    i.AddHelpOption();      // use built-in help output
    i.AddLogLevelOption();  // use built-in logLevel option
    i.AddVersionOption();   // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"input",  "Input file.",  "INPUT"},
        {"output", "Output file.", "OUTPUT"}
    });

    i.AddOptions(
    {
        OptionNames::ForceOutput,
        OptionNames::Zmws,
        OptionNames::MaxLength,
        OptionNames::MinLength,
        OptionNames::MinPasses,
        OptionNames::MinPredictedAccuracy,
        OptionNames::MinIdentity,
        OptionNames::MinZScore,
        OptionNames::MaxDropFraction,
        OptionNames::MinSnr,
        OptionNames::MinReadScore,
        OptionNames::ByStrand,
        OptionNames::NoPolish,
        OptionNames::Polish,
        OptionNames::PolishRepeats,
        OptionNames::RichQVs,
        OptionNames::ReportFile,
        OptionNames::ModelPath,
        OptionNames::ModelSpec,
        OptionNames::NumThreads,
        OptionNames::LogFile,
        OptionNames::ZmwTimings,
        OptionNames::PoaErrSpanVersion,
        OptionNames::PoaIdentityVersion,
        OptionNames::PassesCntVersion,
        OptionNames::ReComputePoaIdentityVersion,
        OptionNames::BamThreads,
        OptionNames::BamCompressionLevel,
        OptionNames::MatchScore,
        OptionNames::MismatchScore,
        OptionNames::InsertionScore,
        OptionNames::DeletionScore
    });

    const std::string id = "geneus.tasks.smc";
    Task tcTask(id);
    tcTask.AddOption(OptionNames::MinSnr);
    tcTask.AddOption(OptionNames::MinReadScore);
    tcTask.AddOption(OptionNames::MaxLength);
    tcTask.AddOption(OptionNames::MinLength);
    tcTask.AddOption(OptionNames::MinPasses);
    tcTask.AddOption(OptionNames::MinPredictedAccuracy);
    tcTask.AddOption(OptionNames::MinIdentity);
    tcTask.AddOption(OptionNames::MinZScore);
    tcTask.AddOption(OptionNames::MaxDropFraction);
    tcTask.AddOption(OptionNames::Polish);
    tcTask.AddOption(OptionNames::ByStrand);
    tcTask.AddOption(OptionNames::ModelPath);
    tcTask.AddOption(OptionNames::ModelSpec);
    tcTask.AddOption(OptionNames::ReportFile);
    tcTask.AddOption(OptionNames::RichQVs);
    tcTask.AddOption(OptionNames::PoaErrSpanVersion);
    tcTask.AddOption(OptionNames::PoaIdentityVersion);
    tcTask.AddOption(OptionNames::PassesCntVersion);
    tcTask.AddOption(OptionNames::ReComputePoaIdentityVersion);
    tcTask.AddOption(OptionNames::BamThreads);
    tcTask.AddOption(OptionNames::BamCompressionLevel);
    tcTask.AddOption(OptionNames::MatchScore);
    tcTask.AddOption(OptionNames::MismatchScore);
    tcTask.AddOption(OptionNames::InsertionScore);
    tcTask.AddOption(OptionNames::DeletionScore);
    tcTask.NumProcessors(Task::MAX_NPROC);

    tcTask.InputFileTypes({
        {
            "subread_set",
            "SubreadSet",
            "Subread DataSet or .bam file",
            "Geneus.DataSet.SubreadSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "bam_output",
            "Consensus Sequences",
            "Consensus sequences generated by SMC",
            "Geneus.DataSet.ConsensusReadSet",
            "smc"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}
}  // ::PacBio::CCS
