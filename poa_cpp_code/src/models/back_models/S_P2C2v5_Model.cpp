// Author: Lance Hepler

#include "../UnanimityInternalConfig.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/data/Read.h>
#include <pacbio/data/internal/BaseEncoding.h>

#include "../ModelFactory.h"
#include "../Recursor.h"
#include "../Simulator.h"
#include "CounterWeight.h"
#include "HelperFunctions.h"

using namespace PacBio::Data;

namespace PacBio {
namespace Consensus {
namespace S_P2C2v5 {
namespace {

static constexpr const size_t CONTEXT_NUMBER = 16;
static constexpr const size_t OUTCOME_NUMBER = 12;

class S_P2C2v5_Model : public ModelConfig
{
public:
    static std::set<std::string> Chemistries() { return {"S/P2-C2/5.0"}; }
    static ModelForm Form() { return ModelForm::PWSNR; }
    S_P2C2v5_Model(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(const MappedRead& mr,
                                                     double scoreDiff) const override;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const override;
    std::pair<Data::Read, std::vector<MoveType>> SimulateRead(
        std::default_random_engine* const rng, const std::string& tpl,
        const std::string& readname) const override;

    double ExpectedLLForEmission(MoveType move, const AlleleRep& prev, const AlleleRep& curr,
                                 MomentType moment) const override;

private:
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
};

// TODO(lhepler) comments regarding the CRTP
class S_P2C2v5_Recursor : public Recursor<S_P2C2v5_Recursor>
{
public:
    S_P2C2v5_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    inline double EmissionPr(MoveType move, uint8_t emission, const AlleleRep& prev,
                             const AlleleRep& curr) const;
    double UndoCounterWeights(size_t nEmissions) const override;

private:
    double counterWeight_;
    double nLgCounterWeight_;
};

static constexpr const double snrRanges[2][4] = {
    {3.91070819, 7.37540293, 3.84641552, 6.27542830},  // minimum
    {8.06702614, 15.2471046, 7.02776623, 11.3469543}   // maximum
};

static constexpr const double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
     {0.3743876,0.0039716447,0.003089646,0.0035468806,0.3311385,0.0026748218,0.0020356874,0.0028444624,0.2698989,0.0023683421,0.0016301404,0.002413403},
{0.009436809,0.3246278,0.0021860588,0.0037720155,0.008225772,0.3434372,0.0015669978,0.0031456032,0.008247825,0.29069397,0.0014768471,0.0031831337},
{0.008979478,0.0023104819,0.42210242,0.0037093742,0.0076876422,0.0014465126,0.31013447,0.0027888336,0.006650736,0.0016090664,0.22982222,0.0027587574},
{0.008542142,0.0052987877,0.0027468759,0.38146475,0.007191374,0.0048200684,0.0021234446,0.3241355,0.0071702152,0.005069819,0.0023180307,0.24911898},
{0.3273029,0.0065320446,0.002292997,0.0031320942,0.32387394,0.006275928,0.0015626487,0.0025058608,0.31587377,0.0068165837,0.0013958178,0.002435421},
{0.0069536483,0.2746294,0.0052677374,0.0031375694,0.005536404,0.30923614,0.0037973074,0.0024950383,0.0052137696,0.37774998,0.003288132,0.002694876},
{0.004054551,0.0065068696,0.36172616,0.0033777722,0.0019943737,0.0062099895,0.31910747,0.002586603,0.0017548806,0.0069950377,0.2833174,0.0023689168},
{0.0028589051,0.0050079883,0.0025377779,0.31664756,0.0020887807,0.004832167,0.0019213151,0.31924525,0.0022057525,0.0057392432,0.0019954215,0.33491984},
{0.340009,0.0021713744,0.00978411,0.003618197,0.32017925,0.0012784195,0.008279505,0.002919115,0.29963332,0.0011985896,0.0081472155,0.0027818833},
{0.00289759,0.26747605,0.0060873865,0.0028682787,0.002314083,0.30967647,0.0050515183,0.002122804,0.0024591289,0.39182794,0.0049466626,0.0022720804},
{0.005159074,0.0045489594,0.35153314,0.005326447,0.0036750333,0.0033829275,0.3114824,0.0042927195,0.0034889304,0.0041420837,0.29849863,0.004469657},
{0.0037127328,0.002214619,0.0071076727,0.31888786,0.002911998,0.0014363208,0.005984731,0.31670412,0.0032599517,0.0019523664,0.0060034897,0.32982412},
{0.26705155,0.0019850007,0.001508162,0.0088831745,0.3013104,0.0014200108,0.0010199049,0.008568936,0.396933,0.0013008012,0.0009577879,0.009061304},
{0.0020560948,0.23314115,0.0017548224,0.00733474,0.0017420627,0.2890201,0.0013185093,0.0068764514,0.00195189,0.44576177,0.0012922809,0.007750141},
{0.0019294656,0.0018505667,0.3000879,0.007320583,0.0014102797,0.0012923954,0.30770272,0.006955482,0.0014802058,0.0014687577,0.36112967,0.007371945},
{0.0029792555,0.0024499511,0.002885215,0.29264787,0.0027836643,0.0020375443,0.0024414614,0.31593922,0.0035728863,0.002806848,0.0029661944,0.36648992}},
    {// branchPmf
{0.29235402,0.019638648,0.021733439,0.04084839,0.2627651,0.01361613,0.017282011,0.040455617,0.2147159,0.015056297,0.016889239,0.044645194},
{0.34281173,0,0.0011725841,0.004137869,0.316282,0,0.0009245374,0.0039461963,0.32548228,0,0.0007779644,0.0044648396},
{0.37170875,0.0016672092,0,0.0017559956,0.3209131,0.0015981533,0,0.0013219293,0.2981641,0.0014896367,0,0.0013811201},
{0.36575347,0.0021407702,0.0025168513,0,0.30673763,0.0017791535,0.0018514768,0,0.31540197,0.002010588,0.0018080828,0},
{0,0.32128757,0.0013652969,0.0019187956,0,0.3062939,0.00088559795,0.0018203958,0,0.36353797,0.0011807972,0.001709696},
{0.06865305,0.16495162,0.041161112,0.070803255,0.05099063,0.16157272,0.028259868,0.061127324,0.049147595,0.18722163,0.030870834,0.085240364},
{0.0052858666,0.3113058,0,0.0035674928,0.00414651,0.29623264,0,0.0031378996,0.003941052,0.369469,0,0.0029137637},
{0.0061308085,0.31980577,0.003743591,0,0.003960611,0.2817459,0.0022515801,0,0.003906356,0.3761224,0.0023329626,0},
{0,0.0018249375,0.35798368,0.0018578193,0,0.0017262923,0.30552083,0.0017920558,0,0.0012330659,0.32623637,0.0018249375},
{0.0064270776,0,0.37029853,0.0059517017,0.0056474614,0,0.29712874,0.0046967105,0.00427838,0,0.3002662,0.005305191},
{0.04367774,0.03743806,0.22499542,0.06735181,0.024958707,0.030831344,0.19434759,0.06074509,0.026426867,0.03468526,0.18186823,0.07267389},
{0.007172532,0.0018071055,0.36663747,0,0.0053467965,0.0017139557,0.3047674,0,0.0042289994,0.0019188852,0.30640686,0},
{0,0.00057767646,0.0008087471,0.31796363,0,0.00042012834,0.00076673424,0.31970719,0,0.00043063157,0.00051465724,0.3588106},
{0.0050365217,0,0.0018226086,0.3079513,0.0034365216,0,0.0011826087,0.30864695,0.0023791303,0,0.0010991305,0.36844522},
{0.0032929129,0.00086964166,0,0.32318425,0.0021203624,0.00081101415,0,0.31332505,0.002149676,0.00081101415,0,0.35343605},
{0.029439924,0.013523217,0.027884156,0.2438966,0.023336524,0.009454285,0.017711824,0.25957397,0.024413595,0.010292006,0.019387268,0.32108665}},
    {// stickPmf
 {0.050029963,0.18459284,0.12791798,0.09762301,0.045821957,0.077661626,0.076290034,0.08170983,0.04495639,0.070630535,0.06723484,0.07553099},
{0.026884483,0.19798343,0.07507877,0.07417713,0.02541083,0.18040623,0.047787096,0.059595715,0.027136555,0.17844781,0.047341123,0.059750836},
{0.029839834,0.087600246,0.20140745,0.13129756,0.026823863,0.04910092,0.1468001,0.06275276,0.026058447,0.04896383,0.1268193,0.0625357},
{0.020332914,0.074926294,0.06729735,0.23021524,0.01838667,0.04231644,0.042899493,0.2017853,0.018337399,0.03715109,0.04226717,0.20408465},
{0.26381382,0.012374412,0.052785765,0.03705558,0.24013396,0.012685633,0.028557897,0.029450966,0.25284666,0.014992727,0.026514664,0.02878793},
{0.12485387,0.031678602,0.14031447,0.10554046,0.09373508,0.03388491,0.09228616,0.08257183,0.09012925,0.03663456,0.0854532,0.08291759},
{0.14031652,0.014349621,0.20443882,0.0849042,0.034721527,0.015196932,0.1690522,0.062309243,0.031787828,0.018695505,0.15893002,0.0652976},
{0.07513966,0.009043322,0.05035353,0.22849199,0.03356902,0.009980858,0.027826607,0.22426006,0.034838602,0.011367631,0.027520606,0.2676081},
{0.25475436,0.081351586,0.01901749,0.041762117,0.21820116,0.032007925,0.015808515,0.03378059,0.22381005,0.031080687,0.016481219,0.031944294},
{0.077661835,0.19923002,0.019596912,0.097090445,0.052763343,0.17922285,0.01644121,0.043811668,0.049775947,0.20790818,0.0154629415,0.04103465},
{0.13717051,0.12721054,0.038805474,0.133667,0.07487488,0.08943944,0.031698365,0.08475142,0.067334004,0.1004004,0.030363698,0.08428428},
{0.08099952,0.06883961,0.014819255,0.22816554,0.030403242,0.026827605,0.0118223345,0.2181758,0.03171224,0.02660714,0.012084134,0.24954358},
{0.22975042,0.041377854,0.04056327,0.020123243,0.22230603,0.025463294,0.024543116,0.018380938,0.3110429,0.023645565,0.024067942,0.018735433},
{0.04505568,0.19280894,0.06797858,0.020171411,0.03329677,0.19468665,0.043466136,0.018191451,0.03499786,0.2860715,0.04327093,0.02000409},
{0.052890528,0.061646275,0.22660549,0.018895255,0.030564187,0.03283,0.2100165,0.016087266,0.029593125,0.03361494,0.26972875,0.017527675},
{0.11146439,0.12373928,0.13871153,0.044343386,0.072606914,0.07543457,0.09301295,0.041502696,0.0761252,0.08265357,0.0950718,0.045333717}}};

static constexpr const double transProbs[CONTEXT_NUMBER][4] = {
    {0.90453374,0.0022890975,0.022505865,0.07067132}, // AA
{0.8598327,0.029506432,0.03431433,0.07634652}, // AC
{0.8662158,0.031438626,0.02714837,0.075197205}, // AG
{0.85688895,0.022383133,0.03942577,0.081302114}, // AT
{0.8597026,0.021593679,0.03925725,0.079446495}, // CA
{0.8611872,0.0021928141,0.020454701,0.11616526}, // CC
{0.8594706,0.017604355,0.036090262,0.08683476}, // CG
{0.85997194,0.011515955,0.047982577,0.080529526}, // CT
{0.86066765,0.019900179,0.035990715,0.08344144}, // GA
{0.8888599,0.014125379,0.025534194,0.07148049}, // GC
{0.87803555,0.0019066241,0.020973215,0.099084616}, // GG
{0.8644595,0.01706724,0.04615185,0.072321385}, // GT
{0.84729385,0.036845375,0.05130891,0.06455187}, // TA
{0.8764729,0.022328349,0.033419676,0.06777908}, // TC
{0.8707742,0.027573135,0.033294357,0.068358295}, // TG
{0.905393,0.0024703476,0.022687819,0.069448814} // TT
    };

inline double CalculateExpectedLLForEmission(const size_t move, const uint8_t row,
                                             const size_t moment)
{
    double expectedLL = 0;
    for (size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[move][row][i];
        double lgCurProb = std::log(curProb);
        if (!std::isfinite(lgCurProb)) continue;
        if (moment == static_cast<uint8_t>(MomentType::FIRST))
            expectedLL += curProb * lgCurProb;
        else if (moment == static_cast<uint8_t>(MomentType::SECOND))
            expectedLL += curProb * (lgCurProb * lgCurProb);
    }
    return expectedLL;
}

S_P2C2v5_Model::S_P2C2v5_Model(const SNR& snr)
    : snr_(ClampSNR(snr, SNR{snrRanges[0]}, SNR{snrRanges[1]}))
{
    for (size_t ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
        const uint8_t bp = ctx & 3;  // (equivalent to % 4)
        const double snr1 = snr_[bp], snr2 = snr1 * snr1, snr3 = snr2 * snr1;
        double sum = 1.0;

        // cached transitions
        for (size_t j = 0; j < 4; ++j) {
            double xb = transProbs[ctx][j];
            ctxTrans_[ctx][j] = xb;
            sum += xb;
        }
        for (size_t j = 0; j < 4; ++j)
            ctxTrans_[ctx][j] /= sum;

        // cached expectations
        for (size_t move = 0; move < 3; ++move)
            for (size_t moment = 0; moment < 2; ++moment)
                cachedEmissionExpectations_[ctx][move][moment] =
                    CalculateExpectedLLForEmission(move, ctx, moment);
    }
}

std::unique_ptr<AbstractRecursor> S_P2C2v5_Model::CreateRecursor(const MappedRead& mr,
                                                                 double scoreDiff) const
{
    const double counterWeight = CounterWeight(
        [this](size_t ctx, MoveType m) { return ctxTrans_[ctx][static_cast<uint8_t>(m)]; },
        [](size_t ctx, MoveType m) {
            double r = 0.0;
            for (size_t o = 0; o < OUTCOME_NUMBER; ++o) {
                const double p = emissionPmf[static_cast<uint8_t>(m)][ctx][o];
                if (p > 0.0) r += p * std::log(p);
            }
            return r;
        },
        CONTEXT_NUMBER);

    return std::make_unique<S_P2C2v5_Recursor>(mr, scoreDiff, counterWeight);
}

std::vector<TemplatePosition> S_P2C2v5_Model::Populate(const std::string& tpl) const
{
    auto rowFetcher = [this](const NCBI2na prev, const NCBI2na curr) -> const double(&)[4]
    {
        const auto row = EncodeContext16(prev, curr);
        const double(&params)[4] = ctxTrans_[row];
        return params;
    };
    return AbstractPopulater(tpl, rowFetcher);
}

double S_P2C2v5_Model::ExpectedLLForEmission(const MoveType move, const AlleleRep& prev,
                                             const AlleleRep& curr, const MomentType moment) const
{
    auto cachedEmissionVisitor = [this](const MoveType move, const NCBI2na prev, const NCBI2na curr,
                                        const MomentType moment) -> double {
        const auto row = EncodeContext16(prev, curr);
        return cachedEmissionExpectations_[row][static_cast<uint8_t>(move)]
                                          [static_cast<uint8_t>(moment)];
    };
    return AbstractExpectedLLForEmission(move, prev, curr, moment, cachedEmissionVisitor);
}

S_P2C2v5_Recursor::S_P2C2v5_Recursor(const MappedRead& mr, double scoreDiff, double counterWeight)
    : Recursor<S_P2C2v5_Recursor>(mr, scoreDiff)
    , counterWeight_{counterWeight}
    , nLgCounterWeight_{-std::log(counterWeight_)}
{
}

std::vector<uint8_t> S_P2C2v5_Recursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;
    result.reserve(read.Length());
    std::vector<uint8_t> boundaries = {18, 46};
    for (size_t i = 0; i < read.Length(); ++i) {
        result.emplace_back(EncodeBaseV2(read.Seq[i], read.PulseWidth[i], boundaries));
    }

    return result;
}

double S_P2C2v5_Recursor::EmissionPr(const MoveType move, const uint8_t emission,
                                     const AlleleRep& prev, const AlleleRep& curr) const
{
    return AbstractEmissionPr(emissionPmf, move, emission, prev, curr) * counterWeight_;
}

double S_P2C2v5_Recursor::UndoCounterWeights(const size_t nEmissions) const
{
    return nLgCounterWeight_ * nEmissions;
}

inline std::pair<Data::SNR, std::vector<TemplatePosition>> S_P2C2v5_InitialiseModel(
    std::default_random_engine* const rng, const std::string& tpl)
{
    Data::SNR snrs{0, 0, 0, 0};
    for (uint8_t i = 0; i < 4; ++i) {
        snrs[i] = std::uniform_real_distribution<double>{snrRanges[0][i], snrRanges[1][i]}(*rng);
    }

    const S_P2C2v5_Model model{snrs};
    std::vector<TemplatePosition> transModel = model.Populate(tpl);

    return {snrs, transModel};
}

BaseData S_P2C2v5_GenerateReadData(std::default_random_engine* const rng, const MoveType state,
                                   const AlleleRep& prev, const AlleleRep& curr)
{
    // distribution is arbitrary at the moment, as
    // IPD is not a covariate of the consensus HMM
    std::uniform_int_distribution<uint8_t> ipdDistrib(1, 5);

    std::array<double, OUTCOME_NUMBER> emissionDist;
    for (size_t i = 0; i < OUTCOME_NUMBER; ++i) {
        emissionDist[i] = AbstractEmissionPr(emissionPmf, state, i, prev, curr);
    }

    std::discrete_distribution<uint8_t> outcomeDistrib(emissionDist.cbegin(), emissionDist.cend());

    const uint8_t event = outcomeDistrib(*rng);
    const std::pair<char, uint8_t> outcome = DecodeEmission(event);

    return {outcome.first, outcome.second, ipdDistrib(*rng)};
}

std::pair<Data::Read, std::vector<MoveType>> S_P2C2v5_Model::SimulateRead(
    std::default_random_engine* const rng, const std::string& tpl,
    const std::string& readname) const
{
    return SimulateReadImpl(rng, tpl, readname, S_P2C2v5_InitialiseModel,
                            S_P2C2v5_GenerateReadData);
}

}  // namespace anonymous
}  // namespace S_P2C2v5

REGISTER_MODEL_IMPL(S_P2C2v5)

}  // namespace Consensus
}  // namespace PacBio
