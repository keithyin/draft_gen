#include <gtest/gtest.h>
#include <pbcopper/cli/Interface.h>
#include <pbcopper/cli/toolcontract/JsonPrinter.h>
#include <pbcopper/utility/PbcopperVersion.h>

using namespace PacBio;
using namespace PacBio::CLI;
using namespace PacBio::CLI::ToolContract;

namespace ToolContractJsonPrinterTests {

static std::string MakeHeader(void)
{
    const auto preVersion = R"({
    "driver": {
        "env": {},
        "exe": "frobber --resolved-tool-contract",
        "serialization": "json"
    },
    "schema_version": "2.0.0",
    "tool_contract": {
        "_comment": "Created by v)";

    const auto postVersion = R"(",
        "description": "Frobb your files in a most delightful, nobbly way",
        "input_types": [
            {
                "description": "Generic Txt file",
                "file_type_id": "PacBio.FileTypes.txt",
                "id": "txt_in",
                "title": "Txt file"
            }
        ],
        "is_distributed": true,
        "name": "frobber",
        "nproc": )";

    return preVersion + PacBio::Utility::LibraryVersionString() + postVersion;
}

static std::string MakeToolContractText(const std::string& nProc)
{
    const auto footer = R"(,
        "output_types": [
            {
                "default_name": "output",
                "description": "Generic Output Txt file",
                "file_type_id": "PacBio.FileTypes.txt",
                "id": "txt_out",
                "title": "Txt outfile"
            }
        ],
        "resource_types": [
            "$tmpfile",
            "$tmpfile",
            "$tmpdir"
        ],
        "schema_options": [
            {
                "choices": [
                    1,
                    2,
                    3
                ],
                "default": 1,
                "description": "Some delta for things",
                "id": "frobber.task_options.delta",
                "name": "Frobbing delta",
                "optionTypeId": "choice_integer"
            },
            {
                "default": 10,
                "description": "Max Number of lines to Copy",
                "id": "frobber.task_options.max_nlines",
                "name": "Max lines",
                "optionTypeId": "integer"
            },
            {
                "choices": [
                    "normal",
                    "super",
                    "secret"
                ],
                "default": "normal",
                "description": "Run frobber in desired mode",
                "id": "frobber.task_options.mode",
                "name": "Frobbing mode",
                "optionTypeId": "choice_string"
            }
        ],
        "task_type": "pbsmrtpipe.task_types.standard",
        "tool_contract_id": "frobber.tasks.dev_txt_app"
    },
    "tool_contract_id": "frobber.tasks.dev_txt_app",
    "version": "3.14"
})";

    return MakeHeader() + nProc + footer;
}

}  // namespace ToolContractJsonPrinterTests

TEST(CLI_ToolContractJsonPrinter, prints_expected_output)
{
    const std::string expectedText = ToolContractJsonPrinterTests::MakeToolContractText("1");

    // --------------------------------------------
    // setup interface & tool contract config
    // --------------------------------------------

    // clang-format off
    Interface i {
        "frobber",
        "Frobb your files in a most delightful, nobbly way",
        "3.14"
    };
    i.AddOptions({
        {"max_nlines", "n", "Max Number of lines to Copy", Option::IntType(10)},
        {"mode", "mode", "Run frobber in desired mode", Option::StringType("normal"), { "normal", "super", "secret" }},
        {"delta", "d", "Some delta for things", Option::IntType(1), { 1, 2, 3 }}
    });

    const std::string id{"frobber.tasks.dev_txt_app"};
    ToolContract::Task tcTask(id);
    tcTask.InputFileTypes({
        { "txt_in", "Txt file", "Generic Txt file", "PacBio.FileTypes.txt" }
    });
    tcTask.OutputFileTypes({
        { "txt_out", "Txt outfile", "Generic Output Txt file", "PacBio.FileTypes.txt", "output" }
    });
    tcTask.ResourceTypes({
        ResourceType::TMP_FILE,
        ResourceType::TMP_FILE,
        ResourceType::TMP_DIR
    });
    tcTask.Options({
        { "max_nlines","Max lines" },
        { "mode", "Frobbing mode" },
        { "delta", "Frobbing delta" }
    });
    // clang-format on

    ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);

    // --------------------------------------------
    // test tool contract output
    // --------------------------------------------

    std::ostringstream s;
    ToolContract::JsonPrinter::Print(i, s, 4);
    EXPECT_EQ(expectedText, s.str());
}

TEST(CLI_ToolContractJsonPrinter, prints_expected_output_with_explicit_nproc)
{
    const std::string expectedText = ToolContractJsonPrinterTests::MakeToolContractText("4");

    // --------------------------------------------
    // setup interface & tool contract config
    // --------------------------------------------

    // clang-format off
    Interface i {
        "frobber",
        "Frobb your files in a most delightful, nobbly way",
        "3.14"
    };
    i.AddOptions({
        {"max_nlines", "n", "Max Number of lines to Copy", Option::IntType(10)},
        {"mode", "mode", "Run frobber in desired mode", Option::StringType("normal"), { "normal", "super", "secret" }},
        {"delta", "d", "Some delta for things", Option::IntType(1), { 1, 2, 3 }}
    });

    const std::string id{"frobber.tasks.dev_txt_app"};
    ToolContract::Task tcTask(id);
    tcTask.InputFileTypes({
        { "txt_in", "Txt file", "Generic Txt file", "PacBio.FileTypes.txt" }
    });
    tcTask.OutputFileTypes({
        { "txt_out", "Txt outfile", "Generic Output Txt file", "PacBio.FileTypes.txt", "output" }
    });
    tcTask.ResourceTypes({
        ResourceType::TMP_FILE,
        ResourceType::TMP_FILE,
        ResourceType::TMP_DIR
    });
    tcTask.Options({
        { "max_nlines","Max lines" },
        { "mode", "Frobbing mode" },
        { "delta", "Frobbing delta" }
    });
    tcTask.NumProcessors(4);
    // clang-format on

    ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);

    // --------------------------------------------
    // test tool contract output
    // --------------------------------------------

    std::ostringstream s;
    ToolContract::JsonPrinter::Print(i, s, 4);
    EXPECT_EQ(expectedText, s.str());
}

TEST(CLI_ToolContractJsonPrinter, prints_expected_output_with_max_nproc)
{
    const std::string expectedText =
        ToolContractJsonPrinterTests::MakeToolContractText("\"$max_nproc\"");

    // --------------------------------------------
    // setup interface & tool contract config
    // --------------------------------------------

    // clang-format off
    Interface i {
        "frobber",
        "Frobb your files in a most delightful, nobbly way",
        "3.14"
    };
    i.AddOptions({
        {"max_nlines", "n", "Max Number of lines to Copy", Option::IntType(10)},
        {"mode", "mode", "Run frobber in desired mode", Option::StringType("normal"), { "normal", "super", "secret" }},
        {"delta", "d", "Some delta for things", Option::IntType(1), { 1, 2, 3 }}
    });

    const std::string id{"frobber.tasks.dev_txt_app"};
    ToolContract::Task tcTask(id);
    tcTask.InputFileTypes({
        { "txt_in", "Txt file", "Generic Txt file", "PacBio.FileTypes.txt" }
    });
    tcTask.OutputFileTypes({
        { "txt_out", "Txt outfile", "Generic Output Txt file", "PacBio.FileTypes.txt", "output" }
    });
    tcTask.ResourceTypes({
        ResourceType::TMP_FILE,
        ResourceType::TMP_FILE,
        ResourceType::TMP_DIR
    });
    tcTask.Options({
        { "max_nlines","Max lines" },
        { "mode", "Frobbing mode" },
        { "delta", "Frobbing delta" }
    });
    tcTask.NumProcessors(Task::MAX_NPROC);
    // clang-format on

    ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);

    // --------------------------------------------
    // test tool contract output
    // --------------------------------------------

    std::ostringstream s;
    ToolContract::JsonPrinter::Print(i, s, 4);
    EXPECT_EQ(expectedText, s.str());
}
