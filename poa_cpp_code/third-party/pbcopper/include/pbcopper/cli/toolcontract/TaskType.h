// Author: Derek Barnett

#ifndef PBCOPPER_CLI_TOOLCONTRACT_TASKTYPE_H
#define PBCOPPER_CLI_TOOLCONTRACT_TASKTYPE_H

#include <pbcopper/PbcopperConfig.h>

namespace PacBio {
namespace CLI {
namespace ToolContract {

enum class TaskType
{
    STANDARD,
    SCATTERED,
    GATHERED
};

}  // namespace ToolContract
}  // namespace CLI
}  // namespace PacBio

#endif  // PBCOPPER_CLI_TOOLCONTRACT_TASKTYPE_H
