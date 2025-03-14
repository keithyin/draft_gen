#include <gtest/gtest.h>
#include <pbcopper/json/JSON.h>
using namespace PacBio;
using namespace PacBio::JSON;

TEST(JSON_Json, default_constructed_object_is_null)
{
    Json j;
    EXPECT_TRUE(j.is_null());
}

// nlohmann tests are __massive__ . Probably no need to check all here, but come
// back and test the main features we use
