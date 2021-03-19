#include <filesystem>
#include <fstream>

auto getTestFiles(std::string filename)
{
    if (auto test_p = std::getenv("PROPOSAL_TEST_FILES")) {
        auto test_path = std::filesystem::path(test_p);
        return std::ifstream { test_path / filename };
    }

    auto default_dir = std::filesystem::path("tests/TestFiles");
    if (std::filesystem::exists(default_dir / filename))
        return std::ifstream { default_dir / filename };

    auto error_msg = std::string("No testfiles found. Please set the "
                                 "PROPOSAL_TEST_FILES environment variable.");
    throw std::invalid_argument(error_msg.c_str());
}
