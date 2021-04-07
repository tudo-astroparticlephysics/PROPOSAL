#include <filesystem>
#include <fstream>

auto getTestFiles(std::string filename)
{
    auto directory = std::filesystem::path("tests/TestFiles");

    if (auto test_p = std::getenv("PROPOSAL_TEST_FILES"))
        directory = std::filesystem::path(test_p);

    if (std::filesystem::exists(directory / filename))
        return std::ifstream { directory / filename };

    std::ostringstream error;
    error << "Unable to find TestFile under path: " << directory/filename <<
    ". Try setting the PROPOSAL_TEST_FILES environment variable.";
    throw std::invalid_argument( error.str() );
}
