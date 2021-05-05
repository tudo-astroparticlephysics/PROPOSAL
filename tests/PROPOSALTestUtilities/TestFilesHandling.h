#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>

auto getTestFiles(std::string filename)
{
    auto directory = boost::filesystem::path("tests/TestFiles");

    if (auto test_p = std::getenv("PROPOSAL_TEST_FILES"))
        directory = boost::filesystem::path(test_p);

    if (boost::filesystem::exists(directory / filename))
        return std::ifstream { directory / filename };

    std::ostringstream error;
    error << "Unable to find TestFile under path: " << directory / filename
          << ". Try setting the PROPOSAL_TEST_FILES environment variable.";
    throw std::invalid_argument(error.str());
}
