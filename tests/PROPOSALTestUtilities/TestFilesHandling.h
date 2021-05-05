#include <boost/filesystem.hpp>
#include <fstream>
#include <sstream>

inline auto getTestFile(std::string filename, std::ifstream& file)
{
    auto directory = boost::filesystem::path("tests/TestFiles");

    if (auto test_p = std::getenv("PROPOSAL_TEST_FILES"))
        directory = boost::filesystem::path(test_p);

    if (!boost::filesystem::exists(directory / filename)) {
        std::ostringstream error;
        error << "Unable to find TestFile under path: " << directory / filename
              << ". Try setting the PROPOSAL_TEST_FILES environment variable.";
        throw std::invalid_argument(error.str());
    }
    auto file_path = directory / filename;
    file.open(file_path.generic_string());
}
