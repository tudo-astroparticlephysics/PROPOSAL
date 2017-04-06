#include <boost/python.hpp>
#include "python/test.h"

BOOST_PYTHON_MODULE(pyPROPOSAL_ext)
{
    using namespace boost::python;
    def( "greet", greet );
}
