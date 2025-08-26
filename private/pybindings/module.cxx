// Currently using trivial pybindings.

#include <icetray/load_project.h>
#include <boost/python.hpp>
#include <pythia-wrapper/I3PythiaWrapperModule.hh>
#include "pybindings.hpp"

namespace bp = boost::python;

void register_pythia_wrapper() {

}

I3_PYTHON_MODULE(pythia_wrapper) {
    load_project("pythia_wrapper",false);

	register_pythia_wrapper();
}
