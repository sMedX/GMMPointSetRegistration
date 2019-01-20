#include <pybind11/pybind11.h>

class Pet {
	public:
		Pet(const std::string &name_){ name = name_;}
		void setName(const std::string &name_) { name = name_; }
		const std::string &getName() const { return name; }

	private:
		std::string name;
};

namespace py = pybind11;

PYBIND11_MODULE(example, m) {
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &>())
        .def("setName", &Pet::setName)
        .def("getName", &Pet::getName);
}
