#include <stdexcept>
#include <string>
#include <vector>
#include <map>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "biolcccexception.h"
#include "chemicalgroup.h"
#include "chemicalbasis.h"
#include "gradientpoint.h"
#include "gradient.h"
#include "chromoconditions.h"
#include "parsing.h"
#include "chain_model.h"
#include "rod_model.h"
#include "biolccc.h"

namespace py = pybind11;

#ifndef VERSION
#define VERSION "dev"
#endif

namespace {

static std::vector<std::string> chemical_group_keys() {
    return {
        "name", "label", "bindEnergy", "bindArea", "averageMass", "monoisotopicMass"
    };
}

static std::vector<std::string> chemical_basis_keys() {
    return {
        "chemicalGroups",
        "firstSolventDensity",
        "firstSolventAverageMass",
        "secondSolventDensity",
        "secondSolventAverageMass",
        "secondSolventBindEnergy",
        "adsorptionLayerWidth",
        "adsorptionLayerFactors",
        "kuhnLength",
        "monomerLength",
        "polymerModel",
        "snyderApproximation",
        "specialRodModel",
    };
}

static std::vector<std::string> gradient_point_keys() {
    return {"time", "concentrationB"};
}

static std::vector<std::string> chromo_conditions_keys() {
    return {
        "columnLength",
        "columnDiameter",
        "columnPoreSize",
        "gradient",
        "secondSolventConcentrationA",
        "secondSolventConcentrationB",
        "delayTime",
        "flowRate",
        "dV",
        "columnRelativeStrength",
        "columnVpToVtot",
        "columnPorosity",
        "temperature",
        "mixingCorrection",
    };
}

static bool equal_chemical_group(const BioLCCC::ChemicalGroup &a, const BioLCCC::ChemicalGroup &b) {
    return a.name() == b.name()
        && a.label() == b.label()
        && a.bindEnergy() == b.bindEnergy()
        && a.bindArea() == b.bindArea()
        && a.averageMass() == b.averageMass()
        && a.monoisotopicMass() == b.monoisotopicMass();
}

static bool equal_gradient_point(const BioLCCC::GradientPoint &a, const BioLCCC::GradientPoint &b) {
    return a.time() == b.time() && a.concentrationB() == b.concentrationB();
}

static bool equal_gradient(const BioLCCC::Gradient &a, const BioLCCC::Gradient &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (!equal_gradient_point(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

static bool equal_chemical_groups_map(
    const std::map<std::string, BioLCCC::ChemicalGroup> &a,
    const std::map<std::string, BioLCCC::ChemicalGroup> &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (const auto &item : a) {
        auto it = b.find(item.first);
        if (it == b.end()) {
            return false;
        }
        if (!equal_chemical_group(item.second, it->second)) {
            return false;
        }
    }
    return true;
}

static bool equal_chemical_basis(const BioLCCC::ChemicalBasis &a, const BioLCCC::ChemicalBasis &b) {
    return equal_chemical_groups_map(a.chemicalGroups(), b.chemicalGroups())
        && a.firstSolventDensity() == b.firstSolventDensity()
        && a.firstSolventAverageMass() == b.firstSolventAverageMass()
        && a.secondSolventDensity() == b.secondSolventDensity()
        && a.secondSolventAverageMass() == b.secondSolventAverageMass()
        && a.secondSolventBindEnergy() == b.secondSolventBindEnergy()
        && a.adsorptionLayerWidth() == b.adsorptionLayerWidth()
        && a.adsorptionLayerFactors() == b.adsorptionLayerFactors()
        && a.kuhnLength() == b.kuhnLength()
        && a.monomerLength() == b.monomerLength()
        && a.polymerModel() == b.polymerModel()
        && a.snyderApproximation() == b.snyderApproximation()
        && a.specialRodModel() == b.specialRodModel();
}

static bool equal_chromo_conditions(const BioLCCC::ChromoConditions &a, const BioLCCC::ChromoConditions &b) {
    return a.columnLength() == b.columnLength()
        && a.columnDiameter() == b.columnDiameter()
        && a.columnPoreSize() == b.columnPoreSize()
        && equal_gradient(a.gradient(), b.gradient())
        && a.secondSolventConcentrationA() == b.secondSolventConcentrationA()
        && a.secondSolventConcentrationB() == b.secondSolventConcentrationB()
        && a.delayTime() == b.delayTime()
        && a.flowRate() == b.flowRate()
        && a.dV() == b.dV()
        && a.columnRelativeStrength() == b.columnRelativeStrength()
        && a.columnVpToVtot() == b.columnVpToVtot()
        && a.columnPorosity() == b.columnPorosity()
        && a.temperature() == b.temperature()
        && a.mixingCorrection() == b.mixingCorrection();
}

static py::dict chemical_group_to_dict(const BioLCCC::ChemicalGroup &group) {
    py::dict d;
    d["name"] = group.name();
    d["label"] = group.label();
    d["bindEnergy"] = group.bindEnergy();
    d["bindArea"] = group.bindArea();
    d["averageMass"] = group.averageMass();
    d["monoisotopicMass"] = group.monoisotopicMass();
    return d;
}

static py::dict chemical_basis_getstate(const BioLCCC::ChemicalBasis &basis) {
    py::dict state;
    state["firstSolventDensity"] = basis.firstSolventDensity();
    state["firstSolventAverageMass"] = basis.firstSolventAverageMass();
    state["secondSolventDensity"] = basis.secondSolventDensity();
    state["secondSolventAverageMass"] = basis.secondSolventAverageMass();
    state["secondSolventBindEnergy"] = basis.secondSolventBindEnergy();
    state["adsorptionLayerWidth"] = basis.adsorptionLayerWidth();
    state["adsorptionLayerFactors"] = basis.adsorptionLayerFactors();
    state["kuhnLength"] = basis.kuhnLength();
    state["monomerLength"] = basis.monomerLength();
    state["polymerModel"] = basis.polymerModel();
    state["snyderApproximation"] = basis.snyderApproximation();
    state["specialRodModel"] = basis.specialRodModel();

    py::dict groups;
    for (const auto &item : basis.chemicalGroups()) {
        groups[item.first.c_str()] = chemical_group_to_dict(item.second);
    }
    state["chemicalGroups"] = groups;
    return state;
}

static void chemical_basis_setstate(BioLCCC::ChemicalBasis &basis, const py::dict &state) {
    basis.clearChemicalGroups();

    py::dict groups = state["chemicalGroups"].cast<py::dict>();
    for (auto item : groups) {
        py::object value = item.second.cast<py::object>();
        if (py::isinstance<py::dict>(value)) {
            py::dict g = value.cast<py::dict>();
            basis.addChemicalGroup(BioLCCC::ChemicalGroup(
                g["name"].cast<std::string>(),
                g["label"].cast<std::string>(),
                g["bindEnergy"].cast<double>(),
                g["averageMass"].cast<double>(),
                g["monoisotopicMass"].cast<double>(),
                g["bindArea"].cast<double>()
            ));
        } else {
            basis.addChemicalGroup(value.cast<BioLCCC::ChemicalGroup>());
        }
    }

    basis.setFirstSolventDensity(state["firstSolventDensity"].cast<double>());
    basis.setFirstSolventAverageMass(state["firstSolventAverageMass"].cast<double>());
    basis.setSecondSolventDensity(state["secondSolventDensity"].cast<double>());
    basis.setSecondSolventAverageMass(state["secondSolventAverageMass"].cast<double>());
    basis.setSecondSolventBindEnergy(state["secondSolventBindEnergy"].cast<double>());
    basis.setAdsorptionLayerWidth(state["adsorptionLayerWidth"].cast<double>());
    basis.setAdsorptionLayerFactors(state["adsorptionLayerFactors"].cast<std::vector<double>>());
    basis.setKuhnLength(state["kuhnLength"].cast<double>());
    basis.setMonomerLength(state["monomerLength"].cast<double>());
    basis.setPolymerModel(state["polymerModel"].cast<BioLCCC::PolymerModel>());
    basis.setSnyderApproximation(state["snyderApproximation"].cast<bool>());
    basis.setSpecialRodModel(state["specialRodModel"].cast<bool>());
}

static py::dict gradient_point_to_dict(const BioLCCC::GradientPoint &point) {
    py::dict d;
    d["time"] = point.time();
    d["concentrationB"] = point.concentrationB();
    return d;
}

static py::dict chromo_conditions_getstate(const BioLCCC::ChromoConditions &conditions) {
    py::dict state;
    state["columnLength"] = conditions.columnLength();
    state["columnDiameter"] = conditions.columnDiameter();
    state["columnPoreSize"] = conditions.columnPoreSize();
    state["secondSolventConcentrationA"] = conditions.secondSolventConcentrationA();
    state["secondSolventConcentrationB"] = conditions.secondSolventConcentrationB();
    state["delayTime"] = conditions.delayTime();
    state["flowRate"] = conditions.flowRate();
    state["dV"] = conditions.dV();
    state["columnRelativeStrength"] = conditions.columnRelativeStrength();
    state["columnVpToVtot"] = conditions.columnVpToVtot();
    state["columnPorosity"] = conditions.columnPorosity();
    state["temperature"] = conditions.temperature();
    state["mixingCorrection"] = conditions.mixingCorrection();

    py::list gradient;
    for (const auto &point : conditions.gradient()) {
        gradient.append(gradient_point_to_dict(point));
    }
    state["gradient"] = gradient;
    return state;
}

static void chromo_conditions_setstate(BioLCCC::ChromoConditions &conditions, const py::dict &state) {
    conditions.setColumnLength(state["columnLength"].cast<double>());
    conditions.setColumnDiameter(state["columnDiameter"].cast<double>());
    conditions.setColumnPoreSize(state["columnPoreSize"].cast<double>());
    conditions.setSecondSolventConcentrationA(state["secondSolventConcentrationA"].cast<double>());
    conditions.setSecondSolventConcentrationB(state["secondSolventConcentrationB"].cast<double>());
    conditions.setDelayTime(state["delayTime"].cast<double>());
    conditions.setFlowRate(state["flowRate"].cast<double>());
    conditions.setDV(state["dV"].cast<double>());
    conditions.setColumnRelativeStrength(state["columnRelativeStrength"].cast<double>());
    conditions.setColumnVpToVtot(state["columnVpToVtot"].cast<double>());
    conditions.setColumnPorosity(state["columnPorosity"].cast<double>());
    conditions.setTemperature(state["temperature"].cast<double>());
    conditions.setMixingCorrection(state["mixingCorrection"].cast<bool>());

    BioLCCC::Gradient gradient;
    for (auto item : state["gradient"].cast<py::list>()) {
        py::dict point = item.cast<py::dict>();
        gradient.addPoint(
            point["time"].cast<double>(),
            point["concentrationB"].cast<double>()
        );
    }
    conditions.setGradient(gradient);
}

} // namespace

PYBIND11_MODULE(_biolccc, m) {
    m.doc() = "pybind11 bindings for libBioLCCC";
    m.attr("VERSION") = VERSION;

    py::register_exception<BioLCCC::BioLCCCException>(m, "BioLCCCException");
    py::register_exception<BioLCCC::ChemicalBasisException>(m, "ChemicalBasisException");
    py::register_exception<BioLCCC::ChromoConditionsException>(m, "ChromoConditionsException");
    py::register_exception<BioLCCC::GradientException>(m, "GradientException");
    py::register_exception<BioLCCC::GradientPointException>(m, "GradientPointException");
    py::register_exception<BioLCCC::ParsingException>(m, "ParsingException");

    py::enum_<BioLCCC::PolymerModel>(m, "PolymerModel")
        .value("ROD", BioLCCC::ROD)
        .value("CHAIN", BioLCCC::CHAIN)
        .export_values();

    py::enum_<BioLCCC::PredefinedChemicalBasis>(m, "PredefinedChemicalBasis")
        .value("RP_ACN_TFA_CHAIN", BioLCCC::RP_ACN_TFA_CHAIN)
        .value("RP_ACN_FA_ROD", BioLCCC::RP_ACN_FA_ROD)
        .export_values();

    py::class_<BioLCCC::ChemicalGroup>(m, "ChemicalGroup")
        .def(py::init<std::string, std::string, double, double, double, double>(),
            py::arg("name") = "",
            py::arg("label") = "",
            py::arg("bindEnergy") = 0.0,
            py::arg("averageMass") = 0.0,
            py::arg("monoisotopicMass") = 0.0,
            py::arg("bindArea") = 1.0)
        .def("name", &BioLCCC::ChemicalGroup::name)
        .def("label", &BioLCCC::ChemicalGroup::label)
        .def("averageMass", &BioLCCC::ChemicalGroup::averageMass)
        .def("monoisotopicMass", &BioLCCC::ChemicalGroup::monoisotopicMass)
        .def("bindEnergy", &BioLCCC::ChemicalGroup::bindEnergy)
        .def("bindArea", &BioLCCC::ChemicalGroup::bindArea)
        .def("isNTerminal", &BioLCCC::ChemicalGroup::isNTerminal)
        .def("isCTerminal", &BioLCCC::ChemicalGroup::isCTerminal)
        .def("isAminoAcid", &BioLCCC::ChemicalGroup::isAminoAcid)
        .def("setName", &BioLCCC::ChemicalGroup::setName)
        .def("setBindEnergy", &BioLCCC::ChemicalGroup::setBindEnergy)
        .def("setBindArea", &BioLCCC::ChemicalGroup::setBindArea)
        .def("setAverageMass", &BioLCCC::ChemicalGroup::setAverageMass)
        .def("setMonoisotopicMass", &BioLCCC::ChemicalGroup::setMonoisotopicMass)
        .def("keys", [](const BioLCCC::ChemicalGroup &) { return chemical_group_keys(); })
        .def("__len__", [](const BioLCCC::ChemicalGroup &) { return static_cast<py::ssize_t>(chemical_group_keys().size()); })
        .def("__iter__", [](const BioLCCC::ChemicalGroup &) {
            auto keys = py::cast(chemical_group_keys());
            return py::iter(keys);
        })
        .def("__contains__", [](const BioLCCC::ChemicalGroup &, const std::string &key) {
            for (const auto &k : chemical_group_keys()) {
                if (k == key) {
                    return true;
                }
            }
            return false;
        })
        .def("__getitem__", [](const BioLCCC::ChemicalGroup &self, const std::string &key) -> py::object {
            if (key == "name") return py::cast(self.name());
            if (key == "label") return py::cast(self.label());
            if (key == "bindEnergy") return py::cast(self.bindEnergy());
            if (key == "bindArea") return py::cast(self.bindArea());
            if (key == "averageMass") return py::cast(self.averageMass());
            if (key == "monoisotopicMass") return py::cast(self.monoisotopicMass());
            throw py::key_error("Unknown ChemicalGroup key: " + key);
        })
        .def("__setitem__", [](BioLCCC::ChemicalGroup &self, const std::string &key, py::object value) {
            if (key == "name") {
                self.setName(value.cast<std::string>());
                return;
            }
            if (key == "label") {
                throw std::runtime_error("Label cannot be set");
            }
            if (key == "bindEnergy") {
                self.setBindEnergy(value.cast<double>());
                return;
            }
            if (key == "bindArea") {
                self.setBindArea(value.cast<double>());
                return;
            }
            if (key == "averageMass") {
                self.setAverageMass(value.cast<double>());
                return;
            }
            if (key == "monoisotopicMass") {
                self.setMonoisotopicMass(value.cast<double>());
                return;
            }
            throw py::key_error("Unknown ChemicalGroup key: " + key);
        })
        .def("__delitem__", [](BioLCCC::ChemicalGroup &, const std::string &) { })
        .def("__eq__", [](const BioLCCC::ChemicalGroup &self, const BioLCCC::ChemicalGroup &other) {
            return equal_chemical_group(self, other);
        })
        .def("__repr__", [](const BioLCCC::ChemicalGroup &self) {
            return py::str(chemical_group_to_dict(self)).cast<std::string>();
        })
        .def("__str__", [](const BioLCCC::ChemicalGroup &self) {
            return py::str(chemical_group_to_dict(self)).cast<std::string>();
        });

    py::class_<BioLCCC::GradientPoint>(m, "GradientPoint")
        .def(py::init<double, double>(), py::arg("time") = 0.0, py::arg("concentrationB") = 0.0)
        .def("time", &BioLCCC::GradientPoint::time)
        .def("concentrationB", &BioLCCC::GradientPoint::concentrationB)
        .def("setTime", &BioLCCC::GradientPoint::setTime)
        .def("setConcentrationB", &BioLCCC::GradientPoint::setConcentrationB)
        .def("keys", [](const BioLCCC::GradientPoint &) { return gradient_point_keys(); })
        .def("__len__", [](const BioLCCC::GradientPoint &) { return static_cast<py::ssize_t>(gradient_point_keys().size()); })
        .def("__iter__", [](const BioLCCC::GradientPoint &) {
            auto keys = py::cast(gradient_point_keys());
            return py::iter(keys);
        })
        .def("__contains__", [](const BioLCCC::GradientPoint &, const std::string &key) {
            for (const auto &k : gradient_point_keys()) {
                if (k == key) {
                    return true;
                }
            }
            return false;
        })
        .def("__getitem__", [](const BioLCCC::GradientPoint &self, const std::string &key) -> py::object {
            if (key == "time") return py::cast(self.time());
            if (key == "concentrationB") return py::cast(self.concentrationB());
            throw py::key_error("Unknown GradientPoint key: " + key);
        })
        .def("__setitem__", [](BioLCCC::GradientPoint &self, const std::string &key, py::object value) {
            if (key == "time") {
                self.setTime(value.cast<double>());
                return;
            }
            if (key == "concentrationB") {
                self.setConcentrationB(value.cast<double>());
                return;
            }
            throw py::key_error("Unknown GradientPoint key: " + key);
        })
        .def("__delitem__", [](BioLCCC::GradientPoint &, const std::string &) { })
        .def("__eq__", [](const BioLCCC::GradientPoint &self, const BioLCCC::GradientPoint &other) {
            return equal_gradient_point(self, other);
        })
        .def("__repr__", [](const BioLCCC::GradientPoint &self) {
            return py::str(gradient_point_to_dict(self)).cast<std::string>();
        })
        .def("__str__", [](const BioLCCC::GradientPoint &self) {
            return py::str(gradient_point_to_dict(self)).cast<std::string>();
        });

    py::class_<BioLCCC::Gradient>(m, "Gradient")
        .def(py::init<>())
        .def(py::init<double, double, double>(),
            py::arg("initialConcentrationB"),
            py::arg("finalConcentrationB"),
            py::arg("time"))
        .def("addPoint", py::overload_cast<BioLCCC::GradientPoint>(&BioLCCC::Gradient::addPoint), py::arg("point"))
        .def("addPoint", py::overload_cast<double, double>(&BioLCCC::Gradient::addPoint), py::arg("time"), py::arg("concentrationB"))
        .def("__len__", [](const BioLCCC::Gradient &self) { return static_cast<py::ssize_t>(self.size()); })
        .def("__iter__", [](BioLCCC::Gradient &self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("__getitem__", [](const BioLCCC::Gradient &self, py::ssize_t idx) -> const BioLCCC::GradientPoint & {
            auto size = static_cast<py::ssize_t>(self.size());
            if (idx < 0) {
                idx += size;
            }
            if (idx < 0 || idx >= size) {
                throw py::index_error("Gradient index out of range");
            }
            return self[static_cast<size_t>(idx)];
        }, py::return_value_policy::copy)
        .def("__eq__", [](const BioLCCC::Gradient &self, const BioLCCC::Gradient &other) {
            return equal_gradient(self, other);
        })
        .def("__repr__", [](const BioLCCC::Gradient &self) {
            py::list points;
            for (const auto &p : self) {
                points.append(gradient_point_to_dict(p));
            }
            return py::str(points).cast<std::string>();
        })
        .def("__str__", [](const BioLCCC::Gradient &self) {
            py::list points;
            for (const auto &p : self) {
                points.append(gradient_point_to_dict(p));
            }
            return py::str(points).cast<std::string>();
        });

    py::class_<std::map<std::string, BioLCCC::ChemicalGroup>>(m, "ChemicalGroupsMap")
        .def(py::init<>())
        .def("__len__", [](const std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            return static_cast<py::ssize_t>(self.size());
        })
        .def("__iter__", [](std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            return py::make_key_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())
        .def("keys", [](std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            py::list keys;
            for (const auto &item : self) {
                keys.append(item.first);
            }
            return keys;
        })
        .def("items", [](std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            py::list items;
            for (auto &item : self) {
                items.append(py::make_tuple(item.first, item.second));
            }
            return items;
        })
        .def("values", [](std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            py::list values;
            for (auto &item : self) {
                values.append(item.second);
            }
            return values;
        })
        .def("__contains__", [](const std::map<std::string, BioLCCC::ChemicalGroup> &self, const std::string &key) {
            return self.find(key) != self.end();
        })
        .def("__getitem__", [](std::map<std::string, BioLCCC::ChemicalGroup> &self, const std::string &key) -> BioLCCC::ChemicalGroup& {
            auto it = self.find(key);
            if (it == self.end()) {
                throw py::key_error("Unknown chemical group label: " + key);
            }
            return it->second;
        }, py::return_value_policy::reference_internal)
        .def("__setitem__", [](std::map<std::string, BioLCCC::ChemicalGroup> &self, const std::string &key, const BioLCCC::ChemicalGroup &value) {
            self[key] = value;
        })
        .def("__delitem__", [](std::map<std::string, BioLCCC::ChemicalGroup> &self, const std::string &key) {
            auto it = self.find(key);
            if (it == self.end()) {
                throw py::key_error("Unknown chemical group label: " + key);
            }
            self.erase(it);
        })
        .def("__repr__", [](const std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            py::dict d;
            for (const auto &item : self) {
                d[item.first.c_str()] = item.second;
            }
            return py::str(d).cast<std::string>();
        })
        .def("__str__", [](const std::map<std::string, BioLCCC::ChemicalGroup> &self) {
            py::dict d;
            for (const auto &item : self) {
                d[item.first.c_str()] = item.second;
            }
            return py::str(d).cast<std::string>();
        });

    py::class_<BioLCCC::ChemicalBasis>(m, "ChemicalBasis")
        .def(py::init<>())
        .def(py::init<BioLCCC::PredefinedChemicalBasis>(), py::arg("predefinedChemicalBasisId"))
        .def("chemicalGroups",
            py::overload_cast<>(&BioLCCC::ChemicalBasis::chemicalGroups),
            py::return_value_policy::reference_internal)
        .def("defaultNTerminus",
            &BioLCCC::ChemicalBasis::defaultNTerminus,
            py::return_value_policy::copy)
        .def("defaultCTerminus",
            &BioLCCC::ChemicalBasis::defaultCTerminus,
            py::return_value_policy::copy)
        .def("addChemicalGroup", &BioLCCC::ChemicalBasis::addChemicalGroup)
        .def("removeChemicalGroup", &BioLCCC::ChemicalBasis::removeChemicalGroup)
        .def("clearChemicalGroups", &BioLCCC::ChemicalBasis::clearChemicalGroups)
        .def("secondSolventBindEnergy", &BioLCCC::ChemicalBasis::secondSolventBindEnergy)
        .def("setSecondSolventBindEnergy", &BioLCCC::ChemicalBasis::setSecondSolventBindEnergy)
        .def("setPolymerModel", &BioLCCC::ChemicalBasis::setPolymerModel)
        .def("setModel", &BioLCCC::ChemicalBasis::setPolymerModel)
        .def("polymerModel", &BioLCCC::ChemicalBasis::polymerModel)
        .def("monomerLength", &BioLCCC::ChemicalBasis::monomerLength)
        .def("setMonomerLength", &BioLCCC::ChemicalBasis::setMonomerLength)
        .def("kuhnLength", &BioLCCC::ChemicalBasis::kuhnLength)
        .def("setKuhnLength", &BioLCCC::ChemicalBasis::setKuhnLength)
        .def("adsorptionLayerWidth", &BioLCCC::ChemicalBasis::adsorptionLayerWidth)
        .def("setAdsorptionLayerWidth", &BioLCCC::ChemicalBasis::setAdsorptionLayerWidth)
        .def("adsorptionLayerFactors", &BioLCCC::ChemicalBasis::adsorptionLayerFactors)
        .def("setAdsorptionLayerFactors", &BioLCCC::ChemicalBasis::setAdsorptionLayerFactors)
        .def("snyderApproximation", &BioLCCC::ChemicalBasis::snyderApproximation)
        .def("setSnyderApproximation", &BioLCCC::ChemicalBasis::setSnyderApproximation)
        .def("specialRodModel", &BioLCCC::ChemicalBasis::specialRodModel)
        .def("setSpecialRodModel", &BioLCCC::ChemicalBasis::setSpecialRodModel)
        .def("neglectPartiallyDesorbedStates", &BioLCCC::ChemicalBasis::neglectPartiallyDesorbedStates)
        .def("setNeglectPartiallyDesorbedStates", &BioLCCC::ChemicalBasis::setNeglectPartiallyDesorbedStates)
        .def("firstSolventDensity", &BioLCCC::ChemicalBasis::firstSolventDensity)
        .def("setFirstSolventDensity", &BioLCCC::ChemicalBasis::setFirstSolventDensity)
        .def("secondSolventDensity", &BioLCCC::ChemicalBasis::secondSolventDensity)
        .def("setSecondSolventDensity", &BioLCCC::ChemicalBasis::setSecondSolventDensity)
        .def("firstSolventAverageMass", &BioLCCC::ChemicalBasis::firstSolventAverageMass)
        .def("setFirstSolventAverageMass", &BioLCCC::ChemicalBasis::setFirstSolventAverageMass)
        .def("secondSolventAverageMass", &BioLCCC::ChemicalBasis::secondSolventAverageMass)
        .def("setSecondSolventAverageMass", &BioLCCC::ChemicalBasis::setSecondSolventAverageMass)
        .def("setPredefinedChemicalBasis", &BioLCCC::ChemicalBasis::setPredefinedChemicalBasis)
        .def("keys", [](const BioLCCC::ChemicalBasis &) { return chemical_basis_keys(); })
        .def("__len__", [](const BioLCCC::ChemicalBasis &) { return static_cast<py::ssize_t>(chemical_basis_keys().size()); })
        .def("__iter__", [](const BioLCCC::ChemicalBasis &) {
            auto keys = py::cast(chemical_basis_keys());
            return py::iter(keys);
        })
        .def("__contains__", [](const BioLCCC::ChemicalBasis &, const std::string &key) {
            for (const auto &k : chemical_basis_keys()) {
                if (k == key) {
                    return true;
                }
            }
            return false;
        })
        .def("__getitem__", [](BioLCCC::ChemicalBasis &self, const std::string &key) -> py::object {
            if (key == "chemicalGroups") return py::cast(self.chemicalGroups(), py::return_value_policy::reference);
            if (key == "firstSolventDensity") return py::cast(self.firstSolventDensity());
            if (key == "firstSolventAverageMass") return py::cast(self.firstSolventAverageMass());
            if (key == "secondSolventDensity") return py::cast(self.secondSolventDensity());
            if (key == "secondSolventAverageMass") return py::cast(self.secondSolventAverageMass());
            if (key == "secondSolventBindEnergy") return py::cast(self.secondSolventBindEnergy());
            if (key == "adsorptionLayerWidth") return py::cast(self.adsorptionLayerWidth());
            if (key == "adsorptionLayerFactors") return py::cast(self.adsorptionLayerFactors());
            if (key == "kuhnLength") return py::cast(self.kuhnLength());
            if (key == "monomerLength") return py::cast(self.monomerLength());
            if (key == "polymerModel") return py::cast(self.polymerModel());
            if (key == "snyderApproximation") return py::cast(self.snyderApproximation());
            if (key == "specialRodModel") return py::cast(self.specialRodModel());
            throw py::key_error("Unknown ChemicalBasis key: " + key);
        })
        .def("__setitem__", [](BioLCCC::ChemicalBasis &self, const std::string &key, py::object value) {
            if (key == "chemicalGroups") {
                self.clearChemicalGroups();
                py::dict groups = value.cast<py::dict>();
                for (auto item : groups) {
                    py::object group_obj = item.second.cast<py::object>();
                    if (py::isinstance<py::dict>(group_obj)) {
                        py::dict group = group_obj.cast<py::dict>();
                        self.addChemicalGroup(BioLCCC::ChemicalGroup(
                            group["name"].cast<std::string>(),
                            group["label"].cast<std::string>(),
                            group["bindEnergy"].cast<double>(),
                            group["averageMass"].cast<double>(),
                            group["monoisotopicMass"].cast<double>(),
                            group["bindArea"].cast<double>()));
                    } else {
                        self.addChemicalGroup(group_obj.cast<BioLCCC::ChemicalGroup>());
                    }
                }
                return;
            }
            if (key == "firstSolventDensity") { self.setFirstSolventDensity(value.cast<double>()); return; }
            if (key == "firstSolventAverageMass") { self.setFirstSolventAverageMass(value.cast<double>()); return; }
            if (key == "secondSolventDensity") { self.setSecondSolventDensity(value.cast<double>()); return; }
            if (key == "secondSolventAverageMass") { self.setSecondSolventAverageMass(value.cast<double>()); return; }
            if (key == "secondSolventBindEnergy") { self.setSecondSolventBindEnergy(value.cast<double>()); return; }
            if (key == "adsorptionLayerWidth") { self.setAdsorptionLayerWidth(value.cast<double>()); return; }
            if (key == "adsorptionLayerFactors") { self.setAdsorptionLayerFactors(value.cast<std::vector<double>>()); return; }
            if (key == "kuhnLength") { self.setKuhnLength(value.cast<double>()); return; }
            if (key == "monomerLength") { self.setMonomerLength(value.cast<double>()); return; }
            if (key == "polymerModel") { self.setPolymerModel(value.cast<BioLCCC::PolymerModel>()); return; }
            if (key == "snyderApproximation") { self.setSnyderApproximation(value.cast<bool>()); return; }
            if (key == "specialRodModel") { self.setSpecialRodModel(value.cast<bool>()); return; }
            throw py::key_error("Unknown ChemicalBasis key: " + key);
        })
        .def("__delitem__", [](BioLCCC::ChemicalBasis &, const std::string &) { })
        .def("__eq__", [](const BioLCCC::ChemicalBasis &self, const BioLCCC::ChemicalBasis &other) {
            return equal_chemical_basis(self, other);
        })
        .def("__repr__", [](const BioLCCC::ChemicalBasis &self) {
            return py::str(chemical_basis_getstate(self)).cast<std::string>();
        })
        .def("__str__", [](const BioLCCC::ChemicalBasis &self) {
            return py::str(chemical_basis_getstate(self)).cast<std::string>();
        })
        .def("__getstate__", [](const BioLCCC::ChemicalBasis &self) {
            return chemical_basis_getstate(self);
        })
        .def("__setstate__", [](BioLCCC::ChemicalBasis &self, const py::dict &state) {
            chemical_basis_setstate(self, state);
        })
        .def("__reduce__", [](const BioLCCC::ChemicalBasis &self) {
            py::object cls = py::module_::import("pyteomics.biolccc").attr("ChemicalBasis");
            return py::make_tuple(cls, py::make_tuple(), chemical_basis_getstate(self));
        });

    py::class_<BioLCCC::ChromoConditions>(m, "ChromoConditions")
        .def(py::init<double, double, double, BioLCCC::Gradient, double, double, double, double, double, double, double, double, double>(),
            py::arg("columnLength") = 150.0,
            py::arg("columnDiameter") = 0.075,
            py::arg("columnPoreSize") = 100.0,
            py::arg("gradient") = BioLCCC::Gradient(0.0, 50.0, 60.0),
            py::arg("secondSolventConcentrationA") = 2.0,
            py::arg("secondSolventConcentrationB") = 80.0,
            py::arg("delayTime") = 0.0,
            py::arg("flowRate") = 0.0003,
            py::arg("dV") = 0.0,
            py::arg("columnRelativeStrength") = 1.0,
            py::arg("columnVpToVtot") = 0.5,
            py::arg("columnPorosity") = 0.9,
            py::arg("temperature") = 293.0)
        .def("columnLength", &BioLCCC::ChromoConditions::columnLength)
        .def("setColumnLength", &BioLCCC::ChromoConditions::setColumnLength)
        .def("columnDiameter", &BioLCCC::ChromoConditions::columnDiameter)
        .def("setColumnDiameter", &BioLCCC::ChromoConditions::setColumnDiameter)
        .def("columnPoreSize", &BioLCCC::ChromoConditions::columnPoreSize)
        .def("setColumnPoreSize", &BioLCCC::ChromoConditions::setColumnPoreSize)
        .def("columnVpToVtot", &BioLCCC::ChromoConditions::columnVpToVtot)
        .def("setColumnVpToVtot", &BioLCCC::ChromoConditions::setColumnVpToVtot)
        .def("columnPorosity", &BioLCCC::ChromoConditions::columnPorosity)
        .def("setColumnPorosity", &BioLCCC::ChromoConditions::setColumnPorosity)
        .def("columnTotalVolume", &BioLCCC::ChromoConditions::columnTotalVolume)
        .def("columnInterstitialVolume", &BioLCCC::ChromoConditions::columnInterstitialVolume)
        .def("columnPoreVolume", &BioLCCC::ChromoConditions::columnPoreVolume)
        .def("temperature", &BioLCCC::ChromoConditions::temperature)
        .def("setTemperature", &BioLCCC::ChromoConditions::setTemperature)
        .def("columnRelativeStrength", &BioLCCC::ChromoConditions::columnRelativeStrength)
        .def("setColumnRelativeStrength", &BioLCCC::ChromoConditions::setColumnRelativeStrength)
        .def("flowRate", &BioLCCC::ChromoConditions::flowRate)
        .def("setFlowRate", &BioLCCC::ChromoConditions::setFlowRate)
        .def("dV", &BioLCCC::ChromoConditions::dV)
        .def("setDV", &BioLCCC::ChromoConditions::setDV)
        .def("delayTime", &BioLCCC::ChromoConditions::delayTime)
        .def("setDelayTime", &BioLCCC::ChromoConditions::setDelayTime)
        .def("secondSolventConcentrationA", &BioLCCC::ChromoConditions::secondSolventConcentrationA)
        .def("setSecondSolventConcentrationA", &BioLCCC::ChromoConditions::setSecondSolventConcentrationA)
        .def("secondSolventConcentrationB", &BioLCCC::ChromoConditions::secondSolventConcentrationB)
        .def("setSecondSolventConcentrationB", &BioLCCC::ChromoConditions::setSecondSolventConcentrationB)
        .def("gradient", &BioLCCC::ChromoConditions::gradient)
        .def("setGradient", &BioLCCC::ChromoConditions::setGradient)
        .def("mixingCorrection", &BioLCCC::ChromoConditions::mixingCorrection)
        .def("setMixingCorrection", &BioLCCC::ChromoConditions::setMixingCorrection)
        .def("SSConcentrations", &BioLCCC::ChromoConditions::SSConcentrations)
        .def("keys", [](const BioLCCC::ChromoConditions &) { return chromo_conditions_keys(); })
        .def("__len__", [](const BioLCCC::ChromoConditions &) { return static_cast<py::ssize_t>(chromo_conditions_keys().size()); })
        .def("__iter__", [](const BioLCCC::ChromoConditions &) {
            auto keys = py::cast(chromo_conditions_keys());
            return py::iter(keys);
        })
        .def("__contains__", [](const BioLCCC::ChromoConditions &, const std::string &key) {
            for (const auto &k : chromo_conditions_keys()) {
                if (k == key) {
                    return true;
                }
            }
            return false;
        })
        .def("__getitem__", [](const BioLCCC::ChromoConditions &self, const std::string &key) -> py::object {
            if (key == "columnLength") return py::cast(self.columnLength());
            if (key == "columnDiameter") return py::cast(self.columnDiameter());
            if (key == "columnPoreSize") return py::cast(self.columnPoreSize());
            if (key == "columnRelativeStrength") return py::cast(self.columnRelativeStrength());
            if (key == "columnVpToVtot") return py::cast(self.columnVpToVtot());
            if (key == "columnPorosity") return py::cast(self.columnPorosity());
            if (key == "delayTime") return py::cast(self.delayTime());
            if (key == "dV") return py::cast(self.dV());
            if (key == "flowRate") return py::cast(self.flowRate());
            if (key == "gradient") return py::cast(self.gradient());
            if (key == "secondSolventConcentrationA") return py::cast(self.secondSolventConcentrationA());
            if (key == "secondSolventConcentrationB") return py::cast(self.secondSolventConcentrationB());
            if (key == "temperature") return py::cast(self.temperature());
            if (key == "mixingCorrection") return py::cast(self.mixingCorrection());
            throw py::key_error("Unknown ChromoConditions key: " + key);
        })
        .def("__setitem__", [](BioLCCC::ChromoConditions &self, const std::string &key, py::object value) {
            if (key == "columnLength") { self.setColumnLength(value.cast<double>()); return; }
            if (key == "columnDiameter") { self.setColumnDiameter(value.cast<double>()); return; }
            if (key == "columnPoreSize") { self.setColumnPoreSize(value.cast<double>()); return; }
            if (key == "columnRelativeStrength") { self.setColumnRelativeStrength(value.cast<double>()); return; }
            if (key == "columnVpToVtot") { self.setColumnVpToVtot(value.cast<double>()); return; }
            if (key == "columnPorosity") { self.setColumnPorosity(value.cast<double>()); return; }
            if (key == "delayTime") { self.setDelayTime(value.cast<double>()); return; }
            if (key == "dV") { self.setDV(value.cast<double>()); return; }
            if (key == "flowRate") { self.setFlowRate(value.cast<double>()); return; }
            if (key == "gradient") { self.setGradient(value.cast<BioLCCC::Gradient>()); return; }
            if (key == "secondSolventConcentrationA") { self.setSecondSolventConcentrationA(value.cast<double>()); return; }
            if (key == "secondSolventConcentrationB") { self.setSecondSolventConcentrationB(value.cast<double>()); return; }
            if (key == "temperature") { self.setTemperature(value.cast<double>()); return; }
            if (key == "mixingCorrection") { self.setMixingCorrection(value.cast<bool>()); return; }
            throw py::key_error("Unknown ChromoConditions key: " + key);
        })
        .def("__delitem__", [](BioLCCC::ChromoConditions &, const std::string &) { })
        .def("__eq__", [](const BioLCCC::ChromoConditions &self, const BioLCCC::ChromoConditions &other) {
            return equal_chromo_conditions(self, other);
        })
        .def("__repr__", [](const BioLCCC::ChromoConditions &self) {
            return py::str(chromo_conditions_getstate(self)).cast<std::string>();
        })
        .def("__str__", [](const BioLCCC::ChromoConditions &self) {
            return py::str(chromo_conditions_getstate(self)).cast<std::string>();
        })
        .def("__getstate__", [](const BioLCCC::ChromoConditions &self) {
            return chromo_conditions_getstate(self);
        })
        .def("__setstate__", [](BioLCCC::ChromoConditions &self, const py::dict &state) {
            chromo_conditions_setstate(self, state);
        })
        .def("__reduce__", [](const BioLCCC::ChromoConditions &self) {
            py::object cls = py::module_::import("pyteomics.biolccc").attr("ChromoConditions");
            return py::make_tuple(cls, py::make_tuple(), chromo_conditions_getstate(self));
        });

    m.def("calculateRT", &BioLCCC::calculateRT,
        py::arg("sequence"),
        py::arg("chemBasis"),
        py::arg("conditions") = BioLCCC::standardChromoConditions,
        py::arg("numInterpolationPoints") = 0,
        py::arg("continueGradient") = true,
        py::arg("backwardCompatibility") = false);

    m.def("calculateKd", &BioLCCC::calculateKd,
        py::arg("sequence"),
        py::arg("secondSolventConcentration"),
        py::arg("chemBasis"),
        py::arg("columnPoreSize") = 100.0,
        py::arg("columnRelativeStrength") = 1.0,
        py::arg("temperature") = 293.0);

    m.def("calculateAverageMass", &BioLCCC::calculateAverageMass,
        py::arg("sequence"), py::arg("chemBasis"));

    m.def("calculateMonoisotopicMass", &BioLCCC::calculateMonoisotopicMass,
        py::arg("sequence"), py::arg("chemBasis"));

    m.def("parseSequence",
        [](const std::string &sequence, const BioLCCC::ChemicalBasis &chemBasis) {
            return BioLCCC::parseSequence(sequence, chemBasis);
        },
        py::arg("sequence"), py::arg("chemBasis") = BioLCCC::rpAcnTfaChain);

    m.def("calculateMonomerEnergyProfile", &BioLCCC::calculateMonomerEnergyProfile,
        py::arg("parsedSequence"),
        py::arg("chemBasis"),
        py::arg("secondSolventConcentration"),
        py::arg("columnRelativeStrength"),
        py::arg("temperature"));

    m.def("calculateSegmentEnergyProfile", &BioLCCC::calculateSegmentEnergyProfile,
        py::arg("monomerEnergyProfile"),
        py::arg("monomerLength"),
        py::arg("kuhnLength"));

    m.attr("standardChromoConditions") = BioLCCC::standardChromoConditions;
    m.attr("rpAcnTfaChain") = BioLCCC::rpAcnTfaChain;
    m.attr("rpAcnFaRod") = BioLCCC::rpAcnFaRod;
}
