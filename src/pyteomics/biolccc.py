from collections.abc import MutableMapping

from pyteomics import _biolccc

VERSION = _biolccc.VERSION

BioLCCCException = _biolccc.BioLCCCException
ChemicalBasisException = _biolccc.ChemicalBasisException
ChromoConditionsException = _biolccc.ChromoConditionsException
GradientException = _biolccc.GradientException
GradientPointException = _biolccc.GradientPointException
ParsingException = _biolccc.ParsingException

ROD = _biolccc.ROD
CHAIN = _biolccc.CHAIN
RP_ACN_TFA_CHAIN = _biolccc.RP_ACN_TFA_CHAIN
RP_ACN_FA_ROD = _biolccc.RP_ACN_FA_ROD

ChemicalGroup = _biolccc.ChemicalGroup
ChemicalBasis = _biolccc.ChemicalBasis
GradientPoint = _biolccc.GradientPoint
Gradient = _biolccc.Gradient
ChromoConditions = _biolccc.ChromoConditions

calculateRT = _biolccc.calculateRT
calculateKd = _biolccc.calculateKd
calculateAverageMass = _biolccc.calculateAverageMass
calculateMonoisotopicMass = _biolccc.calculateMonoisotopicMass
parseSequence = _biolccc.parseSequence
calculateMonomerEnergyProfile = _biolccc.calculateMonomerEnergyProfile
calculateSegmentEnergyProfile = _biolccc.calculateSegmentEnergyProfile

standardChromoConditions = _biolccc.standardChromoConditions
rpAcnTfaChain = _biolccc.rpAcnTfaChain
rpAcnFaRod = _biolccc.rpAcnFaRod


def _chemical_basis_getstate(self):
    state_dict = {}
    for key in self:
        if key != "chemicalGroups":
            state_dict[key] = self[key]
    state_dict["chemicalGroups"] = {}
    for key in self["chemicalGroups"]:
        state_dict["chemicalGroups"][key] = dict(self["chemicalGroups"][key])
    return state_dict


def _chemical_basis_setstate(self, state_dict):
    for key in state_dict:
        self[key] = state_dict[key]


def _chemical_basis_reduce(self):
    return (ChemicalBasis, (), self.__getstate__())


def _chromo_conditions_getstate(self):
    state_dict = {}
    for key in self:
        if key != "gradient":
            state_dict[key] = self[key]
    state_dict["gradient"] = []
    for point in self["gradient"]:
        state_dict["gradient"].append(dict(point))
    return state_dict


def _chromo_conditions_setstate(self, state_dict):
    for key in state_dict:
        if key != "gradient":
            self[key] = state_dict[key]

    gradient = Gradient()
    for point in state_dict.get("gradient", []):
        gradient.addPoint(point["time"], point["concentrationB"])
    self["gradient"] = gradient


def _chromo_conditions_reduce(self):
    return (ChromoConditions, (), self.__getstate__())


def _chromo_conditions_update(self, mapping):
    for key, value in mapping.items():
        self[key] = value


ChemicalBasis.__getstate__ = _chemical_basis_getstate
ChemicalBasis.__setstate__ = _chemical_basis_setstate
ChemicalBasis.__reduce__ = _chemical_basis_reduce

ChromoConditions.__getstate__ = _chromo_conditions_getstate
ChromoConditions.__setstate__ = _chromo_conditions_setstate
ChromoConditions.__reduce__ = _chromo_conditions_reduce
ChromoConditions.update = _chromo_conditions_update

# Preserve historical isinstance checks with explicit MutableMapping registration.
MutableMapping.register(ChemicalGroup)
MutableMapping.register(ChemicalBasis)
MutableMapping.register(GradientPoint)
MutableMapping.register(ChromoConditions)

__all__ = [
    "VERSION",
    "BioLCCCException",
    "ChemicalBasisException",
    "ChromoConditionsException",
    "GradientException",
    "GradientPointException",
    "ParsingException",
    "ROD",
    "CHAIN",
    "RP_ACN_TFA_CHAIN",
    "RP_ACN_FA_ROD",
    "ChemicalGroup",
    "ChemicalBasis",
    "GradientPoint",
    "Gradient",
    "ChromoConditions",
    "calculateRT",
    "calculateKd",
    "calculateAverageMass",
    "calculateMonoisotopicMass",
    "parseSequence",
    "calculateMonomerEnergyProfile",
    "calculateSegmentEnergyProfile",
    "standardChromoConditions",
    "rpAcnTfaChain",
    "rpAcnFaRod",
]
