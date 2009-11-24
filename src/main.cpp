/*! \mainpage TheorChromo 
 *
 * \section intro_sec Introduction
 *
 * TheorChromo is a part of new generation proteomic package. It is 
 * a command-line utility that calculates retention times of peptides in 
 * an HPLC experiment. 
 * \section change_log_sec Changelog.
 * \subsection one_point_zero_one 1.01
 * The project is started.
 * \section main_principles_sec Main principles of TheorChromo
 * - The utility has to calculate retention times of a peptide in an RP HPLC 
 * experiment. Incoming information are peptide sequences and experimental 
 * conditions. 
 * - It should not work with masses or do something alike. It should do the 
 * one thing and do it perfectly.
 * - It must be able to work in a pipeline, with peptide sequences in stdin 
 * and (sequence, RT) pairs in stdout.
 * - It must be as cross-compilable as possible. That's why the language of 
 * choice is a standard C++ and no external dependencies are currently 
 * involved.
 * \section internals_sec Internals of TheorChromo
 * \subsection conventions_subsec Conventions
 * - Main sequence format is 
 * "[X.][N-terminal group-]XXX[-C-terminal group][.X]". If terminal groups are
 * omitted, they are considered to be "H-" (hydrogen) and "-COOH" (free acid).
 * "X." and ".X" are adjacent aminoacids in a protein sequencs, and they are
 * always discarded.  Aminoacids code X could consist of few letters, but 
 * there should be no ambiguity. In the beginning of a calculation any 
 * ambiguity should be found, and in that case program should gave zeros 
 * instead of RTs.
 * - Default chromatographic conditions (standard Dionex conditions) and 
 * energies are incorporated into a program.
 * - Format of entries in .ini file: "property=x". Commenting character is #.
 * It's desirable to have no upper-case letter in a property name.
 * - Aminoacids overloaded or added by entry: "energy of M=0.14"
 * - Terminal groups should be overloaded by entry "energy of H-=-1.70"
 * - By default dV is calculated as: flow_rate / 10
 * \subsection getoptpp_subsec GetOpt_pp
 * Due to absence of getopt() function on Windows, we added one of its' free 
 * implementations called GetOpt_pp (http://code.google.com/p/getoptpp/). It 
 * consists of one header and one .cpp, which are located at the source 
 * directory.
 * \subsection code_format_subsec Code formatting conventions
 * Some formatting rules are adopted in this project. Not all of them are 
 * really rational, but anyway it's desirable to follow them. :)
 * - No "using" statements.
 * - 4 spaces instead of tabs.
 * - Inline comment should state on a new line with a blank line above.
 * - If a function has a lot of arguments they should be located each on a new
 * line with a standard 4-space indent and a closing bracket on a new line.
 * - Naming conventions:
 *     - ClassName
 *     - mPrivateVariable
 *     - privateVariable()
 *     - setPrivateVariable( iPrivateVariable )
 * - Doxygen with Qt style is used in this project.
 * \subsection roadmap_subsec Road map.
 * - to check, whether we could calculate a Kd starting from the first 
 * aminoacid rather than from the last.
 * - find the cause of different RT
 */
 
#include <iostream>

#include "BioLCCC.h"
#include "./boost/foreach.hpp"

int main(int argc, char* argv[]) {
    std::string name;

    BioLCCC::ChemicalBasis standardChemicalBasis;
    BioLCCC::ChromoConditions chromatograph;
    std::vector<std::string> peptideList;
    peptideList.push_back("KYIPGTK");
    peptideList.push_back("YIPGTK");
    peptideList.push_back("IFVQK");
    peptideList.push_back("KTGQAPGFSYTDANK");
    peptideList.push_back("TGQAPGFSYTDANK");
    peptideList.push_back("GEREDLIAYLKK");
    peptideList.push_back("TGPNLHGLFGR");
    peptideList.push_back("MIFAGIK");
    peptideList.push_back("EDLIAYLK");
    peptideList.push_back("IFVQKCAQCHTVEK");
    peptideList.push_back("GITWGEETLMEYLENPKK");
    peptideList.push_back("GITWGEETLMEYLENPK");
    BOOST_FOREACH(std::string peptide, peptideList) {
        std::cout << peptide << ": " << BioLCCC::calculateRT(peptide) << "\n";
    }

    return 1;
}

