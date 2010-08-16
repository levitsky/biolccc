#include "biolccc.h"
#include <gtest/gtest.h>
#include <vector>

// The fixture for testing BioLCCC functions.
class BioLCCCTest: public ::testing::Test {
 protected:
  BioLCCCTest() {
    // You can do set-up work for each test here.
  }

  virtual ~BioLCCCTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }
};

TEST_F(BioLCCCTest, parsesStandardAminoacids) {
    std::vector<BioLCCC::ChemicalGroup> parsedPeptideStructure;
    BioLCCC::ChemicalGroup NTerminus;
    BioLCCC::ChemicalGroup CTerminus;
    std::vector<double> peptideEnergyProfile;
    std::string peptideSource("QWERTYIPASDFGHKLCVNM");

    ASSERT_TRUE(BioLCCC::parseSequence(peptideSource,
        BioLCCC::standardChemicalBasis, &parsedPeptideStructure, 
        &NTerminus, &CTerminus, &peptideEnergyProfile));

    for (unsigned int i=0; i<peptideSource.size(); i++) {
        ASSERT_EQ(parsedPeptideStructure[i].label(),
            peptideSource.substr(i, 1));
    }
}

TEST_F(BioLCCCTest, parsesPhosphoAminoacids) {
    std::vector<BioLCCC::ChemicalGroup> parsedPeptideStructure;
    BioLCCC::ChemicalGroup NTerminus;
    BioLCCC::ChemicalGroup CTerminus;
    std::vector<double> peptideEnergyProfile;
    std::string peptideSource("pSpTpY");

    ASSERT_TRUE(BioLCCC::parseSequence(peptideSource,
        BioLCCC::standardChemicalBasis, &parsedPeptideStructure, 
        &NTerminus, &CTerminus, &peptideEnergyProfile));

    for (unsigned int i=0; i<peptideSource.size()/2; i++) {
        ASSERT_EQ(parsedPeptideStructure[i].label(),
            peptideSource.substr(i*2, 2));
    }
}

TEST_F(BioLCCCTest, parsesStandardTerminalGroups) {
    std::vector<BioLCCC::ChemicalGroup> parsedPeptideStructure;
    BioLCCC::ChemicalGroup NTerminus;
    BioLCCC::ChemicalGroup CTerminus;
    std::vector<double> peptideEnergyProfile;

    ASSERT_TRUE(BioLCCC::parseSequence("GGGG",
        BioLCCC::standardChemicalBasis, &parsedPeptideStructure, 
        &NTerminus, &CTerminus, &peptideEnergyProfile));
    ASSERT_EQ(NTerminus.label(), "H-");
    ASSERT_EQ(CTerminus.label(), "-COOH");

    ASSERT_TRUE(BioLCCC::parseSequence("H-GGGG-COOH",
        BioLCCC::standardChemicalBasis, &parsedPeptideStructure, 
        &NTerminus, &CTerminus, &peptideEnergyProfile));
    ASSERT_EQ(NTerminus.label(), "H-");
    ASSERT_EQ(CTerminus.label(), "-COOH");

    ASSERT_TRUE(BioLCCC::parseSequence("Ac-GGGG-NH2",
        BioLCCC::standardChemicalBasis, &parsedPeptideStructure, 
        &NTerminus, &CTerminus, &peptideEnergyProfile));
    ASSERT_EQ(NTerminus.label(), "Ac-");
    ASSERT_EQ(CTerminus.label(), "-NH2");
}

TEST_F(BioLCCCTest, calculatesMonoisotopicMass) {
    ASSERT_LE(
        (BioLCCC::calculateMonoisotopicMass("QWERTYIPASDFGHKLCVNM")-2394.1248),
            0.0001);
    ASSERT_LE(
        (BioLCCC::calculateMonoisotopicMass("Ac-QWERTYIPASDFGHKLCVNM")-
            2436.1354), 0.0001);
    ASSERT_LE(
        (BioLCCC::calculateMonoisotopicMass("QWERTYIPASDFGHKLCVNM-NH2")-
            2393.1408), 0.0001);
}

TEST_F(BioLCCCTest, calculatesKd) {
    ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 0.0), 0.0);
    ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 50.0), 0.0);
    ASSERT_GT(BioLCCC::calculateKd("QWERTYIPASDFGHKLCVNM", 100.0), 0.0);
}

TEST_F(BioLCCCTest, calculatesRT) {
    ASSERT_GT(BioLCCC::calculateRT("QWERTYIPASDFGHKLCVNM"), 0.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
