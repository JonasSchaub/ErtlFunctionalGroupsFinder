/*
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (c) 2024 Sebastian Fritsch, Stefan Neumann, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 * 
 * Source code is available at <https://github.com/JonasSchaub/ErtlFunctionalGroupsFinder>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.openscience.cdk.tools;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Tests functionalities of ErtlFunctionalGroupsFinderUtility class.
 *
 * @author Jonas Schaub
 * @version 1.2
 */
public class ErtlFunctionalGroupsFinderUtilityTest {
    /**
     * Tests for correct pseudo SMILES generation.
     *
     * @throws Exception if a SMILES code cannot be parsed into a molecule or a pseudo SMILES string cannot be created
     */
    @Test
    public void testPseudoSmilesGeneration() throws Exception {
        HashMap<String, String> tmpTestPairsMap = new HashMap<>(20, 1);
        tmpTestPairsMap.put("*n(*)*", "RN*(R)R");
        tmpTestPairsMap.put("*O*", "ROR");
        tmpTestPairsMap.put("[H]O[c]", "[H]O[C*]");
        tmpTestPairsMap.put("*OC(*)=O", "ROC(R)=O");
        tmpTestPairsMap.put("*o*", "RO*R");
        tmpTestPairsMap.put("[c]=O", "[C*]=O");
        tmpTestPairsMap.put("[As]", "[As]");
        tmpTestPairsMap.put("[Po]", "[Po]");
        //The CDK SmilesParser cannot parse the element Uup, it gets turned into a wildcard ('*')
        tmpTestPairsMap.put("[Uup]", "R");
        tmpTestPairsMap.put("[se]", "[Se*]");
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpSmilesParser.kekulise(false);
        IAtomContainer tmpTestMolecule;
        String tmpPseudoSmilesCode;
        for (String tmpSmilesCode : tmpTestPairsMap.keySet()) {
            tmpTestMolecule = tmpSmilesParser.parseSmiles(tmpSmilesCode);
            tmpPseudoSmilesCode = ErtlFunctionalGroupsFinderUtility.createPseudoSmilesCode(tmpTestMolecule);
            Assertions.assertEquals(tmpTestPairsMap.get(tmpSmilesCode), tmpPseudoSmilesCode);
        }
    }
    //
    /**
     * Test for correct MoleculeHashGenerator settings/performance on some examples.
     *
     * @throws Exception if a SMILES code cannot be parsed into a molecule
     */
    @Test
    public void testMoleculeHashGeneratorSettings() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ErtlFunctionalGroupsFinder tmpGeneralizingEFGF = ErtlFunctionalGroupsFinder.newErtlFunctionalGroupsFinderGeneralizingMode();
        MoleculeHashGenerator tmpHashGenerator = ErtlFunctionalGroupsFinderUtility.getFunctionalGroupHashGenerator();
        /*Chebi70986, Chebi16238 and Chebi57692 all contain the same functional group with pseudo SMILES code
        "O=C1N=C(C(=NR)C(=O)N1R)N(R)R", but different hybridizations in the resulting atom containers. But their hash
        codes should be the same under the given settings. This is tested exemplary for many similar cases*/
        String[] tmpSmilesArray = {"OC[C@@H](O)[C@@H](O)[C@@H](O)CN1CC(CO)N=C2C(=O)NC(=O)N=C12",
                "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C",
                "Cc1cc2nc3c(nc(=O)[n-]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP([O-])(=O)OP([O-])(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C"};
        List<Long> tmpHashCodesList = new LinkedList<>();
        for (String tmpSmilesCode : tmpSmilesArray) {
            IAtomContainer tmpParsedMolecule = tmpSmilesParser.parseSmiles(tmpSmilesCode);
            tmpParsedMolecule = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(tmpParsedMolecule, Aromaticity.cdkLegacy());
            List<IAtomContainer> tmpFunctionalGroups = tmpGeneralizingEFGF.find(tmpParsedMolecule);
            for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroups) {
                if (ErtlFunctionalGroupsFinderUtility.createPseudoSmilesCode(tmpFunctionalGroup).equals("O=C1N=C(C(=NR)C(=O)N1R)N(R)R")) {
                    tmpHashCodesList.add(tmpHashGenerator.generate(tmpFunctionalGroup));
                }
            }
        }
        for (Long tmpHashCode1 : tmpHashCodesList) {
            for (Long tmpHashCode2 : tmpHashCodesList) {
                Assertions.assertEquals(tmpHashCode1.longValue(), tmpHashCode2.longValue());
            }
        }
        /*Functional groups like the tertiary amine or the hydroxyl group appear with aromatic and non-aromatic central
        atoms. These two cases should be discriminated by the MoleculeHashGenerator under the given settings*/
        String tmpTertiaryAmineSmiles = "*N(*)*";
        IAtomContainer tmpAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        IAtomContainer tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("N"))
                tmpAtom.setIsAromatic(true);
        }
        Assertions.assertNotEquals(tmpHashGenerator.generate(tmpAromMol), tmpHashGenerator.generate(tmpNonAromMol));
        String tmpHydroxylGroupSmiles = "[H]O[C]";
        tmpAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("C"))
                tmpAtom.setIsAromatic(true);
        }
        Assertions.assertNotEquals(tmpHashGenerator.generate(tmpAromMol), tmpHashGenerator.generate(tmpNonAromMol));
        /*The following are examples of different (unique!) SMILES codes representing the same functional groups.
        They should be assigned the same hash code*/
        HashMap<String,String> tmpEquivalentSmilesMap = new HashMap<>(20);
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)*", "*N(*)C(=[N](*)*)N(*)*");
        tmpEquivalentSmilesMap.put("*SC1=[N](*)[C]=[C]N1*", "*SC=1N(*)[C]=[C][N]1*");
        tmpEquivalentSmilesMap.put("*[N]1=[C][C]=[C]N1*", "*N1[C]=[C][C]=[N]1*");
        tmpEquivalentSmilesMap.put("*[N](*)=[C]N(*)*", "*N(*)[C]=[N](*)*");
        tmpEquivalentSmilesMap.put("*N(*)[C]=[C][C]=[C][C]=[C][C]=[C][C]=[N](*)*", "*[N](*)=[C][C]=[C][C]=[C][C]=[C][C]=[C]N(*)*");
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)P(=O)(O[H])O[H]", "*N(*)C(=[N](*)*)N(*)P(=O)(O[H])O[H]");
        tmpEquivalentSmilesMap.put("[O]I(=O)=O", "O=I(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Br(=O)=O", "O=Br(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Cl(=O)(=O)=O", "O=Cl(=O)(=O)[O]");
        tmpEquivalentSmilesMap.put("[C]=[C][C]=[C]C#C[C]=[C]C#[C]", "[C]#C[C]=[C]C#C[C]=[C][C]=[C]");
        tmpEquivalentSmilesMap.put("*N1[C]=[C][C]=[N]1*", "*[N]1=[C][C]=[C]N1*");
        tmpEquivalentSmilesMap.put("O=C(*)O*", "*OC(*)=O");
        for (String tmpKeySmiles : tmpEquivalentSmilesMap.keySet()) {
            IAtomContainer tmpKeyMol = tmpSmilesParser.parseSmiles(tmpKeySmiles);
            IAtomContainer tmpValueMol = tmpSmilesParser.parseSmiles(tmpEquivalentSmilesMap.get(tmpKeySmiles));
            Assertions.assertEquals(tmpHashGenerator.generate(tmpKeyMol), tmpHashGenerator.generate(tmpValueMol));
        }
    }
    //
    /**
     * Test for correct preprocessing (neutralization of charges and selection of biggest fragment).
     *
     * @throws Exception if a SMILES code can not be parsed into a molecule
     */
    @Test
    public void testPreprocessing() throws Exception {
        String tmpSmiles = "CC[O-].C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMol = tmpSmilesParser.parseSmiles(tmpSmiles);
        Assertions.assertTrue(ErtlFunctionalGroupsFinderUtility.shouldBePreprocessed(tmpMol));
        Assertions.assertFalse(ErtlFunctionalGroupsFinderUtility.shouldBeFiltered(tmpMol));
        tmpMol = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(tmpMol, Aromaticity.cdkLegacy());
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        Assertions.assertEquals("OCC", tmpGenerator.create(tmpMol));
    }
    //
    /**
     * Tests the restoration of environmental carbon atom objects on one example molecule. Nothing is asserted here, it
     * is meant for visual inspection.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testRestorationOfEnvironmentalCarbons() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        //Adenophostin B, COCONUT ID CNP0214672
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("O=C(OCC1OC(OC2C(OC(N3C=NC=4C(=NC=NC43)N)C2OP(=O)(O)O)CO)C(O)C(OP(=O)(O)O)C1OP(=O)(O)O)C");
        ErtlFunctionalGroupsFinder tmpEFGFFullEnv = ErtlFunctionalGroupsFinder.newErtlFunctionalGroupsFinderFullEnvironmentMode();
        tmpMolecule = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(tmpMolecule, Aromaticity.cdkLegacy());
        List<IAtomContainer> tmpFGList = tmpEFGFFullEnv.find(tmpMolecule, false);
        System.out.println("FGs with full environment returned by EFGF:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
        ErtlFunctionalGroupsFinderUtility.restoreOriginalEnvironmentalCarbons(tmpFGList, tmpMolecule, false, false, SilentChemObjectBuilder.getInstance());
        System.out.println("FGs with full environment, environmental carbons restored:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
        ErtlFunctionalGroupsFinderUtility.restoreOriginalEnvironmentalCarbons(tmpFGList, tmpMolecule, true, true, SilentChemObjectBuilder.getInstance());
        System.out.println("FGs with full environment, environmental carbons restored, explicit hydrogens removed and empty valences filled:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
        tmpMolecule = tmpSmiPar.parseSmiles("O=C(OCC1OC(OC2C(OC(N3C=NC=4C(=NC=NC43)N)C2OP(=O)(O)O)CO)C(O)C(OP(=O)(O)O)C1OP(=O)(O)O)C");
        ErtlFunctionalGroupsFinder tmpEFGFgeneralized = ErtlFunctionalGroupsFinder.newErtlFunctionalGroupsFinderGeneralizingMode();
        tmpMolecule = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(tmpMolecule, Aromaticity.cdkLegacy());
        tmpFGList = tmpEFGFgeneralized.find(tmpMolecule, false);
        System.out.println("FGs with generalized environment returned by EFGF:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
        ErtlFunctionalGroupsFinderUtility.restoreOriginalEnvironmentalCarbons(tmpFGList, tmpMolecule, false, false, SilentChemObjectBuilder.getInstance());
        System.out.println("FGs with generalized environment, environmental carbons restored:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
        ErtlFunctionalGroupsFinderUtility.restoreOriginalEnvironmentalCarbons(tmpFGList, tmpMolecule, true, true, SilentChemObjectBuilder.getInstance());
        System.out.println("FGs with generalized environment, environmental carbons restored, explicit hydrogens removed and empty valences filled:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
    }
    //
    /**
     * Imports a charged molecule with a counter-ion from ChEMBL to test the filtering and preprocessing routines
     * of ErtlFunctionalGroupsFinderUtility.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void testOnMolecule() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        //CHEMBL1201736
        IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles("CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(C[N+]3(C)CCCC3)CS[C@H]12)c1csc(N)n1.Cl");
        Assertions.assertTrue(ErtlFunctionalGroupsFinder.isStructureUnconnected(tmpMolecule));
        Assertions.assertTrue(ErtlFunctionalGroupsFinder.containsChargedAtom(tmpMolecule));
        Assertions.assertFalse(ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(tmpMolecule));
        Assertions.assertFalse(ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(tmpMolecule));
        Assertions.assertFalse(ErtlFunctionalGroupsFinderUtility.shouldBeFiltered(tmpMolecule));
        Assertions.assertTrue(ErtlFunctionalGroupsFinderUtility.shouldBePreprocessed(tmpMolecule));
        Assertions.assertFalse(ErtlFunctionalGroupsFinderUtility.isValidArgumentForFindMethod(tmpMolecule));
        tmpMolecule = ErtlFunctionalGroupsFinderUtility.selectBiggestUnconnectedComponent(tmpMolecule);
        Assertions.assertNotNull(tmpMolecule);
        ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpMolecule);
        ErtlFunctionalGroupsFinderUtility.perceiveAtomTypesAndConfigureAtoms(tmpMolecule);
        ErtlFunctionalGroupsFinderUtility.applyAromaticityDetection(tmpMolecule, Aromaticity.cdkLegacy());
        Assertions.assertTrue(ErtlFunctionalGroupsFinderUtility.isValidArgumentForFindMethod(tmpMolecule));
        ErtlFunctionalGroupsFinder tmpEFGF = ErtlFunctionalGroupsFinder.newErtlFunctionalGroupsFinderGeneralizingMode();
        List<IAtomContainer> tmpFGList = tmpEFGF.find(tmpMolecule);
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(ErtlFunctionalGroupsFinderUtility.createPseudoSmilesCode(tmpFG));
        }
    }
    //
    /**
     * Test charge neutralization.
     */
    @Test
    public void testNeutralization() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpAmmonia = tmpSmiPar.parseSmiles("[NH4+]");
        ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpAmmonia);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical);
        System.out.println(tmpSmiGen.create(tmpAmmonia));
        IAtomContainer tmpNitro = tmpSmiPar.parseSmiles("C[N+](=O)[O-]");
        ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpNitro);
        System.out.println(tmpSmiGen.create(tmpNitro));
    }
}
