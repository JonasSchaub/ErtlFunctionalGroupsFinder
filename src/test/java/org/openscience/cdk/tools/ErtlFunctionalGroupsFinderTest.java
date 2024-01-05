/*
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (c) 2023 Sebastian Fritsch, Stefan Neumann, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Test for ErtlFunctionalGroupsFinder.
 *
 * @author Sebastian Fritsch
 * @version 1.2
 */
public class ErtlFunctionalGroupsFinderTest {

    public ErtlFunctionalGroupsFinderTest() {
        super();
    }

    @Test
    public void testFind1() throws Exception {
        String moleculeSmiles = "Cc1cc(C)nc(NS(=O)(=O)c2ccc(N)cc2)n1";
        String[] expectedFGs = new String[] {"[R]N([R])S(=O)(=O)[R]", "[c]N(H)H", "NarR3", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind2() throws Exception{
        String moleculeSmiles = "NC(=N)c1ccc(\\\\C=C\\\\c2ccc(cc2O)C(=N)N)cc1";
        String[] expectedFGs = new String[] {"[R]N=C-N([R])[R]", "[C]=[C]", "[c]OH", "[R]N=C-N([R])[R]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind3() throws Exception {
        String moleculeSmiles = "CC(=O)Nc1nnc(s1)S(=O)(=O)N";
        String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "[R]S(=O)(=O)N([R])[R]", "NarR3", "NarR3", "SarR2"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind4() throws Exception {
        String moleculeSmiles = "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl";
        String[] expectedFGs = new String[] {"[R]S(=O)(=O)N([R])[R]", "[R]S(=O)(=O)N([R])[C]N([R])[R]", "[R]Cl"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind5() throws Exception {
        String moleculeSmiles = "CNC1=Nc2ccc(Cl)cc2C(=N(=O)C1)c3ccccc3";
        String[] expectedFGs = new String[] {"[R]N([R])[C]=N[R]", "[R]Cl", "[R]N(=O)=[C]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind6() throws Exception {
        String moleculeSmiles = "Cc1onc(c2ccccc2)c1C(=O)N[C@H]3[C@H]4SC(C)(C)[C@@H](N4C3=O)C(=O)O";
        String[] expectedFGs = new String[] {"O=C([R])N([R])[R]",  "O=C([R])N([R])[C]S[R]", "O=C([R])OH", "OarR2", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind7() throws Exception {
        String moleculeSmiles = "Clc1ccccc1C2=NCC(=O)Nc3ccc(cc23)N(=O)=O";
        String[] expectedFGs = new String[] {"[R]Cl", "[R]N=[C]", "[R]C(=O)N([R])[R]", "O=N([R])=O"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind8() throws Exception {
        String moleculeSmiles = "COc1cc(cc(C(=O)NCC2CCCN2CC=C)c1OC)S(=O)(=O)N";
        String[] expectedFGs = new String[] {"[R]O[R]", "[R]N([R])C(=O)[R]", "N([R])([R])[R]", "[C]=[C]", "[R]O[R]", "[R]S(=O)(=O)N([R])[R]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind9() throws Exception {
        String moleculeSmiles = "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl";
        String[] expectedFGs = new String[] {"[R]Cl", "[R]N(H)[R]", "O=C(OH)[R]", "[R]Cl"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind10() throws Exception {
        String moleculeSmiles = "Clc1ccc2Oc3ccccc3N=C(N4CCNCC4)c2c1";
        String[] expectedFGs = new String[] {"[R]Cl", "[R]O[R]", "[R]N([R])[C]=N[R]", "[R]N([H])[R]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind11() throws Exception {
        String moleculeSmiles = "FC(F)(F)CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13";
        String[] expectedFGs = new String[] {"[R]F", "[R]F", "[R]F", "O=C([R])N([R])[R]", "[R]N=[C]", "[R]Cl"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind12() throws Exception {
        String moleculeSmiles = "OC[C@H]1O[C@H](C[C@@H]1O)n2cnc3[C@H](O)CNC=Nc23";;
        String[] expectedFGs = new String[] {"[C]O[H]", "[R]O[R]", "[C]OH", "[C]OH", "[R]N=CN([R])[R]", "NarR3", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind13() throws Exception {
        String moleculeSmiles = "CCN[C@H]1C[C@H](C)S(=O)(=O)c2sc(cc12)S(=O)(=O)N";
        String[] expectedFGs = new String[] {"[R]N([R])H", "O=S(=O)([R])[R]", "[R]S(=O)(=O)N([R])[R]", "SarR2"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind14() throws Exception {
        String moleculeSmiles = "C[C@@H](O)[C@@H]1[C@H]2[C@@H](C)C(=C(N2C1=O)C(=O)O)S[C@@H]3CN[C@@H](C3)C(=O)N(C)C";
        String[] expectedFGs = new String[] {"[C]O[H]", "O=C([R])N([R])C(C(=O)(OH))=[C]S[R]", "[R]N(H)[R]", "[R]N([R])C([R])=O"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind15() throws Exception {
        String moleculeSmiles = "C[C@@H]1CN(C[C@H](C)N1)c2c(F)c(N)c3C(=O)C(=CN(C4CC4)c3c2F)C(=O)O";
        String[] expectedFGs = new String[] {"[R]N([R])[R]", "[R]N([H])[R]", "[R]F", "[c]N(H)H", "[c]=O", "[R]F", "[R]C(=O)OH", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind16() throws Exception {
        String moleculeSmiles = "CC(=CCC1C(=O)N(N(C1=O)c2ccccc2)c3ccccc3)C";
        String[] expectedFGs = new String[] {"[C]=[C]", "[R]C(=O)N([R])N([R])C(=O)[R]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind17() throws Exception {
        String moleculeSmiles = "Clc1ccc2N=C3NC(=O)CN3Cc2c1Cl";
        String[] expectedFGs = new String[] {"Cl[R]", "[R]N=C(N([R])[R])N([R])C(=O)[R]", "Cl[R]"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind18() throws Exception {
        String moleculeSmiles = "CC(=O)N[C@@H]1[C@@H](NC(=N)N)C=C(O[C@H]1[C@H](O)[C@H](O)CO)C(=O)O";
        String[] expectedFGs = new String[] {"[R]N([R])C(=O)[R]", "[R]N([R])C(=N[R])N([R])[R]", "O=C(OH)C(=[C])O[R]" , "[C]OH", "[C]OH", "[C]OH"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind19() throws Exception {
        String moleculeSmiles = "C[C@H](O)[C@H](O)[C@H]1CNc2nc(N)nc(O)c2N1";
        String[] expectedFGs = new String[] {"[C]OH", "[C]OH", "[R]N(H)[R]" , "[c]N(H)H",  "[c]OH", "[R]N(H)[R]", "NarR3", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    @Test
    public void testFind20() throws Exception {
        String moleculeSmiles = "N[C@@H]1CCCCN(C1)c2c(Cl)cc3C(=O)C(=CN(C4CC4)c3c2Cl)C(=O)O";
        String[] expectedFGs = new String[] {"[C]N([H])[H]", "[R]N([R])[R]", "[R]Cl" , "[c]=O", "[R]Cl", "[R]C(=O)OH", "NarR3"};
        testFind(moleculeSmiles, expectedFGs);
    }

    /**
     * Example code to be used in the GitHub wiki of the project.
     *
     * @throws Exception if anything goes wrong
     * @author Jonas Schaub
     */
    @Test
    public void gitHubWikiTest() throws Exception {
        //Prepare input
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpInputMol = tmpSmiPar.parseSmiles("C[C@@H]1CN(C[C@H](C)N1)C2=C(C(=C3C(=C2F)N(C=C(C3=O)C(=O)O)C4CC4)N)F"); //PubChem CID 5257
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpInputMol);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
        tmpAromaticity.apply(tmpInputMol);
        //Identify functional groups
        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(); //default: generalization turned on
        List<IAtomContainer> tmpFunctionalGroupsList = tmpEFGF.find(tmpInputMol);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.UseAromaticSymbols);
        for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroupsList) {
            String tmpSmilesString = tmpSmiGen.create(tmpFunctionalGroup);
            System.out.println(tmpSmilesString);
        }
        //non-generalized functional groups
        System.out.println("----------------");
        tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.NO_GENERALIZATION);
        tmpFunctionalGroupsList = tmpEFGF.find(tmpInputMol);
        for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroupsList) {
            String tmpSmilesString = tmpSmiGen.create(tmpFunctionalGroup);
            System.out.println(tmpSmilesString);
        }
    }

    /**
     * TODO: Investigate code for possible problems with charged atoms?
     *
     * TODO: Test carbon ions.
     *
     * @throws Exception
     */
    @Test
    public void testChargedMolecules() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.UseAromaticSymbols);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

        IAtomContainer tmpChargedASA = tmpSmiPar.parseSmiles("CC(=O)OC1=CC=CC=C1C(=O)[O+]");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpChargedASA);
        tmpAromaticity.apply(tmpChargedASA);

        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.DEFAULT);
        List<IAtomContainer> tmpFGList = tmpEFGF.find(tmpChargedASA);

        System.out.println("Charged ASA:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }

        IAtomContainer tmpNitroPhenol = tmpSmiPar.parseSmiles("C1=CC(=CC=C1[N+](=O)[O-])O");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpNitroPhenol);
        tmpAromaticity.apply(tmpNitroPhenol);

        tmpFGList = tmpEFGF.find(tmpNitroPhenol);

        System.out.println("Nitrophenol:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
    }

    /**
     * TODO: Investigate code for possible problems with disconnected structures?
     *
     * @throws Exception
     */
    @Test
    public void testDisconnectedMolecules() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.UseAromaticSymbols);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

        IAtomContainer tmpChlorhexidineDiacetate = tmpSmiPar.parseSmiles("CC(=O)O.CC(=O)O.C1=CC(=CC=C1NC(=NC(=NCCCCCCN=C(N)N=C(N)NC2=CC=C(C=C2)Cl)N)N)Cl");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpChlorhexidineDiacetate);
        tmpAromaticity.apply(tmpChlorhexidineDiacetate);

        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.DEFAULT);
        List<IAtomContainer> tmpFGList = tmpEFGF.find(tmpChlorhexidineDiacetate);

        System.out.println("Chlorhexidine Diacetate:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }

        IAtomContainer tmpSodiumEdetate = tmpSmiPar.parseSmiles("C(CN(CC(=O)[O-])CC(=O)[O-])N(CC(=O)[O-])CC(=O)[O-].[Na+].[Na+].[Na+].[Na+]");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpSodiumEdetate);
        tmpAromaticity.apply(tmpSodiumEdetate);

        tmpFGList = tmpEFGF.find(tmpSodiumEdetate);

        System.out.println("Sodium edetate:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }
    }

    /**
     *
     *
     * Note: all atoms are marked as hetero atoms by EFGF that are not H or C. So, metals and metalloids get treated like
     * any other hetero atom and should not cause problems.
     *
     * @throws Exception
     */
    @Test
    public void testMetalsMetalloids() throws Exception {
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.UseAromaticSymbols);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

        IAtomContainer tmpTetraethylOrthosilicate = tmpSmiPar.parseSmiles("CCO[Si](OCC)(OCC)OCC");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpTetraethylOrthosilicate);
        tmpAromaticity.apply(tmpTetraethylOrthosilicate);

        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.DEFAULT);
        List<IAtomContainer> tmpFGList = tmpEFGF.find(tmpTetraethylOrthosilicate);

        System.out.println("Tetraethyl Orthosilicate:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }

        IAtomContainer tmpKaolin = tmpSmiPar.parseSmiles("O.O.O=[Al]O[Si](=O)O[Si](=O)O[Al]=O");
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpKaolin);
        tmpAromaticity.apply(tmpKaolin);

        tmpFGList = tmpEFGF.find(tmpKaolin);

        System.out.println("Kaolin:");
        for (IAtomContainer tmpFG : tmpFGList) {
            System.out.println(tmpSmiGen.create(tmpFG));
        }

    }

    //TODO: Clean-up check constraints and add test molecules for these special cases to the testFind#() methods.

    /**
     * TODO: test complete ChEBI?
     *
     * Note: ChEBI lite 3-star subset SDF contains 251 molecules with charges or metal/metalloid atoms or more than one
     * disconnected structure (comment-in checkConstraints in EFGF.find() method to check).
     *
     * @throws Exception
     */
    @Test
    public void readChebiLite3StarSubset() throws Exception {
        IteratingSDFReader tmpChebiSDFReader = new IteratingSDFReader(
                ErtlFunctionalGroupsFinderTest.class.getResourceAsStream("ChEBI_lite_3star_subset.sdf"),
                SilentChemObjectBuilder.getInstance(),
                false);
        Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.DEFAULT);
        int tmpMoleculeCouter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpChebiSDFReader.hasNext()) {
            IAtomContainer tmpMolecule = null;
            tmpMoleculeCouter++;
            try {
                tmpMolecule = tmpChebiSDFReader.next();
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
                tmpAromaticity.apply(tmpMolecule);

                List<IAtomContainer> tmpFGList = tmpEFGF.find(tmpMolecule);
            } catch (Exception anException) {
                tmpExceptionsCounter++;
                if (!Objects.isNull(tmpMolecule)) {
                    System.out.println(tmpMolecule.getProperty("ChEBI ID") + "," + anException.toString() + "," + tmpMoleculeCouter);
                } else {
                    System.out.println("Could not parse molecule! Counter: " + tmpMoleculeCouter);
                }
            }

        }
        System.out.println(tmpMoleculeCouter);
        System.out.println(tmpExceptionsCounter);
    }

    private void testFind(String moleculeSmiles, String[] fGStrings) throws Exception {
        testFind(moleculeSmiles, fGStrings, new Aromaticity(ElectronDonation.daylight(), Cycles.all()));
    }

    private void testFind(String moleculeSmiles, String[] fGStrings, Aromaticity aromaticity) throws Exception {
        // prepare input
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer mol = smilesParser.parseSmiles(moleculeSmiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        aromaticity.apply(mol);

        // find functional groups
        ErtlFunctionalGroupsFinder fgFinder = new ErtlFunctionalGroupsFinder();
        List<IAtomContainer> fGs = fgFinder.find(mol);

        // get expected groups
        List<IAtomContainer> expectedFGs = new LinkedList<>();
        for (String fGString : fGStrings) {
            expectedFGs.add(buildFunctionalGroup(fGString));
        }

        // compare
        this.assertIsomorphism(expectedFGs, fGs);
    }

    /**
     * NOTE: actual and expected functional groups must be in the same order!
     *
     * @param expectedFGs 	list of expected functional groups
     * @param actualFGs			list of actual functional groups
     * @throws Exception	if anything does not work as planned
     */
    private void assertIsomorphism(List<IAtomContainer> expectedFGs, List<IAtomContainer> actualFGs) {
        Assertions.assertEquals(expectedFGs.size(), actualFGs.size(),
                "Number of functional groups does not match the expected number of groups");

        for(int i = 0; i < expectedFGs.size(); i++) {
            IAtomContainer cExp = expectedFGs.get(i);
            IAtomContainer cAct = actualFGs.get(i);

            Assertions.assertEquals(cExp.getAtomCount(), cAct.getAtomCount(),
                    "Groups #" + i + ": different atom count");
            Assertions.assertEquals(cExp.getBondCount(),  cAct.getBondCount(),
                    "Groups #" + i + ": different bond count");

            Pattern pattern = VentoFoggia.findIdentical(cExp);

            Assertions.assertTrue(pattern.matches(cAct), "Groups #" + i + ": not isomorph");

            Mappings mappings = pattern.matchAll(cAct);

            Map<IAtom, IAtom> atomMap = mappings.toAtomMap().iterator().next();
            for (Map.Entry<IAtom, IAtom> e : atomMap.entrySet()) {
                 IAtom atomExp  = e.getKey();
                 IAtom atomAct = e.getValue();
                 Assertions.assertEquals(atomExp.isAromatic(), atomAct.isAromatic(),
                         "Groups #" + i + ": Atom aromaticity does not match"
                                 + atomAct.getSymbol() + atomAct.isAromatic() + atomExp.getSymbol()
                                 + atomExp.isAromatic());
             }

            Map<IBond, IBond> bondMap = mappings.toBondMap().iterator().next();
            for (Map.Entry<IBond, IBond> e : bondMap.entrySet()) {
                 IBond bondExp  = e.getKey();
                 IBond bondAct = e.getValue();
                 Assertions.assertEquals(bondExp.isAromatic(), bondAct.isAromatic(),
                         "Groups #" + i + ": Bond aromaticity does not match");
             }
        }
    }

    private IAtomContainer buildFunctionalGroup(String string) {
        IAtom a1, a2, a3, a4, a5, a6, a7, a8, a9;
        IBond b1, b2, b3, b4, b5, b6, b7, b8, b9;
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        IAtomContainer container;

        // custom templates
        switch(string) {
        case "NarR3":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IPseudoAtom.class, "R");
            a4 = builder.newInstance(IAtom.class, "N");
            a4.setIsAromatic(true);

            b1 = builder.newInstance(IBond.class, a1, a4, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a4, Order.SINGLE);
            b3 = builder.newInstance(IBond.class, a3, a4, Order.SINGLE);

            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3, a4});
            container.setBonds(new IBond[] {b1, b2, b3});
            return container;

        case "SarR2":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "S");
            a3.setIsAromatic(true);

            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);

            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;

        case "OarR2":
            a1 = builder.newInstance(IPseudoAtom.class, "R");
            a2 = builder.newInstance(IPseudoAtom.class, "R");
            a3 = builder.newInstance(IAtom.class, "O");
            a3.setIsAromatic(true);

            b1 = builder.newInstance(IBond.class, a1, a3, Order.SINGLE);
            b2 = builder.newInstance(IBond.class, a2, a3, Order.SINGLE);

            container = new AtomContainer();
            container.setAtoms(new IAtom[] {a1, a2, a3});
            container.setBonds(new IBond[] {b1, b2});
            return container;

            // smiles
        default:
            try {
                SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
                try {
                    if(string.equals("[c]=O"))
                        smilesParser.kekulise(false);
                    container = smilesParser.parseSmiles(string);
                }
                catch(InvalidSmilesException e) {
                    smilesParser.kekulise(false);
                    container = smilesParser.parseSmiles(string);
                }

                for(IAtom a : container.atoms()) {
                    if(a instanceof PseudoAtom) {
                        a.setSymbol("R");
                    }
                }
                return container;
            }
            catch(InvalidSmilesException e) {
                throw new IllegalArgumentException("Input string '" + string + " could not be found as a template " +
                        "and is not a valid SMILES string.");
            }
        }
    }
}
