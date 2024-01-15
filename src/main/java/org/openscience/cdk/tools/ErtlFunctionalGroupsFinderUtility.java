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

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.hash.AtomEncoder;
import org.openscience.cdk.hash.BasicAtomEncoder;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This class gives utility methods for using <a href="https://github.com/zielesny/ErtlFunctionalGroupsFinder">ErtlFunctionalGroupsFinder</a>,
 * a CDK-based implementation, published <a href="https://doi.org/10.1186/s13321-019-0361-8">here</a> of the
 * <a href="https://doi.org/10.1186/s13321-017-0225-z">Ertl algorithm for automated functional groups detection</a>.
 * The methods of this class are basically public static re-implementations of the routines used for testing and
 * evaluating the ErtlFunctionalGroupsFinder, as described in the publication.
 *
 * @author Jonas Schaub
 * @version 1.2
 */
public class ErtlFunctionalGroupsFinderUtility {
    //<editor-fold defaultstate="collapsed" desc="Enum CustomAtomEncoder">
    /**
     * Enumeration of custom atom encoders for seeding atomic hash codes.
     *
     * @author Jonas Schaub
     * @see BasicAtomEncoder
     * @see AtomEncoder
     */
    enum CustomAtomEncoder implements AtomEncoder {
        /**
         * Encode whether an atom is aromatic or not. This specification is necessary to distinguish functional groups with
         * aromatic environments and those without. For example: [H]O[C] and [H]OC* (pseudo SMILES codes) should be
         * assigned different hash codes by the MoleculeHashGenerator.
         *
         * @see IAtom#isAromatic()
         */
        AROMATICITY {
            /**
             *{@inheritDoc}
             */
            @Override
            public int encode(IAtom anAtom, IAtomContainer aContainer) {
                return anAtom.isAromatic()? 3 : 2;
            }
        };
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private static final class constants">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(ErtlFunctionalGroupsFinderUtility.class.getName());
    //</editor-fold>
    //
    private ErtlFunctionalGroupsFinderUtility() {

    }
    //
    //<editor-fold desc="Public static methods">
    //<editor-fold desc="Constants and new instances">
    /**
     * Constructs a CDK MoleculeHashGenerator that is configured to count frequencies of the functional groups
     * returned by ErtlFunctionalGroupsFinder. It takes elements, bond order sum, and aromaticity of the atoms in
     * an atom container into consideration. It does not consider things like isotopes, stereo-chemistry,
     * orbitals, or charges.
     *
     * @return MoleculeHashGenerator object configured for Ertl functional groups
     */
    public static MoleculeHashGenerator getFunctionalGroupHashGenerator() {
        MoleculeHashGenerator tmpHashGenerator = new HashGeneratorMaker()
                .depth(8)
                .elemental()
                /*following line is used instead of .orbital() because the atom hybridizations take more information into
                account than the bond order sum but that is not required here*/
                /*Note: This works here because the ErtlFunctionalGroupsFinder extracts the relevant atoms and bonds only
                resulting in incomplete valences that can be used here in this way*/
                .encode(BasicAtomEncoder.BOND_ORDER_SUM)
                .encode(CustomAtomEncoder.AROMATICITY) //See enum CustomAtomEncoder below
                .molecular();
        return tmpHashGenerator;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Queries for filtering">
    /**
     * Checks whether the atom count or bond count of the given molecule is zero. The ErtlFunctionalGroupsFinder.find()
     * method would still accept these molecules, but it is not recommended to pass them on (simply makes not much sense).
     *
     * @param aMolecule the molecule to check
     * @return true, if the atom or bond count of the molecule is zero
     * @throws NullPointerException if the given molecule is 'null'
     */
    public static boolean isAtomOrBondCountZero(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpAtomCount = aMolecule.getAtomCount();
        int tmpBondCount = aMolecule.getBondCount();
        return (tmpAtomCount == 0 || tmpBondCount == 0);
    }

    /**
     * Checks whether the given molecule represented by an atom container should NOT be passed on to the
     * ErtlFunctionalGroupsFinder.find() method but instead be filtered.
     * <br>In detail, this function returns true if the given atom container contains metal, metalloid, or pseudo atoms
     * or has an atom or bond count equal to zero.
     * <br>If this method returns false, this does NOT mean the molecule can be passed on to find() without a problem. It
     * still might need to be preprocessed first.
     *
     * @param aMolecule the atom container to check
     * @return true if the given atom container should be discarded
     * @throws NullPointerException if parameter is 'null'
     */
    public static boolean shouldBeFiltered(IAtomContainer aMolecule) throws NullPointerException {
        return ErtlFunctionalGroupsFinderUtility.shouldBeFiltered(aMolecule, true);
    }

    /**
     * Checks whether the given molecule represented by an atom container should NOT be passed on to the
     * ErtlFunctionalGroupsFinder.find() method but instead be filtered.
     * <br>In detail, this function returns true if the given atom container contains metal, metalloid, or pseudo atoms
     * or has an atom or bond count equal to zero. If the second parameter is set to "false", single atom molecules
     * (bond count is 0) are accepted and not recommended to be filtered if they fulfill the other requirements.
     * <br>If this method returns false, this does NOT mean the molecule can be passed on to find() without a problem. It
     * still might need to be preprocessed first.
     *
     * @param aMolecule the atom container to check
     * @param areSingleAtomsFiltered if false, molecules with bond count 0 but atom count 1 will return false (do not filter)
     * @return true if the given atom container should be discarded
     * @throws NullPointerException if parameter is 'null'
     */
    public static boolean shouldBeFiltered(IAtomContainer aMolecule, boolean areSingleAtomsFiltered) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered;
        try {
            if (areSingleAtomsFiltered) {
                tmpShouldBeFiltered = (ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(aMolecule)
                        || ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(aMolecule));
            } else {
                tmpShouldBeFiltered = (ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(aMolecule)
                        || aMolecule.getAtomCount() == 0);
            }

        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.WARNING,
                    anException.toString() + " Molecule ID: " + ErtlFunctionalGroupsFinderUtility.getIDForLogging(aMolecule),
                    anException);
            tmpShouldBeFiltered = true;
        }
        return tmpShouldBeFiltered;
    }

    /**
     * Checks whether the given molecule represented by an atom container needs to be preprocessed before it is passed
     * on to the ErtlFunctionalGroupsFinder.find() method because it is unconnected or contains charged atoms.
     * <br>It is advised to check via shouldBeFiltered() whether the given molecule should be discarded anyway before
     * calling this function.
     *
     * @param aMolecule the atom container to check
     * @return true is the given molecule needs to be preprocessed
     * @throws NullPointerException if parameter is 'null'
     */
    public static boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpNeedsPreprocessing;
        try {
            tmpNeedsPreprocessing = (ErtlFunctionalGroupsFinder.containsChargedAtom(aMolecule)
                    || ErtlFunctionalGroupsFinder.isStructureUnconnected(aMolecule));
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.WARNING,
                    anException.toString() + " Molecule ID: " + ErtlFunctionalGroupsFinderUtility.getIDForLogging(aMolecule),
                    anException);
            throw new NullPointerException("An unknown error occurred.");
        }
        return tmpNeedsPreprocessing;
    }

    /**
     * Checks whether the given molecule represented by an atom container can be passed on to the
     * ErtlFunctionalGroupsFinder.find() method without problems.
     * <br>This method will return false if the molecule contains any metal, metalloid, pseudo, or charged atoms, contains
     * multiple unconnected parts, or has an atom or bond count of zero.
     *
     * @param aMolecule the molecule to check
     * @return true if the given molecule is a valid parameter for ErtlFunctionalGroupsFinder.find() method
     * @throws NullPointerException if parameter is 'null'
     */
    public static boolean isValidArgumentForFindMethod(IAtomContainer aMolecule) throws NullPointerException {
        return ErtlFunctionalGroupsFinderUtility.isValidArgumentForFindMethod(aMolecule, true);
    }

    /**
     * Checks whether the given molecule represented by an atom container can be passed on to the
     * ErtlFunctionalGroupsFinder.find() method without problems.
     * <br>This method will return false if the molecule contains any metal, metalloid, pseudo, or charged atoms, contains
     * multiple unconnected parts, or has an atom or bond count of zero. If the second parameter is set to "false", single atom molecules
     * (bond count is 0) are accepted and not recommended to be filtered if they fulfill the other requirements.
     *
     * @param aMolecule the molecule to check
     * @param areSingleAtomsFiltered if false, molecules with bond count 0 but atom count 1 will return true (do not filter)
     * @return true if the given molecule is a valid parameter for ErtlFunctionalGroupsFinder.find() method
     * @throws NullPointerException if parameter is 'null'
     */
    public static boolean isValidArgumentForFindMethod(IAtomContainer aMolecule, boolean areSingleAtomsFiltered) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpIsValid;
        try {
            if (areSingleAtomsFiltered) {
                tmpIsValid = !(ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(aMolecule)
                        || ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(aMolecule)
                        || ErtlFunctionalGroupsFinder.containsChargedAtom(aMolecule)
                        || ErtlFunctionalGroupsFinder.isStructureUnconnected(aMolecule));
            } else {
                tmpIsValid = !(ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(aMolecule)
                        || aMolecule.getAtomCount() == 0
                        || ErtlFunctionalGroupsFinder.containsChargedAtom(aMolecule)
                        || ErtlFunctionalGroupsFinder.isStructureUnconnected(aMolecule));
            }
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE,
                    anException.toString() + " Molecule ID: " + ErtlFunctionalGroupsFinderUtility.getIDForLogging(aMolecule),
                    anException);
            tmpIsValid = false;
        }
        return tmpIsValid;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Preprocessing methods">
    /**
     * Returns the biggest unconnected component/structure of the given atom container, judging by the atom count. To
     * pre-check whether the atom container consists of multiple unconnected components, use isStructureUnconnected().
     * All set properties of aMolecule will be set as properties of the returned atom container.
     * <br>NOTE: The atom, bond etc. objects of the given atom container are re-used in the returned atom container but
     * the former remains unchanged
     * <br>Iterates through all unconnected components in the given atom container, so the method scales linearly with
     * O(n) with n: number of unconnected components.
     *
     * @param aMolecule the molecule whose biggest unconnected component should be found
     * @return the biggest (judging by the atom count) unconnected component of the given atom container
     * @throws NullPointerException if aMolecule is null or the biggest component
     */
    public static IAtomContainer selectBiggestUnconnectedComponent(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecules is 'null'.");
        IAtomContainerSet tmpUnconnectedComponentsSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestComponent = null;
        for (IAtomContainer tmpComponent : tmpUnconnectedComponentsSet.atomContainers()) {
            if (Objects.isNull(tmpBiggestComponent) || tmpBiggestComponent.getAtomCount() < tmpComponent.getAtomCount()) {
                tmpBiggestComponent = tmpComponent;
            }
        }
        Objects.requireNonNull(tmpBiggestComponent, "The resulting biggest component is 'null'.");
        tmpBiggestComponent.setProperties(aMolecule.getProperties());
        return tmpBiggestComponent;
    }

    /**
     * Neutralizes charged atoms in the given atom container by zeroing the formal atomic charges and filling up free
     * valences with implicit hydrogen atoms (according to the CDK atom types). This procedure allows a more general
     * charge treatment than a pre-defined transformation list but may produce "wrong" structures, e.g. it turns a
     * nitro NO2 group into a structure represented by the SMILES code "[H]O[N](=O)*" with an uncharged four-bonded
     * nitrogen atom (other examples are "*[N](*)(*)*", "[C]#[N]*" or "*S(*)(*)*"). Thus, an improved charge
     * neutralization scheme is desirable for future implementations.
     * <br>NOTE: This method changes major properties and the composition of the given IAtomContainer object! If you
     * want to retain your object unchanged for future calculations, use the IAtomContainer's
     * clone() method.
     * <br>Iterates through all atoms in the given atom container, so the method scales linearly with
     * O(n) with n: number of atoms.
     *
     * @param aMolecule the molecule to be neutralized
     * @throws NullPointerException if aMolecule is 'null' or one of its atoms
     * @throws CDKException if no matching atom type can be determined for one atom or there is a problem with adding
     * the implicit hydrogen atoms.
     */
    public static void neutralizeCharges(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Iterable<IAtom> tmpAtoms = aMolecule.atoms();
        for (IAtom tmpAtom : tmpAtoms) {
            ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpAtom, aMolecule);
        }
    }

    /**
     * Neutralizes a charged atom in the given parent atom container by zeroing the formal atomic charge and filling up free
     * valences with implicit hydrogen atoms (according to the CDK atom types).
     * <br>NOTE: This method changes major properties and the composition of the given IAtom and IAtomContainer object!
     * If you want to retain your objects unchanged for future calculations, use the IAtomContainer's
     * clone() method.
     *
     * @param anAtom the atom to be neutralized
     * @param aParentMolecule the molecule the atom belongs to
     * @throws NullPointerException if anAtom or aParentMolecule is 'null'
     * @throws CDKException if the atom is not part of the molecule or no matching atom type can be determined for the
     * atom or there is a problem with adding the implicit hydrogen atoms.
     */
    public static void neutralizeCharges(IAtom anAtom, IAtomContainer aParentMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Objects.requireNonNull(aParentMolecule, "Given parent molecule is 'null'.");
        boolean tmpIsAtomInMolecule = aParentMolecule.contains(anAtom);
        if (!tmpIsAtomInMolecule) {
            throw new CDKException("Given atom is not part of the given atom container.");
        }
        Integer tmpFormalChargeObject = anAtom.getFormalCharge();
        if (Objects.isNull(tmpFormalChargeObject)) {
            return;
        }
        int tmpFormalCharge = tmpFormalChargeObject.intValue();
        if (tmpFormalCharge != 0) {
            anAtom.setFormalCharge(0);
            IChemObjectBuilder tmpBuilder = aParentMolecule.getBuilder();
            if (Objects.isNull(tmpBuilder)) {
                throw new CDKException("Builder of the given atom container is 'null'.");
            }
            CDKHydrogenAdder tmpHAdder = CDKHydrogenAdder.getInstance(tmpBuilder);
            CDKAtomTypeMatcher tmpMatcher = CDKAtomTypeMatcher.getInstance(tmpBuilder);
            //Can throw CDKException
            IAtomType tmpMatchedType = tmpMatcher.findMatchingAtomType(aParentMolecule, anAtom);
            if (Objects.isNull(tmpMatchedType)) {
                throw new CDKException("Matched atom type is 'null'.");
            }
            AtomTypeManipulator.configure(anAtom, tmpMatchedType);
            //Can throw CDKException
            tmpHAdder.addImplicitHydrogens(aParentMolecule, anAtom);
        }
    }

    /**
     * Convenience method to perceive atom types for all IAtoms in the IAtomContainer, using the
     * CDK AtomContainerManipulator or rather the CDKAtomTypeMatcher. If the matcher finds a matching atom type, the
     * IAtom will be configured to have the same properties as the IAtomType. If no matching atom type is found, no
     * configuration is performed.
     * <br>Calling this method is equal to calling AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule).
     * It has been given its own method here because it is a necessary step in the preprocessing for
     * ErtlFunctionalGroupsFinder.
     * <br>NOTE: This method changes major properties of the given IAtomContainer object! If you
     * want to retain your object unchanged for future calculations, use the IAtomContainer's
     * clone() method.
     *
     * @param aMolecule the molecule to configure
     * @throws NullPointerException is aMolecule is 'null'
     * @throws CDKException when something went wrong with going through the AtomType options
     * @see AtomContainerManipulator#percieveAtomTypesAndConfigureAtoms(IAtomContainer)
     * @see CDKAtomTypeMatcher#findMatchingAtomType(IAtomContainer, IAtom)
     */
    public static void perceiveAtomTypesAndConfigureAtoms(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        //Might throw CDKException but it is unclear in what case
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
    }

    /**
     * Convenience method for applying the given aromaticity model to the given molecule. Any existing aromaticity flags
     * are removed - even if no aromatic bonds were found. This follows the idea of applying an aromaticity model to a
     * molecule such that the result is the same irrespective of existing aromatic flags.
     * <br>Calling this method is equal to calling Aromaticity.apply(aMolecule).
     * It has been given its own method here because it is a necessary step in the preprocessing for
     * ErtlFunctionalGroupsFinder.
     * <br>NOTE: This method changes major properties and the composition of the given IAtomContainer object! If you
     * want to retain your object unchanged for future calculations, use copy() in this class or the IAtomContainer's
     * clone() method.
     *
     * @param aMolecule the molecule to apply the model to
     * @param anAromaticityModel the model to apply; Note that the choice of electron donation model and cycle finder
     *                           algorithm has a heavy influence on the functional group detection of
     *                           ErtlFunctionalGroupsFinder
     * @return true if the molecule (or parts of it) is determined to be aromatic
     * @throws NullPointerException if a parameter is 'null'
     * @throws CDKException if a problem occurred with the cycle perception (see CDK docs)
     * @see Aromaticity#apply(IAtomContainer)
     */
    public static boolean applyAromaticityDetection(IAtomContainer aMolecule, Aromaticity anAromaticityModel) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Objects.requireNonNull(anAromaticityModel, "Given aromaticity model is 'null'.");
        boolean tmpIsAromatic = false;
        try {
            //throws CDKException if a problem occurred with the cycle perception (see CDK docs)
            //Note: Contrary to the docs, an Intractable exception might be thrown
            tmpIsAromatic = anAromaticityModel.apply(aMolecule);
        } catch (Intractable anIntractableException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE,
                    anIntractableException.toString() + " Molecule ID: " + ErtlFunctionalGroupsFinderUtility.getIDForLogging(aMolecule),
                    anIntractableException);
            String tmpMessage = anIntractableException.getMessage();
            Throwable tmpCause = anIntractableException.getCause();
            throw new CDKException(tmpMessage, tmpCause);
        }
        return tmpIsAromatic;
    }

    /**
     * Checks whether the given molecule represented by an atom container should be filtered instead of being passed
     * on to the ErtlFunctionalGroupsFinder.find() method and if not, applies necessary preprocessing steps.
     * In the second case, this method applies preprocessing
     * to the given atom container that is always needed (setting atom types and applying an aromaticity model) and
     * preprocessing steps that are only needed in specific cases (selecting the biggest unconnected component, neutralizing
     * charges). Molecules processed by this method can be passed on to find() without problems (Caution: The return value
     * of this method is 'null' if the molecule should be filtered!).
     * <br>NOTE: This method changes major properties and the composition of the given IAtomContainer object! If you
     * want to retain your object unchanged for future calculations, use the IAtomContainer's
     * clone() method.
     * <br>NOTE2: The returned IAtomContainer object is the same as the one given as parameter!
     *
     * @param aMolecule the molecule to check and process
     * @param anAromaticityModel the aromaticity model to apply to the molecule in preprocessing; Note: The chosen
     * ElectronDonation model can massively influence the extracted function groups of a molecule when using
     * ErtlFunctionGroupsFinder!
     * @return the preprocessed atom container or 'null' if the molecule should be discarded
     * @throws NullPointerException if a parameter is 'null'; Note: All other exceptions are caught and logged by this
     * class' logger
     */
    public static IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule, Aromaticity anAromaticityModel) throws NullPointerException {
        return ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(aMolecule, anAromaticityModel, true);
    }

    /**
     * Checks whether the given molecule represented by an atom container should be filtered instead of being passed
     * on to the ErtlFunctionalGroupsFinder.find() method and if not, applies necessary preprocessing steps.
     * In the second case, this method applies preprocessing
     * to the given atom container that is always needed (setting atom types and applying an aromaticity model) and
     * preprocessing steps that are only needed in specific cases (selecting the biggest unconnected component, neutralizing
     * charges). Molecules processed by this method can be passed on to find() without problems (Caution: The return value
     * of this method is 'null' if the molecule should be filtered!).
     * <br>NOTE: This method changes major properties and the composition of the given IAtomContainer object! If you
     * want to retain your object unchanged for future calculations, use the IAtomContainer's
     * clone() method.
     * <br>NOTE2: The returned IAtomContainer object is the same as the one given as parameter!
     *
     * @param aMolecule the molecule to check and process
     * @param anAromaticityModel the aromaticity model to apply to the molecule in preprocessing; Note: The chosen
     * ElectronDonation model can massively influence the extracted functional groups of a molecule when using
     * ErtlFunctionGroupsFinder!
     * @param areSingleAtomsFiltered if false, molecules with bond count 0 but atom count 1 will be processed and
     *                               not return null
     * @return the preprocessed atom container or 'null' if the molecule should be discarded
     * @throws NullPointerException if a parameter is 'null'; Note: All other exceptions are caught and logged by this
     * class' logger
     */
    public static IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule, Aromaticity anAromaticityModel, boolean areSingleAtomsFiltered) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'.");
        Objects.requireNonNull(anAromaticityModel, "Given aromaticity model is 'null'.");
        try {
            ErtlFunctionalGroupsFinderUtility.perceiveAtomTypesAndConfigureAtoms(aMolecule);
            //Filter
            if (areSingleAtomsFiltered) {
                boolean tmpIsAtomOrBondCountZero = ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(aMolecule);
                if (tmpIsAtomOrBondCountZero) {
                    return null;
                }
            } else {
                boolean tmpIsAtomCountZero = aMolecule.getAtomCount() == 0;
                if (tmpIsAtomCountZero) {
                    return null;
                }
            }
            //From structures containing two or more unconnected structures (e.g. ions) choose the largest structure
            boolean tmpIsUnconnected = ErtlFunctionalGroupsFinder.isStructureUnconnected(aMolecule);
            if (tmpIsUnconnected) {
                aMolecule = ErtlFunctionalGroupsFinderUtility.selectBiggestUnconnectedComponent(aMolecule);
            }
            //Filter
            boolean tmpContainsInvalidAtoms = ErtlFunctionalGroupsFinder.containsMetalMetalloidOrPseudoAtom(aMolecule);
            if (tmpContainsInvalidAtoms) {
                return null;
            }
            //Neutralize charges if there are any
            boolean tmpIsCharged = ErtlFunctionalGroupsFinder.containsChargedAtom(aMolecule);
            if (tmpIsCharged) {
                ErtlFunctionalGroupsFinderUtility.neutralizeCharges(aMolecule);
            }
            //Application of aromaticity model
            ErtlFunctionalGroupsFinderUtility.applyAromaticityDetection(aMolecule, anAromaticityModel);
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE,
                    anException.toString() + " Molecule ID: " + ErtlFunctionalGroupsFinderUtility.getIDForLogging(aMolecule),
                    anException);
            return null;
        }
        return aMolecule;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Additional functionalities">

    /**
     * Replaces the environmental carbon or pseudo-atoms (new IAtom objects) inserted by the EFGF in an identified
     * functional group with the carbon IAtom objects from the original molecule object.
     * <br>Important note: This method only works if the atom container has not been cloned for the extraction of
     * functional groups by ErtlFunctionalGroupsFinder. Use the method
     * "List{@literal <}IAtomContainer{@literal >} find(IAtomContainer container, boolean clone)" with clone set to false for this purpose.
     * <br>Also note that the result differs if the environment has been generalized by the EFGF or not. In the former
     * case, only environmental carbon atoms replaced by R-atoms in the generalized FG are restored.
     *
     * @param aListOfFunctionalGroups functional groups of the molecule identified by EFGF
     * @param aMolecule original structure in which the groups were identified
     * @param aConvertExplicitHydrogens should explicit hydrogen atoms in the functional groups be converted to implicit
     *                                  hydrogens
     * @param aFillEmptyValences should empty valences on the restored environmental carbon atoms be filled with
     *                           implicit hydrogen atoms
     * @param aBuilder a chem object builder instance
     * @throws NullPointerException if a parameter is null
     * @throws IllegalArgumentException if one of the functional groups does not originate from the given molecule
     *                                  or the molecule has been cloned for the extraction of functional groups
     * @author Michael Wenk, Jonas Schaub
     */
    public static void restoreOriginalEnvironmentalCarbons(
            List<IAtomContainer> aListOfFunctionalGroups,
            IAtomContainer aMolecule,
            boolean aConvertExplicitHydrogens,
            boolean aFillEmptyValences,
            IChemObjectBuilder aBuilder)
            throws NullPointerException, IllegalArgumentException {
        //<editor-fold desc="Parameter checks">
        Objects.requireNonNull(aListOfFunctionalGroups, "Given list of functional groups is null.");
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        Objects.requireNonNull(aBuilder, "Given chem object builder is null.");
        if (aListOfFunctionalGroups.isEmpty()) {
            return;
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Given molecule is empty.");
        }
        for (IAtomContainer tmpFG : aListOfFunctionalGroups) {
            boolean tmpIsFGofMolecule = false;
            for (IAtom tmpAtom : tmpFG.atoms()) {
                if (aMolecule.contains(tmpAtom)) {
                    tmpIsFGofMolecule = true;
                }
            }
            if (!tmpIsFGofMolecule) {
                throw new IllegalArgumentException("At least one functional group has been given that does not originate " +
                        "from the given molecule or the molecule has been cloned for the extraction of functional groups.");
            }
        }
        //</editor-fold>
        CDKHydrogenAdder tmpHadder = CDKHydrogenAdder.getInstance(aBuilder);
        for (int i = 0; i < aListOfFunctionalGroups.size(); i++) {
            IAtomContainer tmpFG = aListOfFunctionalGroups.get(i);
            //convert explicit hydrogens to implicit
            if (aConvertExplicitHydrogens) {
                AtomContainerManipulator.suppressHydrogens(tmpFG);
            }
            //create a list of all atoms of the group because of atom removals and additions in group atom container
            List<IAtom> tmpListofFGatoms = new ArrayList<>();
            for (IAtom tmpAtom : tmpFG.atoms()) {
                tmpListofFGatoms.add(tmpAtom);
            }
            for (IAtom tmpAtom : tmpListofFGatoms) {
                //technically, all elements except carbon should be excluded but this way, it is the easiest to also include
                // pseudo atoms
                if (tmpAtom.getAtomicNumber().equals(1)) {
                    continue;
                }
                //detect whether the current atom is an "unknown" one, inserted as new environmental IAtom object
                if (!aMolecule.contains(tmpAtom)) {
                    //environmental carbon and pseudo-atoms (carbon or hydrogen) added by the EFGF can only have one bond partner in the FG.
                    // identify its bond partner in the FG that should be part of the original molecule
                    IAtom tmpConnectedAtomInGroup = tmpFG.getConnectedAtomsList(tmpAtom).get(0);
                    //remove the inserted atom and the bond to it
                    tmpFG.removeBond(tmpAtom, tmpConnectedAtomInGroup);
                    tmpFG.removeAtom(tmpAtom);
                    //starting from the parent atom search for neighboring carbons which are not already in the group and add them
                    for (IAtom tmpConnectedAtomInOriginalStructure : aMolecule.getConnectedAtomsList(tmpConnectedAtomInGroup)) {
                        if (tmpConnectedAtomInOriginalStructure.getSymbol().equals("C")
                                && !tmpFG.contains(tmpConnectedAtomInOriginalStructure)) {
                            tmpFG.addAtom(tmpConnectedAtomInOriginalStructure);
                            tmpFG.addBond(aMolecule.getBond(tmpConnectedAtomInGroup, tmpConnectedAtomInOriginalStructure));
                        }
                    }
                }
            }
            if (aFillEmptyValences) {
                try {
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpFG);
                    tmpHadder.addImplicitHydrogens(tmpFG);
                } catch (CDKException aCDKException) {
                    ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.WARNING, aCDKException.toString(), aCDKException);
                    continue;
                }
            }
        }
    }

    /**
     * Gives the pseudo SMILES code for a given molecule / functional group. In this notation, aromatic atoms are marked
     * by asterisks (*) and pseudo atoms are indicated by 'R'.
     * <br>The function generates the SMILES string of the given molecule using CDK's SmilesGenerator and then
     * replaces lowercase c, n, o etc. by C*, N*, O* etc. and wildcards ('*') by 'R' in the resulting string.
     * For that, the function iterates through all characters in the generated SMILES string.
     * <br>Note: All pseudo atoms or atoms that are represented by a wildcard ('*') in the generated SMILES string
     * (e.g. the element [Uup] is interpreted by the CDK SmilesGenerator as a wildcard) are turned into an 'R' atom.
     *
     * @param aMolecule the molecule whose pseudo SMILES code to generate
     * @return the pseudo SMILES representation as a string
     * @throws NullPointerException if aMolecule is 'null'
     * @throws CDKException if the SMILES code of aMolecule cannot be generated
     */
    public static String createPseudoSmilesCode(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        String tmpPseudoSmilesCode;
        try {
            //Might throw CDKException if the SMILES string cannot be created or NullPointerException if an atom has an
            //  undefined number of implicit hydrogen atoms in the SMILES string
            tmpPseudoSmilesCode = tmpSmilesGenerator.create(aMolecule);
        } catch (NullPointerException anException) {
            throw new CDKException(anException.getMessage(), anException);
        }
        tmpPseudoSmilesCode = tmpPseudoSmilesCode.replaceAll("\\*", "R");
        tmpPseudoSmilesCode = tmpPseudoSmilesCode.replaceAll("\\[se", "[Se*");
        StringBuilder tmpStringBuilder = new StringBuilder(tmpPseudoSmilesCode);
        int tmpLength = tmpStringBuilder.length();
        for (int tmpIndex = 0; tmpIndex < tmpLength; tmpIndex++) {
            char tmpChar = tmpStringBuilder.charAt(tmpIndex);
            char tmpPrevChar = '_';
            char tmpPrevPrevChar = '_';
            if (tmpIndex > 0) {
                tmpPrevChar = tmpStringBuilder.charAt(tmpIndex - 1);
            }
            if (tmpIndex > 1) {
                tmpPrevPrevChar = tmpStringBuilder.charAt(tmpIndex - 2);
            }
            switch (tmpChar) {
                case 'c':
                    //c in [Sc], [Tc], and [Ac] should not be replaced
                    if ((tmpPrevChar == 'S' || tmpPrevChar == 'T' || tmpPrevChar == 'A') && tmpPrevPrevChar == '[') {
                        break;
                    } else {
                        tmpStringBuilder.setCharAt(tmpIndex, 'C');
                        tmpStringBuilder.insert(tmpIndex + 1, '*');
                    }
                    break;
                case 'n':
                    //n in [Mn], [Zn], [Cn], [In], [Sn], and [Rn] should not be replaced
                    if ((tmpPrevChar == 'M' || tmpPrevChar == 'Z' || tmpPrevChar == 'C' || tmpPrevChar == 'I' || tmpPrevChar == 'S' || tmpPrevChar == 'R') && tmpPrevPrevChar == '[') {
                        break;
                    } else {
                        tmpStringBuilder.setCharAt(tmpIndex, 'N');
                        tmpStringBuilder.insert(tmpIndex + 1, '*');
                    }
                    break;
                case 's':
                    //s in [Cs], [Os], [As], [Es], [Hs], and [Uus] should not be replaced
                    if ((tmpPrevChar == 'C' || tmpPrevChar == 'O' || tmpPrevChar == 'A' || tmpPrevChar == 'E' || tmpPrevChar == 'H') && tmpPrevPrevChar == '[') {
                        break;
                    } else if (tmpPrevChar == 'u' && tmpPrevPrevChar == 'U') {
                        break;
                    } else {
                        tmpStringBuilder.setCharAt(tmpIndex, 'S');
                        tmpStringBuilder.insert(tmpIndex + 1, '*');
                    }
                    break;
                case 'o':
                    //o in [Mo], [Co], [Po], [Uuo], [Ho], and [No] should not be replaced
                    if ((tmpPrevChar == 'M' || tmpPrevChar == 'C' || tmpPrevChar == 'P' || tmpPrevChar == 'H' || tmpPrevChar == 'N') && tmpPrevPrevChar == '[') {
                        break;
                    } else if (tmpPrevChar == 'u' && tmpPrevPrevChar == 'U') {
                        break;
                    } else {
                        tmpStringBuilder.setCharAt(tmpIndex, 'O');
                        tmpStringBuilder.insert(tmpIndex + 1, '*');
                    }
                    break;
                case 'p':
                    //p in [Uup] and [Np] should not be replaced
                    if (tmpPrevChar == 'N' && tmpPrevPrevChar == '[') {
                        break;
                    } else if (tmpPrevChar == 'u' && tmpPrevPrevChar == 'U') {
                        break;
                    } else {
                        tmpStringBuilder.setCharAt(tmpIndex, 'P');
                        tmpStringBuilder.insert(tmpIndex + 1, '*');
                    }
                    break;
                default:
                    break;
            }
            tmpLength = tmpStringBuilder.length();
        }
        tmpPseudoSmilesCode = tmpStringBuilder.toString();
        return tmpPseudoSmilesCode;
    }
    //</editor-fold>
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Returns the CDK title or ID of the given molecule.
     *
     * @param aMolecule the molecule to determine the title or ID of
     * @return the CDK title, title or ID of the given molecule (depending on which property is set)
     * @throws NullPointerException if aMolecule is 'null'
     */
    private static String getIDForLogging(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        String tmpCdkTitle = aMolecule.getProperty(CDKConstants.TITLE);
        String tmpTitle = aMolecule.getTitle();
        String tmpID = aMolecule.getID();
        if (!Objects.isNull(tmpCdkTitle) && !tmpCdkTitle.isEmpty()) {
            return "CDK title: " + tmpCdkTitle;
        } else if (!Objects.isNull(tmpTitle) && !tmpTitle.isEmpty()) {
            return "Title: " + tmpTitle;
        } else if (!Objects.isNull(tmpID) && !tmpID.isEmpty()) {
            return "ID: " + tmpID;
        } else {
            return "No title or id could be determined.";
        }
    }
    //</editor-fold>
}
