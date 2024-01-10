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

import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.ISingleElectron;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Queue;
import java.util.Set;

/**
 * Finds and extracts a molecule's functional groups in a purely rule-based manner.
 *
 * This class implements Peter Ertl's algorithm for the automated detection and extraction
 * of functional groups in organic molecules
 * (<a href="https://doi.org/10.1186/s13321-017-0225-z">[Ertl P. An algorithm to identify functional groups in organic molecules. J Cheminform. 2017; 9:36.])</a>.
 * <p>
 *     Note: this implementation is not thread-safe. Each parallel thread should have its own instance of this class.
 * </p>
 *
 * @author Sebastian Fritsch, Jonas Schaub
 * @version 1.2.1
 */
public class ErtlFunctionalGroupsFinder {

    /**
     * Defines the mode for generalizing functional group environments (default) or keeping them whole.
     */
    public static enum Mode {
        /**
         * Default mode including the generalization step.
         */
        DEFAULT,
        /**
         * Skips the generalization step. Functional groups will keep their full environment.
         */
        NO_GENERALIZATION;
    }
    //
    /**
     * Defines whether an environmental carbon atom is aromatic or aliphatic. Only for internal use for caching this
     * info in the EnvironmentalC instances (see private class below).
     */
    private static enum EnvironmentalCType {
        /**
         * Aromatic environmental carbon.
         */
        C_AROMATIC,
        /**
         * Aliphatic environmental carbon.
         */
        C_ALIPHATIC;
    }
    //
    /**
     * Describes one carbon atom in the environment of a marked atom. It can either be aromatic
     * or aliphatic and also contains a clone of its connecting bond.
     */
    private class EnvironmentalC {

        /**
         * Indicates whether carbon atom is aromatic or aliphatic.
         */
        private EnvironmentalCType type;
        //
        /**
         * Bond index of the original C atom.
         */
        private int bondIndex;
        //
        /**
         * Order of the bond connecting this environmental C atom to the marked functional group atom.
         */
        private IBond.Order bondOrder;
        //
        /**
         * Stereo information of the bond connecting this environmental C atom to the marked functional group atom.
         */
        private IBond.Stereo bondStereo;
        //
        /**
         * Flags of the bond connecting this environmental C atom to the marked functional group atom. IChemObjecflags
         * are properties defined by an integer value (array position) and a boolean value.
         */
        private boolean[] bondFlags;
        //
        /**
         * Default constructor defining all fields. Order, stereo, and flags are taken from the IBond object directly.
         *
         * @param aType aromatic or aliphatic
         * @param aConnectingBond bond instance connecting to the marked atom
         * @param anIndexInBond index of the atom in the connecting bond
         */
        public EnvironmentalC(EnvironmentalCType aType, IBond aConnectingBond, int anIndexInBond) {
            this.type = aType;
            this.bondIndex = anIndexInBond;
            this.bondOrder = aConnectingBond.getOrder();
            this.bondStereo = aConnectingBond.getStereo();
            this.bondFlags = aConnectingBond.getFlags();
        }
        //
        /**
         * Returns the type, i.e. whether this carbon atom is aromatic or aliphatic.
         *
         * @return EnvironmentalCType enum constant
         */
        public EnvironmentalCType getType() {
            return this.type;
        }
        //
        /**
         * Method for translating this instance back into a "real" IAtom instance when expanding the functional group
         * environment, transferring all the cached properties, except the type(!).
         *
         * @param aTargetAtom marked functional group atom
         * @param anEnvCAtom new carbon atom instance that should receive all the cached properties except the type(!);
         *                   element, atom type "C" and implicit hydrogen count = 0 should be set already; type can later
         *                   be set via .setIsAromatic(boolean);
         * @return new bond connecting marked FG atom and environment atom in the correct order and with the cached properties
         */
        public IBond createBond(IAtom aTargetAtom, IAtom anEnvCAtom) {
            IBond tmpBond = aTargetAtom.getBuilder().newInstance(IBond.class);
            if (this.bondIndex == 0) {
                tmpBond.setAtoms(new IAtom[] {anEnvCAtom, aTargetAtom});
            }
            else {
                tmpBond.setAtoms(new IAtom[] {aTargetAtom, anEnvCAtom});
            }
            tmpBond.setOrder(this.bondOrder);
            tmpBond.setStereo(this.bondStereo);
            tmpBond.setFlags(this.bondFlags);
            return tmpBond;
        }
    }
    //
    /**
     * CDK logging tool instance for this class.
     */
    public static final ILoggingTool LOGGING_TOOL = LoggingToolFactory.createLoggingTool(ErtlFunctionalGroupsFinder.class);
    //
    /**
     * Property name for marking carbonyl carbon atoms via IAtom properties.
     */
    public static final String CARBONYL_C_MARKER = "EFGF-Carbonyl-C";
    //
    /**
     * Set of atomic numbers that are accepted in the input molecule if the strict input restrictions are activated
     * (excludes metal and metalloid elements, only organic elements included).
     */
    public static final Set<Integer> NONMETAL_ATOMIC_NUMBERS = Set.of(1, 2, 6, 7, 8, 9, 10, 15, 16, 17, 18, 34, 35, 36, 53, 54, 86);
    //
    /**
     * Environment mode setting, defining whether environments should be generalized (default) or kept as whole.
     */
    private Mode envMode;
    //
    /**
     * Map of bonds in the input molecule, cache(!).
     */
    private EdgeToBondMap bondMap;
    //
    /**
     * Adjacency list representation of input molecule, cache(!).
     */
    private int[][] adjList;
    //
    /**
     * Set for atoms marked as being part of a functional group, represented by an internal index based on the atom
     * count in the input molecule, cache(!).
     */
    private HashSet<Integer> markedAtoms;
    //
    /**
     * HashMap for storing aromatic hetero-atom indices and whether they have already been assigned to a larger functional
     * group. If false, they form single-atom FG by themselves, cache(!).
     *
     * key: atom idx, value: isInGroup
     */
    private HashMap<Integer, Boolean> aromaticHeteroAtomIndicesToIsInGroupBoolMap;
    //
    /**
     * HashMap for storing marked atom to connected environmental carbon atom relations, cache(!).
     */
    private HashMap<IAtom, List<EnvironmentalC>> markedAtomToConnectedEnvCMap;
    //
    /**
     * Default constructor for ErtlFunctionalGroupsFinder with functional group generalization turned ON.
     */
    public ErtlFunctionalGroupsFinder() {
        this(Mode.DEFAULT);
    }
    //
    /**
     * Constructor for ErtlFunctionalGroupsFinder that allows setting the treatment of environments in the identified
     * functional groups. Default: environments will be generalized; no generalization: environments will be kept as whole.
     *
     * @param anEnvMode mode for treating functional group environments (see {@link ErtlFunctionalGroupsFinder.Mode}).
     */
    public ErtlFunctionalGroupsFinder(Mode anEnvMode) {
        Objects.requireNonNull(anEnvMode, "Given environment mode cannot be null.");
        this.envMode = anEnvMode;
    }
    //
    /**
     * Allows setting the treatment of functional group environments after extraction. Default: environments will be
     * generalized; no generalization: environments will be kept as whole.
     *
     * @param anEnvMode mode for treating functional group environments (see {@link ErtlFunctionalGroupsFinder.Mode}).
     */
    public void setEnvMode(Mode anEnvMode) {
        Objects.requireNonNull(anEnvMode, "Given environment mode cannot be null.");
        this.envMode = anEnvMode;
    }
    //
    /**
     * Returns the current setting for the treatment of functional group environments after extraction.
     *
     * @return currently set environment mode
     */
    public Mode getEnvMode() {
        return this.envMode;
    }
    //
    /**
     * Find all functional groups in a molecule. The input atom container instance is cloned before processing to leave
     * the input container intact.
     * <p>
     *     Note: The strict input restrictions from previous versions (no charged atoms, metals, metalloids or
     *     unconnected components) do not apply anymore by default. They can be turned on again in another variant of
     *     this method below.
     * </p>
     *
     * @param aMolecule the molecule to identify functional groups in
     * @throws CloneNotSupportedException if cloning is not possible
     * @return a list with all functional groups found in the molecule
     */
    public List<IAtomContainer> find(IAtomContainer aMolecule) throws CloneNotSupportedException {
        return this.find(aMolecule, true, false);
    }
    //
    /**
     * Find all functional groups in a molecule.
     * <p>
     *     Note: The strict input restrictions from previous versions (no charged atoms, metals, metalloids or
     *     unconnected components) do not apply anymore by default. They can be turned on again in another variant of
     *     this method below.
     * </p>
     *
     * @param aMolecule the molecule to identify functional groups in
     * @param aShouldInputBeCloned use 'false' to reuse the input container's bonds and atoms in the extraction of the functional
     *                             groups; this may speed up the extraction and lower the memory consumption for processing large
     *                             amounts of data but corrupts the original input container; use 'true' to work with a clone and
     *                             leave the input container intact
     * @throws CloneNotSupportedException if cloning is not possible
     * @return a list with all functional groups found in the molecule
     */
    public List<IAtomContainer> find (IAtomContainer aMolecule, boolean aShouldInputBeCloned) throws CloneNotSupportedException {
        return this.find(aMolecule, aShouldInputBeCloned, false);
    }

    /**
     * Find all functional groups in a molecule.
     *
     * @param aMolecule the molecule to identify functional groups in
     * @param aShouldInputBeCloned use 'false' to reuse the input container's bonds and atoms in the extraction of the functional
     *                             groups; this may speed up the extraction and lower the memory consumption for processing large
     *                             amounts of data but corrupts the original input container; use 'true' to work with a clone and
     *                             leave the input container intact
     * @param anAreInputRestrictionsApplied if true, the input must consist of one connected structure and may not
     *                                      contain charged atoms, metals or metalloids; an IllegalArgumentException will
     *                                      be thrown otherwise
     * @throws CloneNotSupportedException if cloning is not possible
     * @throws IllegalArgumentException if input restrictions are applied and the given molecule does not fulfill them
     * @return a list with all functional groups found in the molecule
     */
    public List<IAtomContainer> find(IAtomContainer aMolecule, boolean aShouldInputBeCloned, boolean anAreInputRestrictionsApplied)
            throws CloneNotSupportedException, IllegalArgumentException {
        IAtomContainer tmpMolecule;
        if (aShouldInputBeCloned) {
            tmpMolecule = aMolecule.clone();
        } else {
            tmpMolecule = aMolecule;
        }
        if (anAreInputRestrictionsApplied) {
            this.checkConstraints(tmpMolecule);
        }
        for (IAtom tmpAtom : tmpMolecule.atoms()) {
            if(Objects.isNull(tmpAtom.getImplicitHydrogenCount())) {
                tmpAtom.setImplicitHydrogenCount(0);
            }
        }
        this.bondMap = EdgeToBondMap.withSpaceFor(tmpMolecule);
        this.adjList = GraphUtil.toAdjList(tmpMolecule, this.bondMap);
        this.markAtoms(tmpMolecule);
        // extract raw groups
        List<IAtomContainer> tmpFunctionalGroupsList = this.extractGroups(tmpMolecule);
        // handle environment
        if (this.envMode == Mode.DEFAULT) {
            this.expandGeneralizedEnvironments(tmpFunctionalGroupsList);
        } else if (this.envMode == Mode.NO_GENERALIZATION) {
            this.expandFullEnvironments(tmpFunctionalGroupsList);
        } else {
            throw new IllegalArgumentException("Unknown mode.");
        }
        this.clearCache();
        return tmpFunctionalGroupsList;
    }

    /**
     * Clear caches related to the input molecule.
     */
    private void clearCache() {
        this.bondMap = null;
        this.adjList = null;
        this.markedAtoms = null;
        this.aromaticHeteroAtomIndicesToIsInGroupBoolMap = null;
        this.markedAtomToConnectedEnvCMap = null;
    }

    /**
     * Mark all atoms and store them in a set for further processing.
     *
     * @param aMolecule molecule with atoms to mark
     */
    private void markAtoms(IAtomContainer aMolecule) {
        if (this.isDbg()) {
            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug("########## Starting search for atoms to mark ... ##########");
        }
        // store marked atoms
        this.markedAtoms = new HashSet<>((int) ((aMolecule.getAtomCount() / 0.75f) + 2), 0.75f);
        // store aromatic heteroatoms
        this.aromaticHeteroAtomIndicesToIsInGroupBoolMap = new HashMap<>((int) ((aMolecule.getAtomCount() / 0.75f) + 2), 0.75f);
        //TODO set and use a more explicit and trustworthy index?
        for (int idx = 0; idx < aMolecule.getAtomCount(); idx++) {
            // skip atoms that were already marked in a previous iteration
            if (this.markedAtoms.contains(idx)) {
                continue;
            }
            IAtom tmpAtom = aMolecule.getAtom(idx);
            // skip aromatic atoms but add aromatic HETERO-atoms to map for later processing
            if (tmpAtom.isAromatic()) {
                if (ErtlFunctionalGroupsFinder.isHeteroatom(tmpAtom)) {
                    this.aromaticHeteroAtomIndicesToIsInGroupBoolMap.put(idx, false);
                }
                continue;
            }
            int tmpAtomicNr = tmpAtom.getAtomicNumber();
            // if C...
            if (tmpAtomicNr == 6) {
                // to detect if for loop ran with or without marking the C atom
                boolean tmpIsMarked = false;
                // count for the number of connected O, N & S atoms to detect acetal carbons
                int tmpConnectedONSatomsCounter = 0;
                for (int tmpConnectedIdx : this.adjList[idx]) {
                    IAtom tmpConnectedAtom = aMolecule.getAtom(tmpConnectedIdx);
                    IBond tmpConnectedBond = this.bondMap.get(idx, tmpConnectedIdx);

                    // if connected to heteroatom or C in aliphatic double or triple bond... [CONDITIONS 2.1 & 2.2]
                    if (tmpConnectedAtom.getAtomicNumber() != 1
                            && ((tmpConnectedBond.getOrder() == Order.DOUBLE || tmpConnectedBond.getOrder() == Order.TRIPLE)
                            && !tmpConnectedBond.isAromatic())) {

                        // set the *connected* atom as marked (add() true if this set did not already contain the specified element)
                        if (this.markedAtoms.add(tmpConnectedIdx)) {
                            if (ErtlFunctionalGroupsFinder.isDbg()) {
                                ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                        "Marking Atom #%d (%s) - Met condition %s",
                                        tmpConnectedIdx,
                                        tmpConnectedAtom.getSymbol(),
                                        tmpConnectedAtom.getAtomicNumber() == 6 ? "2.1/2.2" : "1"));
                            }
                        }
                        // set the *current* atom as marked and break out of connected atoms
                        tmpIsMarked = true;
                        if (ErtlFunctionalGroupsFinder.isDbg()) {
                            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                    "Marking Atom #%d (%s) - Met condition 2.1/2.2",
                                    idx,
                                    tmpAtom.getSymbol()));
                        }
                        // but check for carbonyl-C before break
                        if (tmpConnectedAtom.getAtomicNumber() == 8
                                && tmpConnectedBond.getOrder() == Order.DOUBLE
                                && this.adjList[idx].length == 3) {
                            tmpAtom.setProperty(CARBONYL_C_MARKER, true);
                            if (ErtlFunctionalGroupsFinder.isDbg())  {
                                ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug("- was flagged as Carbonly-C");
                            }
                        }
                        // break out of connected atoms
                        break;
                    } else if ((tmpConnectedAtom.getAtomicNumber() == 7
                            || tmpConnectedAtom.getAtomicNumber() == 8
                            || tmpConnectedAtom.getAtomicNumber() == 16)
                            && tmpConnectedBond.getOrder() == Order.SINGLE) {
                        // if connected to O/N/S in single bond...
                        // if connected O/N/S is not aromatic...
                        if (!tmpConnectedAtom.isAromatic()) {
                            // set the connected O/N/S atom as marked
                            this.markedAtoms.add(tmpConnectedIdx);
                            if (ErtlFunctionalGroupsFinder.isDbg()) {
                                ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                        "Marking Atom #%d (%s) - Met condition 1",
                                        tmpConnectedIdx,
                                        tmpConnectedAtom.getSymbol()));
                            }
                            // if "acetal C" (2+ O/N/S in single bonds connected to sp3-C)... [CONDITION 2.3]
                            boolean tmpIsAllSingleBonds = true;
                            for (int tmpConnectedInSphere2Idx : this.adjList[tmpConnectedIdx]) {
                                IBond tmpSphere2Bond = this.bondMap.get(tmpConnectedIdx, tmpConnectedInSphere2Idx);
                                if (tmpSphere2Bond.getOrder() != Order.SINGLE) {
                                    tmpIsAllSingleBonds = false;
                                    break;
                                }
                            }
                            if (tmpIsAllSingleBonds) {
                                tmpConnectedONSatomsCounter++;
                                if (tmpConnectedONSatomsCounter > 1 && this.adjList[idx].length + tmpAtom.getImplicitHydrogenCount() == 4) {
                                    // set as marked and break out of connected atoms
                                    tmpIsMarked = true;
                                    if (ErtlFunctionalGroupsFinder.isDbg()) {
                                        ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                                "Marking Atom #%d (%s) - Met condition 2.3",
                                                idx,
                                                tmpAtom.getSymbol()));
                                    }
                                    break;
                                }
                            }
                        }
                        // if part of oxirane, aziridine, or thiirane ring... [CONDITION 2.4]
                        for (int tmpConnectedInSphere2Idx : this.adjList[tmpConnectedIdx]) {
                            IAtom tmpConnectedInSphere2Atom = aMolecule.getAtom(tmpConnectedInSphere2Idx);
                            if (tmpConnectedInSphere2Atom.getAtomicNumber() == 6) {
                                for (int tmpConnectedInSphere3Idx : this.adjList[tmpConnectedInSphere2Idx]) {
                                    IAtom tmpConnectedInSphere3Atom = aMolecule.getAtom(tmpConnectedInSphere3Idx);
                                    if (tmpConnectedInSphere3Atom.equals(tmpAtom)) {
                                        // set connected atoms as marked
                                        this.markedAtoms.add(tmpConnectedInSphere2Idx);
                                        this.markedAtoms.add(tmpConnectedInSphere3Idx);
                                        if (ErtlFunctionalGroupsFinder.isDbg()) {
                                            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                                    "Marking Atom #%d (%s) - Met condition 2.4",
                                                    tmpConnectedInSphere2Idx,
                                                    tmpConnectedInSphere2Atom.getSymbol()));
                                            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                                    "Marking Atom #%d (%s) - Met condition 2.4",
                                                    tmpConnectedInSphere3Idx,
                                                    tmpConnectedInSphere3Atom.getSymbol()));
                                        }
                                        // set current atom as marked and break out of connected atoms
                                        tmpIsMarked = true;
                                        if (ErtlFunctionalGroupsFinder.isDbg()) {
                                            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                                                    "Marking Atom #%d (%s) - Met condition 2.4",
                                                    idx,
                                                    tmpAtom.getSymbol()));
                                        }
                                        break;
                                    }
                                }
                            }
                        } //end of for loop iterating over second sphere atoms
                    } // end of else if connected to O/N/S in single bond
                } //end of for loop that iterates over all connected atoms of the carbon atom
                if (tmpIsMarked) {
                    this.markedAtoms.add(idx);
                    continue;
                }
                // if none of the conditions 2.X apply, we have an unmarked C (not relevant here)
            } else if (tmpAtomicNr == 1){
                // if H...
                // convert to implicit H
                IAtom tmpConnectedAtom;
                try {
                    tmpConnectedAtom = aMolecule.getAtom(this.adjList[idx][0]);
                } catch(ArrayIndexOutOfBoundsException anException) {
                    //TODO: do sth here?
                    break;
                }
                if (Objects.isNull(tmpConnectedAtom.getImplicitHydrogenCount())) {
                    tmpConnectedAtom.setImplicitHydrogenCount(1);
                } else {
                    tmpConnectedAtom.setImplicitHydrogenCount(tmpConnectedAtom.getImplicitHydrogenCount() + 1);
                }
                continue;
            } else {
                // if heteroatom... (CONDITION 1)
                this.markedAtoms.add(idx);
                if (ErtlFunctionalGroupsFinder.isDbg()) {
                    ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                            "Marking Atom #%d (%s) - Met condition 1",
                            idx,
                            tmpAtom.getSymbol()));
                }
                continue;
            }
        } //end of for loop that iterates over all atoms in the mol
        if (ErtlFunctionalGroupsFinder.isDbg()) {
            ErtlFunctionalGroupsFinder.LOGGING_TOOL.debug(String.format(
                    "########## End of search. Marked %d/%d atoms. ##########",
                    this.markedAtoms.size(),
                    aMolecule.getAtomCount()));
        }
    }
    //
    /**
     * Searches the molecule for groups of connected marked atoms and extracts each as a new functional group.
     * The extraction process includes marked atom's "environments". Connected H's are captured implicitly.
     *
     * @param molecule the molecule which contains the functional groups
     * @return a list of all functional groups (including "environments") extracted from the molecule
     */
    private List<IAtomContainer> extractGroups(IAtomContainer molecule) {
        if(isDbg()) LOGGING_TOOL.debug("########## Starting identification & extraction of functional groups... ##########");

        markedAtomToConnectedEnvCMap = new HashMap<IAtom, List<org.openscience.cdk.tools.ErtlFunctionalGroupsFinder.EnvironmentalC>>(molecule.getAtomCount());//Maps.newHashMapWithExpectedSize(molecule.getAtomCount());
        int[] atomIdxToFGMap = new int[molecule.getAtomCount()];
        Arrays.fill(atomIdxToFGMap, -1);
        int fGroupIdx = -1;

        while(!markedAtoms.isEmpty()) {
            // search for another functional group
            fGroupIdx++;

            // get next markedAtom as the starting node for the search
            int beginIdx = markedAtoms.iterator().next();
            if(isDbg()) LOGGING_TOOL.debug(String.format("Searching new functional group from atom #%d (%s)...", beginIdx,  molecule.getAtom(beginIdx).getSymbol()));

            // do a BFS from there
            Queue<Integer> queue = new ArrayDeque<>();
            queue.add(beginIdx);

            while(!queue.isEmpty()) {
                int currentIdx = queue.poll();

                // we are only interested in marked atoms that are not yet included in a group
                if(!markedAtoms.contains(currentIdx)){
                    continue;
                }

                // if it isn't...
                IAtom currentAtom = molecule.getAtom(currentIdx);
                if(isDbg()) LOGGING_TOOL.debug(String.format("	visiting marked atom: #%d (%s)", currentIdx, currentAtom.getSymbol()));

                // add its index to the functional group
                atomIdxToFGMap[currentIdx] = fGroupIdx;
                // also scratch the index from markedAtoms
                markedAtoms.remove(currentIdx);

                // and take look at the connected atoms
                List<EnvironmentalC> currentEnvironment = new ArrayList<>();
                for(int connectedIdx : adjList[currentIdx]) {
                    // add connected marked atoms to queue
                    if(markedAtoms.contains(connectedIdx)) {
                        queue.add(connectedIdx);
                        continue;
                    }

                    // ignore already handled connected atoms
                    if(atomIdxToFGMap[connectedIdx] >= 0){
                        continue;
                    }

                    // add unmarked connected aromatic heteroatoms
                    IAtom connectedAtom = molecule.getAtom(connectedIdx);
                    if(isHeteroatom(connectedAtom) && connectedAtom.isAromatic()) {
                        if(isDbg()) LOGGING_TOOL.debug("	   added connected aromatic heteroatom " + connectedAtom.getSymbol());
                        atomIdxToFGMap[connectedIdx] = fGroupIdx;
                        // note that this aromatic heteroatom has been added to a group
                        aromaticHeteroAtomIndicesToIsInGroupBoolMap.put(connectedIdx, true);
                    }

                    // add unmarked connected atoms to current marked atom's environment
                    IBond connectedBond = bondMap.get(currentIdx, connectedIdx);

                    EnvironmentalCType type;
                    if (connectedAtom.getAtomicNumber() == 6) {
                        if(connectedAtom.isAromatic())
                            type = EnvironmentalCType.C_AROMATIC;
                        else
                            type = EnvironmentalCType.C_ALIPHATIC;
                    }
                    else {
                        // aromatic heteroatom, so just ignore
                        continue;
                    }
                    currentEnvironment.add(new EnvironmentalC(type, connectedBond, connectedBond.getBegin() == connectedAtom ? 0 : 1));
                }
                markedAtomToConnectedEnvCMap.put(currentAtom, currentEnvironment);

                // debug logging
                if(isDbg()) {
                    int cAromCount = 0, cAliphCount = 0;
                    for(EnvironmentalC comp : currentEnvironment) {
                        if(comp.getType() == EnvironmentalCType.C_AROMATIC)
                            cAromCount++;
                        else if(comp.getType() == EnvironmentalCType.C_ALIPHATIC)
                            cAliphCount++;
                    }
                    LOGGING_TOOL.debug(String.format("	   logged marked atom's environment: C_ar:%d, C_al:%d (and %d implicit hydrogens)", cAromCount, cAliphCount, currentAtom.getImplicitHydrogenCount()));
                }
            }

            if(isDbg()) LOGGING_TOOL.debug("	search completed.");
        }

        // also create FG for lone aromatic heteroatoms, not connected to a FG yet.
        for(int atomIdx : aromaticHeteroAtomIndicesToIsInGroupBoolMap.keySet()) {
            if(!aromaticHeteroAtomIndicesToIsInGroupBoolMap.get(atomIdx)) {
                fGroupIdx++;
                atomIdxToFGMap[atomIdx] = fGroupIdx;
                if(isDbg()) LOGGING_TOOL.debug("Created FG for lone aromatic heteroatom: " + molecule.getAtom(atomIdx).getSymbol());
            }
        }

        List<IAtomContainer> fGs = partitionIntoGroups(molecule, atomIdxToFGMap, fGroupIdx + 1);

        if(isDbg()) LOGGING_TOOL.debug(String.format("########## Found & extracted %d functional groups. ##########", fGroupIdx + 1));
        return fGs;
    }

    /**
     * Generalizes the full environments of functional groups, providing a good balance between preserving
     * meaningful detail and generalization.
     *
     * @param fGroups the list of functional groups including "environments"
     */
    private void expandGeneralizedEnvironments(List<IAtomContainer> fGroups){
        if(isDbg()) LOGGING_TOOL.debug("########## Starting generalization of functional groups... ##########");

        for(IAtomContainer fGroup : fGroups) {
            int atomCount = fGroup.getAtomCount();

            if(isDbg()) LOGGING_TOOL.debug(String.format("Generalizing functional group (%d atoms)...", atomCount));

            // prechecking for special cases...
            if(fGroup.getAtomCount() == 1) {
                IAtom atom = fGroup.getAtom(0);
                List<EnvironmentalC> environment = markedAtomToConnectedEnvCMap.get(atom);

                if(environment != null) {
                    int envCCount = environment.size();

                    // for H2N-C_env & HO-C_env -> do not replace H & C_env by R!
                    if((atom.getAtomicNumber() == 8 && envCCount == 1)
                        || (atom.getAtomicNumber() == 7 && envCCount == 1)){
                            if(isDbg()) LOGGING_TOOL.debug(String.format("   - found single atomic N or O FG with one env. C. Expanding environment...", atom.getSymbol()));
                            expandEnvironment(atom, fGroup);

                            int hCount = atom.getImplicitHydrogenCount();
                            if(hCount != 0) {
                                if(isDbg()) LOGGING_TOOL.debug(String.format("   - adding %d hydrogens...", hCount));
                                addHydrogens(atom, hCount, fGroup);
                                atom.setImplicitHydrogenCount(0);
                            }
                            continue;
                    }
                    // for HN-(C_env)-C_env & HS-C_env -> do not replace H by R! (only C_env!)
                    if((atom.getAtomicNumber() == 7 && envCCount == 2)
                        || (atom.getAtomicNumber() == 16 && envCCount == 1)) {
                        if(isDbg()) LOGGING_TOOL.debug("   - found sec. amine or simple thiol");
                        int hCount = atom.getImplicitHydrogenCount();
                        if(hCount != 0) {
                            if(isDbg()) LOGGING_TOOL.debug(String.format("   - adding %d hydrogens...", hCount));
                            addHydrogens(atom, hCount, fGroup);
                            atom.setImplicitHydrogenCount(0);
                        }
                        if(isDbg()) LOGGING_TOOL.debug("   - expanding environment...");
                        expandEnvironmentGeneralized(atom, fGroup);
                        continue;
                    }
                }
                else if(isHeteroatom(atom)) {
                    int rAtomCount = atom.getValency();
                    Integer hCount = atom.getImplicitHydrogenCount();
                    if(hCount != null && hCount != 0) {
                        atom.setImplicitHydrogenCount(0);
                    }
                    String atomTypeName = atom.getAtomTypeName();
                    if(isDbg()) LOGGING_TOOL.debug(String.format("   - found single aromatic heteroatom (%s, Atomtype %s). Adding %d R-Atoms...", atom.getSymbol(), atomTypeName, rAtomCount));
                    addRAtoms(atom, rAtomCount, fGroup);
                    continue;
                }
            }

            // get atoms to process
            List<IAtom> fGroupAtoms = new ArrayList<IAtom>(fGroup.getAtomCount());//Lists.newArrayList(fGroup.atoms());
            fGroup.atoms().forEach(fGroupAtoms::add);

            // process atoms...
            for(IAtom atom : fGroupAtoms) {
                List<EnvironmentalC> environment = markedAtomToConnectedEnvCMap.get(atom);

                if(environment == null) {
                    if(atom.getImplicitHydrogenCount() != 0) {
                        atom.setImplicitHydrogenCount(0);
                    }
                    int rAtomCount = atom.getValency() - 1;
                    if(isDbg()) LOGGING_TOOL.debug(String.format("   - found connected aromatic heteroatom (%s). Adding %d R-Atoms...", atom.getSymbol(), rAtomCount));
                    addRAtoms(atom, rAtomCount, fGroup);
                }

                // processing carbons...
                if(atom.getAtomicNumber() == 6) {
                    if(atom.getProperty(CARBONYL_C_MARKER) == null) {
                        if(atom.getImplicitHydrogenCount() != 0) {
                            atom.setImplicitHydrogenCount(0);
                        }
                        if(isDbg()) LOGGING_TOOL.debug("   - ignoring environment for marked carbon atom");
                        continue;
                    }
                    else {
                        if(isDbg()) LOGGING_TOOL.debug("   - found carbonyl-carbon. Expanding environment...");
                        expandEnvironmentGeneralized(atom, fGroup);
                        continue;
                    }
                }
                // processing heteroatoms...
                else {
                    if(isDbg()) LOGGING_TOOL.debug(String.format("   - found heteroatom (%s). Expanding environment...", atom.getSymbol()));
                    expandEnvironmentGeneralized(atom, fGroup);
                    continue;
                }
            }
        }

        if(isDbg()) LOGGING_TOOL.debug("########## Generalization of functional groups completed. ##########");
    }

    /**
     * Expands the full environments of functional groups, converted into atoms and bonds.
     *
     * @param fGroups the list of functional groups including "environments"
     */
    private void expandFullEnvironments(List<IAtomContainer> fGroups) {
        if(isDbg()) LOGGING_TOOL.debug("########## Starting expansion of full environments for functional groups... ##########");

        for(IAtomContainer fGroup : fGroups) {
            int atomCount = fGroup.getAtomCount();
            if(isDbg()) LOGGING_TOOL.debug(String.format("Expanding environment on functional group (%d atoms)...", atomCount));

            for(int i = 0; i < atomCount; i++) {
                IAtom atom = fGroup.getAtom(i);

                if(isDbg()) LOGGING_TOOL.debug(String.format(" - Atom #%d   - Expanding environment...", i));
                expandEnvironment(atom, fGroup);

                int hCount = atom.getImplicitHydrogenCount();
                if(hCount != 0) {
                    if(isDbg()) LOGGING_TOOL.debug(String.format("   - adding %d hydrogens...", hCount));
                    addHydrogens(atom, hCount, fGroup);
                    atom.setImplicitHydrogenCount(0);
                }
            }
        }

        if(isDbg()) LOGGING_TOOL.debug("########## Expansion of full environments for functional groups completed. ##########");
    }

    private void expandEnvironment(IAtom atom, IAtomContainer container) {
        List<EnvironmentalC> environment = markedAtomToConnectedEnvCMap.get(atom);

        if(environment == null || environment.isEmpty()) {
            if(isDbg()) LOGGING_TOOL.debug("		found no environment to expand.");
            return;
        }

        int cAromCount = 0, cAliphCount = 0;
        for(EnvironmentalC envC : environment) {
            IAtom cAtom = atom.getBuilder().newInstance(IAtom.class, "C");
            cAtom.setAtomTypeName("C");
            cAtom.setImplicitHydrogenCount(0);
            if(envC.getType() == EnvironmentalCType.C_AROMATIC) {
                cAtom.setIsAromatic(true);
                cAromCount++;
            }
            else {
                cAliphCount++;
            }

            IBond bond = envC.createBond(atom, cAtom);

            container.addAtom(cAtom);
            container.addBond(bond);
        }

        if(isDbg()) LOGGING_TOOL.debug(String.format("		expanded environment: %dx C_ar and %dx C_al", cAromCount, cAliphCount));
    }

    // only call this on marked heteroatoms / carbonyl-C's!
    private void expandEnvironmentGeneralized(IAtom atom, IAtomContainer container) {

        List<EnvironmentalC> environment = markedAtomToConnectedEnvCMap.get(atom);

        if(environment == null) {
            if(isDbg()) LOGGING_TOOL.debug("		found no environment to expand.");
            return;
        }

        int rAtomCount = environment.size();
        int rAtomsForCCount = rAtomCount;
        if(atom.getAtomicNumber() == 8 && atom.getImplicitHydrogenCount() == 1) {
            addHydrogens(atom, 1, container);
            atom.setImplicitHydrogenCount(0);
            if(isDbg()) LOGGING_TOOL.debug("		expanded hydrogen on connected OH-Group");
        }
        else if(isHeteroatom(atom)) rAtomCount += atom.getImplicitHydrogenCount();
        addRAtoms(atom, rAtomCount, container);

        if(atom.getImplicitHydrogenCount() != 0) {
            atom.setImplicitHydrogenCount(0);
        }

        if(isDbg()) LOGGING_TOOL.debug(String.format("		expanded environment: %dx R-atom (incl. %d for H replacement)", rAtomCount, rAtomCount - rAtomsForCCount));
    }

    private static final boolean isHeteroatom(IAtom atom) {
        int atomicNr = atom.getAtomicNumber();
        return atomicNr != 1 && atomicNr != 6;
    }

    private final boolean  isNonmetal(IAtom atom) {
        return NONMETAL_ATOMIC_NUMBERS.contains(atom.getAtomicNumber());
    }

    private void addHydrogens(IAtom atom, int number, IAtomContainer container) {
        for(int i = 0; i < number; i++) {
            IAtom hydrogen = atom.getBuilder().newInstance(IAtom.class, "H");
            hydrogen.setAtomTypeName("H");
            hydrogen.setImplicitHydrogenCount(0);

            container.addAtom(hydrogen);
            container.addBond(atom.getBuilder().newInstance(IBond.class, atom, hydrogen, Order.SINGLE));
        }
    }

    private void addRAtoms(IAtom atom, int number, IAtomContainer container) {
        for(int i = 0; i < number; i++) {
            IPseudoAtom rAtom = atom.getBuilder().newInstance(IPseudoAtom.class, "R");
            rAtom.setAttachPointNum(1);
            rAtom.setImplicitHydrogenCount(0);

            container.addAtom(rAtom);
            container.addBond(atom.getBuilder().newInstance(IBond.class, atom, rAtom, Order.SINGLE));
        }
    }

    private List<IAtomContainer> partitionIntoGroups(IAtomContainer sourceContainer, int[] atomIdxToFGMap, int fGroupCount) {
        List<IAtomContainer> groups = new ArrayList<>(fGroupCount);
        for(int i = 0; i < fGroupCount; i++) {
            groups.add(sourceContainer.getBuilder().newInstance(IAtomContainer.class));
        }

        Map<IAtom, IAtomContainer> atomtoFGMap = new HashMap<IAtom, IAtomContainer>(sourceContainer.getAtomCount());//Maps.newHashMapWithExpectedSize(sourceContainer.getAtomCount());

        // atoms
        for(int atomIdx = 0; atomIdx < sourceContainer.getAtomCount(); atomIdx++) {
            int fGroupId = atomIdxToFGMap[atomIdx];

            if(fGroupId == -1) {
                continue;
            }

            IAtom atom = sourceContainer.getAtom(atomIdx);
            IAtomContainer myGroup = groups.get(fGroupId);
            myGroup.addAtom(atom);
            atomtoFGMap.put(atom, myGroup);
        }

        // bonds
        for(IBond bond : sourceContainer.bonds()) {
            IAtomContainer beginGroup = atomtoFGMap.get(bond.getBegin());
            IAtomContainer endGroup = atomtoFGMap.get(bond.getEnd());

            if(beginGroup == null || endGroup == null || beginGroup != endGroup)
                continue;

            beginGroup.addBond(bond);
        }

        // single electrons
        for (ISingleElectron electron : sourceContainer.singleElectrons()) {
            IAtomContainer group = atomtoFGMap.get(electron.getAtom());
            if(group != null)
                group.addSingleElectron(electron);
        }

        // lone pairs
        for (ILonePair lonePair : sourceContainer.lonePairs()) {
            IAtomContainer group = atomtoFGMap.get(lonePair.getAtom());
            if(group != null)
                group.addLonePair(lonePair);
        }

        return groups;
    }

    private static boolean isDbg() {
        return ErtlFunctionalGroupsFinder.LOGGING_TOOL.isDebugEnabled();
    }

    private boolean checkConstraints(IAtomContainer molecule) {
        for(IAtom atom : molecule.atoms()) {
            if(atom.getFormalCharge() != null && atom.getFormalCharge() != 0) {
                throw new IllegalArgumentException("Input molecule must not contain any charges.");
            }
            if(!isNonmetal(atom)) {
                throw new IllegalArgumentException("Input molecule must not contain metals or metalloids.");
            }
        }

        ConnectedComponents cc = new ConnectedComponents(adjList);
        if(cc.nComponents() != 1) {
            throw new IllegalArgumentException("Input molecule must consist of only a single connected structure.");
        }

        return true;
    }
}
