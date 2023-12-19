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

package org.openscience.cdk.tools.efgf.app;

/**
 * Main class for starting application.
 *
 * @author Jonas Schaub
 * @version 1.2
 */
public class Main {

    private Main() {

    }

    /**
     * Starts the application. Command line arguments must be the name of an SD-file to read (must be located in the
     * same directory as the application's .jar file) and the number of different threads to use for calculation.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            ErtlFunctionalGroupsFinderPerformanceSnapshotApp tmpApplication = new ErtlFunctionalGroupsFinderPerformanceSnapshotApp(args);
        } catch (Exception anException) {
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }

}
