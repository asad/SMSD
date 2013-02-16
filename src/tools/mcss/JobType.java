/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tools.mcss;

/**
 *
 * @author Asad
 */
public enum JobType implements Comparable<JobType> {

    /**
     * Default MCS algorithm.
     */
    MCS(0, "MCS search"),
    /**
     * Substructure search algorithm.
     */
    Substructure(1, "Substructure search");
    private final int type;
    private final String description;

    JobType(int aStatus, String desc) {
        this.type = aStatus;
        this.description = desc;
    }

    /**
     * Returns type of algorithm.
     *
     * @return type of algorithm
     */
    public int type() {
        return this.type;
    }

    /**
     * Returns short description of the algorithm.
     *
     * @return description of the algorithm
     */
    public String description() {
        return this.description;
    }
}
