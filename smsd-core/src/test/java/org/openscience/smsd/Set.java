/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd;

import java.util.TreeSet;

/**
 *
 * @author Asad
 */
public class Set {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        TreeSet<Integer> test1 = new TreeSet<>();

        test1.add(1);
        test1.add(2);
        test1.add(3);
        test1.add(6);

        TreeSet<Integer> test2 = new TreeSet<>();
        test2.add(1);
        test2.add(2);
        test2.add(3);
        test2.add(4);
        test2.add(5);
        TreeSet<Integer> common = (TreeSet) test1.clone();
        TreeSet<Integer> diff = (TreeSet) test1.clone();
        ////System.out.println("common " + common.retainAll(test2));
        ////System.out.println("diff " + diff.removeAll(test2));
        ////System.out.println("diff elements" + diff);
        ////System.out.println("common elements" + common);

    }

}
