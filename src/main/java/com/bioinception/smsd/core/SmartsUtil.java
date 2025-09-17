/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

public class SmartsUtil {
    public static List<int[]> matchAll(String smarts, IAtomContainer target) throws Exception {
        try {
            Class<?> sp = Class.forName("org.openscience.cdk.smarts.SmartsPattern");
            Method create = sp.getMethod("create", String.class);
            Object patt = create.invoke(null, smarts);
            Method matchAll = sp.getMethod("matchAll", IAtomContainer.class);
            Object mappings = matchAll.invoke(patt, target);
            Method uniqueAtoms = mappings.getClass().getMethod("uniqueAtoms");
            Object uniq = uniqueAtoms.invoke(mappings);
            Method toArray = uniq.getClass().getMethod("toArray");
            return toIntArrays((int[][]) toArray.invoke(uniq));
        } catch (ClassNotFoundException e) {
            try {
                Class<?> sp = Class.forName("org.openscience.cdk.smiles.smarts.SmartsPattern");
                Method create = sp.getMethod("create", String.class);
                Object patt = create.invoke(null, smarts);
                Method matchAll = sp.getMethod("matchAll", IAtomContainer.class);
                Object mappings = matchAll.invoke(patt, target);
                Method uniqueAtoms = mappings.getClass().getMethod("uniqueAtoms");
                Object uniq = uniqueAtoms.invoke(mappings);
                Method toArray = uniq.getClass().getMethod("toArray");
                return toIntArrays((int[][]) toArray.invoke(uniq));
            } catch (ClassNotFoundException e2) {
                Class<?> sqt = Class.forName("org.openscience.cdk.smiles.smarts.SMARTSQueryTool");
                Object tool = sqt.getConstructor(String.class).newInstance(smarts);
                Method matches = sqt.getMethod("matches", IAtomContainer.class);
                boolean ok = (Boolean) matches.invoke(tool, target);
                if (!ok) return new ArrayList<int[]>();
                Method getUniqueMatchingAtoms = sqt.getMethod("getUniqueMatchingAtoms");
                @SuppressWarnings("unchecked")
                List<List<Integer>> idxs = (List<List<Integer>>) getUniqueMatchingAtoms.invoke(tool);
                List<int[]> out = new ArrayList<int[]>(idxs.size());
                for (List<Integer> lst : idxs) {
                    int[] arr = new int[lst.size()];
                    for (int i=0;i<lst.size();i++) arr[i] = lst.get(i).intValue();
                    out.add(arr);
                }
                return out;
            }
        }
    }
    private static List<int[]> toIntArrays(int[][] arr) {
        List<int[]> out = new ArrayList<int[]>(arr.length);
        for (int[] a : arr) out.add(a);
        return out;
    }
}
