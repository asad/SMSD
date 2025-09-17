/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd.core;

/** Expand $name tokens to recursive SMARTS using a registry, leaving existing $() intact. */
public class SmartPredicateProcessor {
    public static String expandNamedPredicates(String smarts, PredicateRegistry reg) {
        StringBuilder out = new StringBuilder();
        int i=0, n=smarts.length();
        while (i<n) {
            char c = smarts.charAt(i);
            if (c=='[') {
                int depth=1, j=i+1;
                while (j<n && depth>0) {
                    char d = smarts.charAt(j);
                    if (d=='[') depth++; else if (d==']') depth--;
                    j++;
                }
                String content = smarts.substring(i+1, j-1);
                out.append('[').append(rewriteContent(content, reg)).append(']');
                i = j;
            } else {
                out.append(c); i++;
            }
        }
        return out.toString();
    }

    private static String rewriteContent(String content, PredicateRegistry reg) {
        StringBuilder sb = new StringBuilder();
        int i=0, n=content.length();
        while (i<n) {
            char c = content.charAt(i);
            if (c=='$') {
                if (i+1<n && content.charAt(i+1)=='(') {
                    // copy $(...)
                    int depth=1, j=i+2;
                    while (j<n && depth>0) {
                        char d = content.charAt(j);
                        if (d=='(') depth++; else if (d==')') depth--;
                        j++;
                    }
                    sb.append(content, i, j);
                    i = j;
                } else {
                    int j=i+1;
                    while (j<n) {
                        char d = content.charAt(j);
                        if (Character.isLetterOrDigit(d) || d=='_') j++; else break;
                    }
                    String name = content.substring(i+1, j);
                    if (!name.isEmpty() && reg.contains(name)) {
                        sb.append("$(").append(reg.get(name)).append(")");
                        i = j;
                    } else {
                        sb.append(c); i++;
                    }
                }
            } else { sb.append(c); i++; }
        }
        return sb.toString();
    }
}
