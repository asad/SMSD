/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.core;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Named SMARTS predicate registry: $name -> SMARTS.
 * Minimal JSON loader for a flat object: { "name": "SMARTS", ... }.
 */
public class PredicateRegistry {
    private final Map<String, String> map = new LinkedHashMap<>();

    public PredicateRegistry() {
        map.put("isAmideN",       "[$([NX3;H2,H1;!$(NC=O)]),$([NX3;$(NC=O)])]");
        map.put("isCarboxylC",    "[CX3](=O)[OX2H1,OX1-]");
        map.put("isSulfonamideN", "[NX3;$(NS(=O)=O)]");
        map.put("isEster",        "[#6][CX3](=O)[OX2][#6]");
        map.put("isKetone",       "[#6][CX3](=O)[#6]");
        map.put("isHalogen",      "[F,Cl,Br,I]");
        map.put("isPositive",     "[*+]");
        map.put("isNegative",     "[*-]");
        map.put("isHBDonor",      "[N!H0,O!H0,S!H0]");
        map.put("isHBAcceptor",   "[N,O,S;H0;!+0]");
    }

    public void put(String name, String smarts) { map.put(name, smarts); }
    public boolean contains(String name) { return map.containsKey(name); }
    public String get(String name) { return map.get(name); }
    public Map<String, String> snapshot() { return new LinkedHashMap<>(map); }

    /** Load a tiny JSON object from file: {"name":"SMARTS", ...}. */
    public void loadJsonFile(String path) throws Exception {
        String txt = new String(Files.readAllBytes(Paths.get(path)), StandardCharsets.UTF_8).trim();
        if (txt.isEmpty() || txt.charAt(0) != '{')
            throw new IllegalArgumentException("Predicate pack must be a JSON object {name: SMARTS}");
        int i = 1, n = txt.length();

        while (i < n) {
            i = skipWs(txt, i);
            if (i >= n) break;
            if (txt.charAt(i) == '}') { i++; break; }

            if (txt.charAt(i) != '\"') throw new IllegalArgumentException("Invalid JSON (key must start with \")");
            Parse keyP = readJsonString(txt, i + 1);
            String key = keyP.text;
            i = skipWs(txt, keyP.next);
            if (i >= n || txt.charAt(i) != ':') throw new IllegalArgumentException("Invalid JSON (missing :)");
            i = skipWs(txt, i + 1);

            if (i >= n || txt.charAt(i) != '\"') throw new IllegalArgumentException("Invalid JSON (value must start with \")");
            Parse valP = readJsonString(txt, i + 1);
            String value = valP.text;
            i = skipWs(txt, valP.next);

            map.put(key, value);

            if (i < n && txt.charAt(i) == ',') { i++; continue; }
            if (i < n && txt.charAt(i) == '}') { i++; break; }
        }
    }

    private static int skipWs(String s, int i) {
        int n = s.length();
        while (i < n) {
            char c = s.charAt(i);
            if (c == ' ' || c == '\n' || c == '\r' || c == '\t') i++;
            else break;
        }
        return i;
    }

    private static class Parse {
        final String text; final int next;
        Parse(String text, int next) { this.text = text; this.next = next; }
    }

    private static Parse readJsonString(String s, int i) {
        StringBuilder out = new StringBuilder();
        int n = s.length();
        while (i < n) {
            char c = s.charAt(i);
            if (c == '\\') {
                if (i + 1 >= n) throw new IllegalArgumentException("Invalid JSON escape");
                char e = s.charAt(i + 1);
                switch (e) {
                    case '\"': out.append('\"'); break;
                    case '\\': out.append('\\'); break;
                    case 'n':  out.append('\n'); break;
                    case 'r':  out.append('\r'); break;
                    case 't':  out.append('\t'); break;
                    default:   out.append(e);   break;
                }
                i += 2; continue;
            }
            if (c == '\"') return new Parse(out.toString(), i + 1);
            out.append(c); i++;
        }
        throw new IllegalArgumentException("Unterminated JSON string");
    }
}
