/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.cmd;

/**
 * Legacy entry point kept for backward compatibility with older scripts.
 * Prints a deprecation notice and forwards to the new CLI:
 *   com.bioinception.smsd.cli.SMSDcli
 *
 * No imports of the old class name are used to avoid compile errors.
 */
public class SMSDcmd {
    public static void main(String[] args) {
        System.err.println("[DEPRECATED] uk.ac.ebi.smsd.cmd.SMSDcmd is legacy. "
                + "Please use com.bioinception.smsd.cli.SMSDcli (same options).");
        // Forward to the new CLI
        com.bioinception.smsd.cli.SMSDcli.main(args);
    }
}

