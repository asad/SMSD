/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import static org.junit.jupiter.api.Assertions.*;

class ChemistryProfileRegressionTest {
    @ParameterizedTest
    @ValueSource(strings = {"strict", "STRICT", "compat-substruct", "Compat-FMCS", "default"})
    void basicProfilesKeepChargeAndRingMatchingOptional(String name) {
        ChemOptions options = ChemOptions.profile(name);
        assertFalse(options.matchFormalCharge);
        assertFalse(options.ringMatchesRingOnly);
    }

    @Test
    void specialisedProfilesApplyTheirDeclaredConstraints() {
        ChemOptions pharma = ChemOptions.profile("pharma");
        assertTrue(pharma.matchFormalCharge);
        assertTrue(pharma.ringMatchesRingOnly);
        assertTrue(pharma.completeRingsOnly);
        assertEquals(ChemOptions.BondOrderMode.LOOSE, ChemOptions.profile("reaction").matchBondOrder);
    }

    @Test
    void unknownProfileIsRejected() {
        assertThrows(IllegalArgumentException.class, () -> ChemOptions.profile("strcit"));
    }
}
