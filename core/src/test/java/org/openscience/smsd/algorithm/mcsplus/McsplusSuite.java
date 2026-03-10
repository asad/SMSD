/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.mcsplus;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.platform.suite.api.SelectClasses;
import org.junit.platform.suite.api.Suite;


/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
@Suite
@SelectClasses({org.openscience.smsd.algorithm.mcsplus.MCSPlusHandlerTest.class,
    org.openscience.smsd.algorithm.mcsplus.IsomorphismMCSPlusTest.class})
public class McsplusSuite {

    @BeforeAll
    public static void setUpClass() throws Exception {
    }

    @AfterAll
    public static void tearDownClass() throws Exception {
    }

    @BeforeEach 

    public void setUp() throws Exception {
    }

    @AfterEach 

    public void tearDown() throws Exception {
    }
}





