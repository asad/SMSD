/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.tools;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * test-smsd
 */
public class TimeManagerTest {

    @Test
    public void testTimeManager() throws Exception {
        TimeManager tMan = new TimeManager();
        Assert.assertNotNull(tMan);
    }

    /**
     * Test of getElapsedTimeInHours method, of class TimeManager.
     */
    @Test
    public void testGetElapsedTimeInHours() {
        //////System.out.println("getElapsedTimeInHours");
        TimeManager instance = new TimeManager();
        double expResult = 0.0001;
        myMethod(360);
        double result = instance.getElapsedTimeInHours();
        Assert.assertEquals(expResult, result, 0.0001);
    }

    /**
     * Test of getElapsedTimeInMinutes method, of class TimeManager.
     */
    @Test
    public void testGetElapsedTimeInMinutes() {
        //////System.out.println("getElapsedTimeInMinutes");
        TimeManager instance = new TimeManager();
        double expResult = 0.006;
        myMethod(360);
        double result = instance.getElapsedTimeInMinutes();
        Assert.assertEquals(expResult, result, 0.006);
    }

    /**
     * Test of getElapsedTimeInSeconds method, of class TimeManager.
     */
    @Test
    public void testGetElapsedTimeInSeconds() {
        //////System.out.println("getElapsedTimeInSeconds");
        TimeManager instance = new TimeManager();
        double expResult = 0.36;
        myMethod(360);
        double result = instance.getElapsedTimeInSeconds();
        Assert.assertEquals(expResult, result, 0.36);
    }

    /**
     * Test of getElapsedTimeInMilliSeconds method, of class TimeManager.
     */
    @Test
    public void testGetElapsedTimeInMilliSeconds() {
        //////System.out.println("getElapsedTimeInMilliSeconds");
        TimeManager instance = new TimeManager();
        double expResult = 360;
        myMethod(360);
        double result = instance.getElapsedTimeInMilliSeconds();
        Assert.assertEquals(expResult, result, 360);
    }

    public void myMethod(long timeMillis) {
        //////System.out.println("Starting......");

        // pause for a while
        Thread thisThread = Thread.currentThread();
        try {
            thisThread.sleep(timeMillis);
        } catch (Throwable t) {
            throw new OutOfMemoryError("An Error has occured");
        }
        //////System.out.println("Ending......");

    }
}
