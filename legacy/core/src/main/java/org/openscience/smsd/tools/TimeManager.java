/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.tools;

import java.text.SimpleDateFormat;
import java.util.TimeZone;

/**
 * Class that handles execution time of the MCS search.
 *
 * long diffSeconds = time / 1000; long diffMinutes = time / (60 * 1000); long
 * diffHours = time / (60 * 60 * 1000); long diffDays = time / (24 * 60 * 60 *
 * 1000);
 *
 *
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class TimeManager {

    private final double startTime;
    private final SimpleDateFormat dateFormat;

    /**
     * Constructor for storing execution time
     */
    public TimeManager() {

        dateFormat = new SimpleDateFormat("HH:mm:ss");
        dateFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
        startTime = System.currentTimeMillis();
    }

    /**
     * Returns Elapsed Time In Hours
     *
     * @return Elapsed Time In Hours
     */
    public synchronized double getElapsedTimeInHours() {
        double currentTime = System.currentTimeMillis();
        return (currentTime - startTime) / (60 * 60 * 1000);

    }

    /**
     * Returns Elapsed Time In Minutes
     *
     * @return Elapsed Time In Minutes
     */
    public synchronized double getElapsedTimeInMinutes() {
        double currentTime = System.currentTimeMillis();
        return (currentTime - startTime) / (60 * 1000);

    }

    /**
     * Return Elapsed Time In Seconds
     *
     * @return Elapsed Time In Seconds
     */
    public synchronized double getElapsedTimeInSeconds() {
        double currentTime = System.currentTimeMillis();
        return ((currentTime - startTime) / 1000);

    }

    /**
     * Returns Elapsed Time In Mill Seconds
     *
     * @return Elapsed Time In Mill Seconds
     */
    public synchronized double getElapsedTimeInMilliSeconds() {
        double currentTime = System.currentTimeMillis();
        return (currentTime - startTime);

    }
}
