/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

/**
 * Thrown where an exact ring basis is required but could not be established, for example by a
 * fingerprint term or SMARTS key that counts ring memberships. It is an
 * {@link IllegalStateException}, so existing handlers of that type still catch it.
 */
public class RingBasisUnavailableException extends IllegalStateException {
  private static final long serialVersionUID = 1L;

  /** The message, the same in C++, Java and Python. */
  static final String MESSAGE =
      "Exact ring basis unavailable: ring perception reached its work limit before the minimum cycle basis was complete";

  /** An exception with the standard message. */
  public RingBasisUnavailableException() {
    super(MESSAGE);
  }
}
