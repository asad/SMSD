/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import static org.junit.jupiter.api.Assertions.*;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.jar.Attributes;
import java.util.jar.JarEntry;
import java.util.jar.JarOutputStream;
import java.util.jar.Manifest;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import picocli.CommandLine;

/** The CLI reports the release version recorded in its JAR manifest. */
class CLIVersionRegressionTest {
  @TempDir Path directory;

  @Test
  void versionLineNamesTheManifestVersion() {
    assertEquals("SMSD Pro 9.8.7 — com.bioinceptionlabs", SMSDcli.versionLine("9.8.7"));
    assertEquals("SMSD Pro 9.8.7 — com.bioinceptionlabs", SMSDcli.versionLine(" 9.8.7 "));
    assertEquals("SMSD Pro (development build) — com.bioinceptionlabs", SMSDcli.versionLine(null));
    assertEquals("SMSD Pro (development build) — com.bioinceptionlabs", SMSDcli.versionLine(""));
  }

  @Test
  void versionOptionPrintsTheProviderLine() {
    StringWriter out = new StringWriter();
    int status = new CommandLine(new SMSDcli()).setOut(new PrintWriter(out, true))
        .setErr(new PrintWriter(new StringWriter())).execute("--version");
    assertEquals(0, status);
    assertEquals(new SMSDcli.ManifestVersion().getVersion()[0], out.toString().trim());
  }

  /** Loads the CLI classes from a JAR whose manifest names a version, as the release JAR does. */
  @Test
  void versionIsReadFromTheJarManifest() throws Exception {
    Path classes = Paths.get(SMSDcli.class.getProtectionDomain().getCodeSource().getLocation().toURI());
    Path jar = directory.resolve("cli.jar");
    Manifest manifest = new Manifest();
    manifest.getMainAttributes().put(Attributes.Name.MANIFEST_VERSION, "1.0");
    manifest.getMainAttributes().put(Attributes.Name.IMPLEMENTATION_VERSION, "9.8.7");
    Path cli = classes.resolve("com/bioinception/smsd/cli");
    List<Path> files;
    try (Stream<Path> listing = Files.list(cli)) {
      files = listing.filter(path -> path.toString().endsWith(".class")).collect(Collectors.toList());
    }
    try (OutputStream file = Files.newOutputStream(jar); JarOutputStream out = new JarOutputStream(file, manifest)) {
      for (Path path : files) {
        out.putNextEntry(new JarEntry("com/bioinception/smsd/cli/" + path.getFileName()));
        out.write(Files.readAllBytes(path));
        out.closeEntry();
      }
    }
    try (URLClassLoader loader = new CliFirstLoader(jar.toUri().toURL(), getClass().getClassLoader())) {
      Class<?> provider = loader.loadClass("com.bioinception.smsd.cli.SMSDcli$ManifestVersion");
      assertNotSame(SMSDcli.ManifestVersion.class, provider);
      CommandLine.IVersionProvider version =
          (CommandLine.IVersionProvider) provider.getDeclaredConstructor().newInstance();
      assertArrayEquals(new String[] {"SMSD Pro 9.8.7 — com.bioinceptionlabs"}, version.getVersion());
    }
  }

  /**
   * Loads the CLI package from the JAR and every other class from the parent. The CLI package is
   * the one this loader defines from the JAR manifest: Java 8 looks a package up in the parent
   * first, which defines it from the class directory, without a manifest.
   */
  private static final class CliFirstLoader extends URLClassLoader {
    private static final String CLI_PACKAGE = "com.bioinception.smsd.cli";
    private Package cliPackage;

    CliFirstLoader(URL jar, ClassLoader parent) {
      super(new URL[] {jar}, parent);
    }

    @Override
    @SuppressWarnings("deprecation")
    protected Package getPackage(String name) {
      return CLI_PACKAGE.equals(name) ? cliPackage : super.getPackage(name);
    }

    @Override
    protected Package definePackage(String name, Manifest manifest, URL url) {
      Package defined = super.definePackage(name, manifest, url);
      if (CLI_PACKAGE.equals(name)) cliPackage = defined;
      return defined;
    }

    @Override
    protected Class<?> loadClass(String name, boolean resolve) throws ClassNotFoundException {
      synchronized (getClassLoadingLock(name)) {
        Class<?> found = findLoadedClass(name);
        if (found == null && name.startsWith("com.bioinception.smsd.cli.")) found = findClass(name);
        if (found == null) return super.loadClass(name, resolve);
        if (resolve) resolveClass(found);
        return found;
      }
    }
  }
}
