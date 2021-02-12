package edu.unc.genomics.visualization;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.io.IOException;
import java.nio.file.Path;

import org.apache.log4j.Logger;

import com.beust.jcommander.Parameter;

import edu.ucsc.genome.TrackHeader;
import edu.unc.genomics.CommandLineTool;
import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.ReadablePathValidator;
import edu.unc.genomics.BedEntry;
import edu.unc.genomics.io.BedFileReader;
import edu.unc.genomics.Interval;
/**
 * Test class for looping over a file of intervals
*/

public class CatBEDlike extends CommandLineTool {
  private static final Logger log = Logger.getLogger(CatBEDlike.class);

  @Parameter(names = { "-i", "--input" }, description = "BED or BED-like file", required = true, validateWith = ReadablePathValidator.class)
  public Path intervalFile;


  private List<BedEntry> intervals;

  @Override
  public final void run() throws IOException {
    log.debug("Here we go!!!");
    TrackHeader header = TrackHeader.newWiggle();
    header.setName("Processed " + intervalFile.getFileName());
    header.setDescription("Processed " + intervalFile.getFileName());
    try (BedFileReader reader = new BedFileReader(intervalFile)) {
        intervals = reader.loadAll();
        int i = 0;
        for (BedEntry interval : intervals ) {
            System.out.printf("line %1$5d: %2$s%n", ++i, interval.toBed());
        }
    }

  }
  public static void main(String[] args) {
    new CatBEDlike().instanceMain(args);
  }
}
