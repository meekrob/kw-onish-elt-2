package edu.unc.genomics.wigmath;

import java.io.IOException;
import java.nio.file.Path;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.beust.jcommander.Parameter;

import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.Interval;
import edu.unc.genomics.ReadablePathValidator;
import edu.unc.genomics.FloatValidator;
import edu.unc.genomics.WigMathTool;
import edu.unc.genomics.io.WigFileReader;
import edu.unc.genomics.io.WigFileException;

/**
 * Custom Equation 1.
 * This is equivalent to log(signal - input), but it uses an absolute min value
 * in place of missing data or cases where signal - input <= 0
 */

/**
 * Based on:
 * Subtract two (Big)Wig files
 * by
 * @author timpalpant
 *
 */
public class CustomEquation_1 extends WigMathTool {

  private static final Logger log = Logger.getLogger(CustomEquation_1.class);

  @Parameter(names = { "-m", "--minuend" }, description = "Minuend (top - file 1)", required = true, validateWith = ReadablePathValidator.class)
  public Path minuendFile;
  @Parameter(names = { "-s", "--subtrahend" }, description = "Subtrahend (bottom - file 2)", required = true, validateWith = ReadablePathValidator.class)
  public Path subtrahendFile;
  @Parameter(names = { "-z", "--absolute-min"}, description = "Use this value in place of NaN for log(minuend - subtrahend <= 0)", required = true, validateWith = FloatValidator.class)
  public float absolute_min;


  public boolean assumeZero = false;

  WigFileReader minuendReader, subtrahendReader;

  @Override
  public void setup() {
    log.setLevel((Level)Level.INFO); // defaults to DEBUG
    log.debug("Initializing input files");
    try {
      minuendReader = WigFileReader.autodetect(minuendFile);
      subtrahendReader = WigFileReader.autodetect(subtrahendFile);
    } catch (IOException e) {
      throw new CommandLineToolException(e);
    }
    inputs.add(minuendReader);
    inputs.add(subtrahendReader);
    log.debug("Initialized " + inputs.size() + " input files");
    log.info("Using " + absolute_min + " as absolute minimum");
  }

  @Override
  public float[] compute(Interval chunk) throws IOException, WigFileException {
    float observed_min = 0;
    log.setLevel((Level)Level.INFO); // defaults to DEBUG
    float[] minuend = minuendReader.query(chunk).getValues();
    float[] subtrahend = subtrahendReader.query(chunk).getValues();

    // Fill missing data with zeros
    if (assumeZero) {
      for (int i = 0; i < minuend.length; i++) {
        if (Float.isNaN(minuend[i])) {
          minuend[i] = 0;
        }
        if (Float.isNaN(subtrahend[i])) {
          subtrahend[i] = 0;
        }
      }
    }

    for (int i = 0; i < minuend.length; i++) {
      if (minuend[i] <= subtrahend[i]) {
        minuend[i] = absolute_min;
      }
      else { 
        minuend[i] = (float)Math.log(minuend[i] - subtrahend[i]);
        if (minuend[i] < observed_min) {
          observed_min = minuend[i];
        }
      }

      if (minuend[i] < absolute_min) {
        log.warn("Encountered a value less than absolute minimum (" + absolute_min + "): " + minuend[i]);
      }
    }
    log.info("Compute chunk: observed min was: " + observed_min);
    return minuend;
  }

  /**
   * @param args
   * @throws WigFileException
   * @throws IOException
   */
  public static void main(String[] args) throws IOException, WigFileException {
    new CustomEquation_1().instanceMain(args);
  }

}
