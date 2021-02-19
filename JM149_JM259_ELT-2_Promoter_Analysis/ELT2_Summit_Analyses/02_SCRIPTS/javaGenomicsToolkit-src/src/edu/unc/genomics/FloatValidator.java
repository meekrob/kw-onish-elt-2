/**
 * 
 */
package edu.unc.genomics;

import java.nio.file.Files;
import java.nio.file.Path;

import com.beust.jcommander.IParameterValidator;
import com.beust.jcommander.ParameterException;

/**
 * @author meekrob based on @author timpalpant
 *
 */
public class FloatValidator implements IParameterValidator {

  /*
   * (non-Javadoc)
   * 
   * @see com.beust.jcommander.IParameterValidator#validate(java.lang.String,
   * java.lang.String)
   */
  @Override
  public void validate(String name, String value) throws ParameterException {
    try {
        Float.parseFloat(value);
    }
    catch (NumberFormatException ex) {
        throw new ParameterException("Parameter " + name + " should be a parseable floating point value. [" + value + "] given.");
    }
  }

}
