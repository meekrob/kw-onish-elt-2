package edu.unc.genomics.ngs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*; // for Path, Paths, newBufferedWriter

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import com.beust.jcommander.Parameter;

import edu.unc.genomics.CommandLineTool;
import edu.unc.genomics.CommandLineToolException;
import edu.unc.genomics.Interval;
import edu.unc.genomics.Contig;
import edu.unc.genomics.ReadablePathValidator;
import edu.unc.genomics.io.IntervalFileReader;
import edu.unc.genomics.io.WigFileReader;
import edu.unc.genomics.io.WigFileWriter;
import edu.unc.genomics.io.WigFileException;
import edu.ucsc.genome.TrackHeader;

/**
 * For each interval in Loci file (Bed), output min, mean, max, N of data overlapping from Input file (Wig) to a BedGraph that includes the interval name/id.
 * @author davidcking
 * modified from timothypalpant
 */
public class SplitWigIntervalsToDataFrame extends CommandLineTool {

	private static final Logger log = Logger.getLogger(SplitWigIntervalsToDataFrame.class);

	@Parameter(description = "Input files (Wig)", required = true)
	public List<String> inputWigs = new ArrayList<String>(); 

	@Parameter(names = {"-l", "--loci"}, description = "Loci file (Bed)", 
             required = true, validateWith = ReadablePathValidator.class)
	public Path lociFile;
	@Parameter(names = {"-o", "--output"}, description = "Output R-style dataframe filename")
	public Path outputFile;
    @Parameter(names = { "-s", "--skip-missing" }, description = "Skip missing data during aggregation")
    public boolean skipMissing = false;

    // as per WigMathTool (extends WigAnalysisTool)
    private void prepare() throws CommandLineToolException {
        setup();
    }
    // as per end tool (extends WigMathTool)
    private void setup() throws CommandLineToolException {
        log.debug("Initializing input files");
        for (String inputFile : inputWigs) {
          try {
            addInputFile(WigFileReader.autodetect(Paths.get(inputFile)));
          } catch (IOException e) {
            throw new CommandLineToolException(e);
          }
        }
        log.debug("Initialized " + inputs.size() + " input files");
    }
    private List<WigFileReader> inputs = new ArrayList<WigFileReader>();
    private void addInputFile(WigFileReader wig) {
        inputs.add(wig);
    }
	
	@Override
	public void run() throws IOException {
        log.setLevel((Level)Level.INFO); // defaults to DEBUG
        prepare(); // as per WigAnalysisTool
        log.debug("Initializing output file");
        Charset charset = Charset.forName("US-ASCII");
        BufferedWriter writer = Files.newBufferedWriter(outputFile, charset);
		log.debug("Initializing input file");
		int count = 0, skipped = 0;
		try (
               //WigFileReader wig = WigFileReader.autodetect(inputFile);
               IntervalFileReader<? extends Interval> intervals = IntervalFileReader.autodetect(lociFile) // no ending semicolon in try-with-resources block
            )
        {
            // Output the DataFrame header line
            writer.write("chrom\tchromStart\tchromEnd\tID");
            for (String inputFile : inputWigs) {
                writer.write( String.format("\tmin_%s", inputFile) );
                writer.write( String.format("\tmean_%s", inputFile) );
                writer.write( String.format("\tmax_%s", inputFile) );
                writer.write( String.format("\tN_%s", inputFile) );
            }
            writer.write('\n');

			log.debug("Iterating over all intervals and writing Wig min,mean,max,N for each");
			for (Interval interval : intervals) {
                log.trace("encountered interval: " + interval);
                final String chr = interval.getChr();
                String id = interval.getId();
                if (id == null) {
                    id = "interval_" + count;
                }

                int intervalStart = interval.low();
                int intervalEnd = interval.high();
                String interval_line = String.format("%1$s\t%2$d\t%3$d\t%4$s", chr, intervalStart, intervalEnd, id); // no newline or trailing delim
                writer.write(interval_line);

                int wignum = 0;
                for ( WigFileReader wig : inputs ) {
                    ++wignum;
                    try {
                        Contig query = wig.query(interval);
                        float[] values = query.getValues();
                        float i_max = Float.MIN_VALUE;
                        float i_min = Float.MAX_VALUE;
                        float i_sum = 0;
                        int i_N = 0;
                        for (int startPos = interval.low(); startPos < interval.high(); startPos++) {
                            int i = startPos - interval.low();
                            if (Float.isNaN(values[i]) && skipMissing) continue;
                            i_max = Math.max(i_max, values[i]);
                            i_min = Math.min(i_min, values[i]);
                            i_sum += values[i];
                            i_N++;
                        }
                        float i_mean = i_sum / i_N;

                        String data_line = String.format("\t%1$f\t%2$f\t%3$f\t%4$d", i_min, i_mean, i_max, i_N);
                        writer.write(data_line);
                        
                    } catch (WigFileException e) {
                        log.info("Skipping interval "+interval+" which has no data for input file " + wignum);
                        String error_line  = String.format("\tNA\tNA\tNA\tNA");
                        writer.write(error_line);
                        skipped++;
                    }
                }
            writer.write('\n');
            count++;
            }
		}
        writer.flush();
        writer.close();		
		log.info(count + " intervals processed");
		log.info(skipped + " interval skipped");
	}
	
	public static void main(String[] args) {
		new SplitWigIntervalsToDataFrame().instanceMain(args);
	}

}
