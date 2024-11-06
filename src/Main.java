import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
public class Main {
    private static Logger logger = LoggerFactory.getLogger(Main.class);
    public static void main(String[] args) {
        ArgumentParser parser = ArgumentParsers.newFor("Main").build().defaultHelp(true).description("Read Simulator");
        parser.addArgument("-length").required(true).help("read length").metavar("<Length>").type(Integer.class);
        parser.addArgument("-frlength").required(true).help("fragment length").metavar("<fragment length>").type(Integer.class);
        parser.addArgument("-SD").required(true).help("Standard Deviation").metavar("<standard deviation>").type(Integer.class);
        parser.addArgument("-readcounts").required(true).help("Path to Read Counts").metavar("<path to read counts>").type(String.class);
        parser.addArgument("-mutationrate").required(true).help("mutation rate between 0.0 and 100.0").metavar("<mutation rate>").type(Double.class);
        parser.addArgument("-fasta").required(true).help("path to fasta").metavar("<fasta>").type(String.class);
        parser.addArgument("-fidx").required(true).help("path to fidx").metavar("<fidx>").type(String.class);
        parser.addArgument("-gtf").required(true).help("path to gtf").metavar("<gtf>").type(String.class);
        parser.addArgument("-od").required(true).help("output directory").metavar("<od>").type(String.class);

        try {
            Namespace res = parser.parseArgs(args);
            var start = System.currentTimeMillis();
            var length = res.get("length");
            var fr_length = res.get("frlength");
            var sd = res.get("SD");
            var read_counts = res.get("readcounts");
            var mutation_rate = res.get("mutationrate");
            var fasta = res.get("fasta");
            var fidx = res.get("fidx");
            var gtf = res.get("gtf");
            var od = res.get("od");
            var data = GtfParser.parse(gtf.toString());
            logger.info(String.format("Time needed for parsing: %s seconds", (System.currentTimeMillis() - start) / 1000.0));


            logger.info(String.format("Time needed for whole program: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
        } catch(ArgumentParserException e){
            parser.printHelp();
        } catch (Exception e) {
            logger.error("Error while executing main", e);
        }
    }
}