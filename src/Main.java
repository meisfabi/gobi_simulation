import model.Genes;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import parser.FidxParser;
import parser.GtfParser;
import parser.ReadCountParser;
import utils.Constants;

import java.io.File;
import java.util.Random;

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
            int length = res.get("length");
            int frLength = res.get("frlength");
            int sd = res.get("SD");
            var read_counts = res.get("readcounts");
            var mutationRate = res.get("mutationrate");
            var fasta = res.get("fasta");
            var fidx = res.get("fidx");
            var gtf = res.get("gtf");
            var od = res.get("od");
            var start = System.currentTimeMillis();
            var rcData = ReadCountParser.parse(read_counts.toString());
            var fidxData = FidxParser.parse(fidx.toString());
            var gtfData = GtfParser.parse(gtf.toString(), rcData);
            logger.info(String.format("Time needed for parsing: %s seconds", (System.currentTimeMillis() - start) / 1000.0));

            var seqExtractor = new GenomeSequenceExtractor(new File(fasta.toString()), fidxData);
            var geneSeq = new StringBuilder();
            for(var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()){
                var gene = entry.getValue();
                geneSeq.setLength(0);
                geneSeq.append(seqExtractor.getSequence(gene.getSeqName(), gene.getStart(), gene.getStop()));

                for(var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()){
                    var transcriptSeq = new StringBuilder();
                    for(var exon : transcript.getTranscriptEntry().getPositions().values()){
                        if(exon.getStrand() == '+'){
                            transcriptSeq.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                        } else{
                            var revComp = seqExtractor.getReverseComplement(geneSeq.toString(), exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                            transcriptSeq.append(revComp);
                        }
                    }

                    var random = new Random();
                    var fragmentLength = Math.round(Math.max(random.nextGaussian(frLength, sd), length));
                    var position = random.nextLong(0, transcriptSeq.length() - fragmentLength);
                    var fragment = transcriptSeq.substring((int)position, (int)(position + fragmentLength));
                    var a = "";
                }
            }

            logger.info(String.format("Time needed for whole program: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
        } catch(ArgumentParserException e){
            parser.printHelp();
        } catch (Exception e) {
            logger.error("Error while executing main", e);
        }
    }
}