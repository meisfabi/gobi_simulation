import augmentedTree.IntervalTree;
import model.FeatureRecord;
import model.FidxEntry;
import model.Genes;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import utils.Constants;
import utils.GenomeSequenceExtractor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.SplittableRandom;

public class ReadSimulation {
    private static final SplittableRandom random = new SplittableRandom();
    private static final StringBuilder transcriptSeq = new StringBuilder();
    private static final Logger logger = LoggerFactory.getLogger(ReadSimulation.class);

    public static void simulate(Genes gtfData, String fasta, Map<String, FidxEntry> fidxData, int frLength, int sd, int readLength, double mutationRate) throws FileNotFoundException, IOException {
        var seqExtractor = new GenomeSequenceExtractor(new File(fasta), fidxData);
        var geneSeq = new StringBuilder();
        for (var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()) {
            var gene = entry.getValue();
            geneSeq.setLength(0);
            geneSeq.append(seqExtractor.getSequence(gene.getSeqName(), gene.getStart(), gene.getStop()));

            for (var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()) {
                transcriptSeq.setLength(0);
                var exons = transcript.getTranscriptEntry().getPositions().values();
                for (var exon : exons) {
                    if (exon.getStrand() == '+') {
                        transcriptSeq.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                    } else {
                        var revComp = GenomeSequenceExtractor.getReverseComplement(geneSeq.toString(), exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                        transcriptSeq.append(revComp);
                    }
                }
                int fragmentLength;

                fragmentLength = (int)Math.round(Math.max(random.nextGaussian(frLength, sd), readLength));

                if(fragmentLength > transcriptSeq.length()){
                    readLength = transcriptSeq.length();
                    fragmentLength = readLength;
                }

                var fragmentStartPos = random.nextLong(0, transcriptSeq.length() - fragmentLength + 1);
                var fragment = transcriptSeq.substring((int) fragmentStartPos, (int) (fragmentStartPos + fragmentLength));

                var startPos = random.nextInt(0, fragment.length() - readLength + 1);
                var readSeqs = getReadSeq(fragment, startPos, readLength);

                var mutatedSeq = simulateRead(readSeqs[0], mutationRate);
                var mutatedRevSeq = simulateRead(readSeqs[1], mutationRate);

                var fwStart = fragmentStartPos + startPos;
                var fwEnd = fragmentStartPos + startPos + readLength - 1;
                var rwStart = fragmentStartPos + startPos + fragmentLength - readLength;
                var rwEnd = fragmentStartPos + startPos + fragmentLength - 1;

                // map back to reality

                var currentPos = 0;

                for (var exon : exons) {
                    var exonStart = exon.getStart();
                    var exonEnd = exon.getStop();
                    var exonLength = exonEnd - exonStart + 1;
                    if (currentPos + exonLength - 1 >= fwStart && currentPos <= fwEnd) {
                        var startTranscript = Math.max(fwStart, currentPos);
                        var stopTranscript = Math.min(fwEnd, currentPos + exonLength - 1);

                    }
                    currentPos += exonLength;

                }
            }
        }
    }

    private static String[] getReadSeq(String seq, int startPos, int length) {
        return new String[]{
            seq.substring(startPos, startPos + length),
            GenomeSequenceExtractor.getReverseComplement(seq, startPos, startPos + length)
        };
    }

    private static String simulateRead(String seqs, double mutationRate) {
        var bd = new BinomialDistribution(seqs.length(), mutationRate / 100);
        var numMutations = bd.sample();

        var seqArr = seqs.toCharArray();

        var mutatedPositions = new HashSet<Integer>();

        for(int i = 0; i < numMutations; i++){
            int pos;
            do{
                pos = random.nextInt(seqArr.length);
            } while(mutatedPositions.contains(pos));

            mutatedPositions.add(pos);
            var nuc = seqArr[pos];
            char newNuc;

            do{
                newNuc = Constants.NUCLEOTIDES[random.nextInt(Constants.NUCLEOTIDES.length)];
            } while(newNuc == nuc);

            seqArr[pos] = newNuc;
        }

        return new String(seqArr);
    }
}
