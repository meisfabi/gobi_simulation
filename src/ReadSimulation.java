import augmentedTree.IntervalTree;
import model.FeatureRecord;
import model.FidxEntry;
import model.Genes;
import org.apache.commons.math3.distribution.BinomialDistribution;
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
    private static final StringBuilder stringBuilder = new StringBuilder();

    public static void simulate(Genes gtfData, String fasta, Map<String, FidxEntry> fidxData, int frLength, int sd, int length, double mutationRate) throws FileNotFoundException, IOException {
        var seqExtractor = new GenomeSequenceExtractor(new File(fasta), fidxData);
        var geneSeq = new StringBuilder();
        for (var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()) {
            var gene = entry.getValue();
            geneSeq.setLength(0);
            geneSeq.append(seqExtractor.getSequence(gene.getSeqName(), gene.getStart(), gene.getStop()));

            for (var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()) {
                stringBuilder.setLength(0);
                var exons = transcript.getTranscriptEntry().getPositions().values();
                var exonTree = new IntervalTree<FeatureRecord>();
                for (var exon : exons) {
                    exonTree.add(exon);
                    if (exon.getStrand() == '+') {
                        stringBuilder.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                    } else {
                        var revComp = GenomeSequenceExtractor.getReverseComplement(geneSeq.toString(), exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                        stringBuilder.append(revComp);
                    }
                }

                var fragmentLength = Math.round(Math.max(random.nextGaussian(frLength, sd), length));
                var position = random.nextLong(0, stringBuilder.length() - fragmentLength);
                var fragment = stringBuilder.substring((int) position, (int) (position + fragmentLength));

                var startPos = random.nextInt(0, 75 - length + 1);
                var readSeqs = getReadSeq(fragment, startPos, length);

                var mutatedSeq = simulateRead(readSeqs[0], mutationRate);
                var mutatedRevSeq = simulateRead(readSeqs[1], mutationRate);

                var a = "";
            }
        }
    }

    private static String[] getReadSeq(String seq, int startPos, int length) {
        var seqLength = seq.length();

        var fragment = seq.substring(startPos, length);

        return new String[]{
            fragment,
            GenomeSequenceExtractor.getReverseComplement(fragment)
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
