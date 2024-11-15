import model.FeatureRecord;
import model.FidxEntry;
import model.Genes;
import model.SimulationOutputEntry;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import utils.Constants;
import utils.GenomeSequenceExtractor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class ReadSimulation {
    private static final SplittableRandom random = new SplittableRandom();
    private static final StringBuilder transcriptSeq = new StringBuilder();
    private static final Logger logger = LoggerFactory.getLogger(ReadSimulation.class);

    public static void simulate(Genes gtfData, String fasta, Map<String, FidxEntry> fidxData, int frLength, int sd, int readLength, double mutationRate, String od) throws FileNotFoundException, IOException {
        var seqExtractor = new GenomeSequenceExtractor(new File(fasta), fidxData);
        var geneSeq = new StringBuilder();
        var output = new ArrayList<SimulationOutputEntry>();
        int readId = 0;
        for (var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()) {
            var gene = entry.getValue();
            var strand = gene.getStrand();

            geneSeq.setLength(0);
            geneSeq.append(seqExtractor.getSequence(gene.getSeqName(), gene.getStart(), gene.getStop()));

            for (var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()) {
                transcriptSeq.setLength(0);
                var exons = transcript.getTranscriptEntry().getPositions().values();
                for (var exon : exons) {
                    transcriptSeq.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                }

                if(strand == '-'){
                    var revComp = GenomeSequenceExtractor.getReverseComplement(transcriptSeq.toString());
                    transcriptSeq.setLength(0);
                    transcriptSeq.append(revComp);
                }

                var fragmentLength = (int) Math.round(Math.max(random.nextGaussian(frLength, sd), readLength));

                if (fragmentLength > transcriptSeq.length()) {
                    readLength = transcriptSeq.length();
                    fragmentLength = readLength;
                }

                var fragmentStartPos = random.nextLong(0, transcriptSeq.length() - fragmentLength + 1);

                /*fragmentStartPos = 0;
                readLength = 75;
                fragmentLength = 90;
*/
                var fragment = transcriptSeq.substring((int) fragmentStartPos, (int) (fragmentStartPos + fragmentLength));

                var fwRead = fragment.substring(0, readLength);
                var rvRead = GenomeSequenceExtractor.getReverseComplement(fragment.substring(fragmentLength - readLength, fragmentLength));
                var mutatedPos = new ArrayList<Integer>();
                var mutatedSeq = simulateRead(fwRead, mutationRate, mutatedPos);
                var mutatedRevPos = new ArrayList<Integer>();
                var mutatedRevSeq = simulateRead(rvRead, mutationRate, mutatedRevPos);

                var fwStart = fragmentStartPos;
                var fwEnd = fragmentStartPos + readLength;
                var fwTranscriptVector = new long[]{fwStart, fwEnd};
                var rvStart = fragmentStartPos + fragmentLength - readLength;
                var rvEnd = fragmentStartPos + fragmentLength;
                var rwTranscriptVector = new long[]{rvStart, rvEnd};
                // map back to reality

                var fwGenomicRegion = getGenomicRegion(exons, fwStart, fwEnd, strand);
                var rvGenomicRegion = getGenomicRegion(exons, rvStart, rvEnd, strand);

                var outputEntry = new SimulationOutputEntry(
                        readId++,
                        mutatedSeq,
                        mutatedRevSeq,
                        gene.getSeqName(),
                        gene.getGeneId(),
                        transcript.getTranscriptId(),
                        fwTranscriptVector,
                        rwTranscriptVector,
                        fwGenomicRegion,
                        rvGenomicRegion,
                        mutatedPos,
                        mutatedRevPos
                );

                output.add(outputEntry);
            }
        }

        Writer.writeMapping(output, od);
        Writer.writeFastqFw(output, od + "/fw.fastq");
        Writer.writeFastqRv(output, od + "/rw.fastq");
    }

    private static String[] getReadSeq(String seq, int length) {
        return new String[]{
                seq.substring(0, length),
                GenomeSequenceExtractor.getReverseComplement(seq, 0,  length)
        };
    }

    private static String simulateRead(String seq, double mutationRate, List<Integer> mutatedPositions) {
        var bd = new BinomialDistribution(seq.length(), mutationRate / 100);
        var numMutations = bd.sample();

        var seqArr = seq.toCharArray();

        for (int i = 0; i < numMutations; i++) {
            int pos;
            do {
                pos = random.nextInt(seqArr.length);
            } while (mutatedPositions.contains(pos));

            mutatedPositions.add(pos);
            var nuc = seqArr[pos];
            char newNuc;

            do {
                newNuc = Constants.NUCLEOTIDES[random.nextInt(Constants.NUCLEOTIDES.length)];
            } while (newNuc == nuc);

            seqArr[pos] = newNuc;
        }

        return new String(seqArr);
    }

    private static List<Long[]> getGenomicRegion(Collection<FeatureRecord> exons, long start, long stop, char strand) {
        var currentPos = 0L;
        var genomicRegions = new ArrayList<Long[]>();

        var exonList = new ArrayList<>(exons);
        if (strand == '-') {
            Collections.reverse(exonList);
        }

        for (var exon : exonList) {
            var exonStart = exon.getStart();
            var exonEnd = exon.getStop();
            var exonLength = exonEnd - exonStart + 1;

            var exonTranscriptStart = currentPos;
            var exonTranscriptEnd = currentPos + exonLength - 1;


            if (exonTranscriptEnd >= start && exonTranscriptStart <= stop) {
                long genomicStart, genomicEnd;

                var overlapStart = Math.max(start, exonTranscriptStart);
                var overlapEnd = Math.min(stop, exonTranscriptEnd + 1);

                if (strand == '+') {
                    genomicStart = exonStart + (overlapStart - exonTranscriptStart);
                    genomicEnd = exonStart + (overlapEnd - exonTranscriptStart);
                } else {
                    genomicStart = exonEnd - (overlapEnd - exonTranscriptStart) + 1;
                    genomicEnd = exonEnd - (overlapStart - exonTranscriptStart) + 1;
                }

                genomicRegions.add(new Long[]{genomicStart, genomicEnd});
            }
            currentPos += exonLength;
        }

        if(strand == '-'){
            Collections.reverse(genomicRegions);
        }

        return genomicRegions;
    }
}
