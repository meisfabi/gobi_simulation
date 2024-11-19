import model.FeatureRecord;
import model.FidxEntry;
import model.Genes;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import pooling.SimulationOutputFactory;
import utils.Constants;
import utils.GenomeSequenceExtractor;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class ReadSimulation {
    private static final SplittableRandom random = new SplittableRandom();
    private static final StringBuilder transcriptSeq = new StringBuilder();
    private static final Logger logger = LoggerFactory.getLogger(ReadSimulation.class);
    private static final Map<Integer, BinomialDistribution> binomialDistributionCache = new HashMap<>();

    public static void simulate(Genes gtfData, String fasta, Map<String, FidxEntry> fidxData, int frLength, int sd, int readLength, double mutationRate, String od) throws FileNotFoundException, IOException {
        var seqExtractor = new GenomeSequenceExtractor(new File(fasta), fidxData);
        var geneSeq = new StringBuilder();
        var simulationOutputFactory = new SimulationOutputFactory();
        var readId = new AtomicInteger(0);

        try (BufferedWriter mappingWriter = new BufferedWriter(new FileWriter(od + "/read.mappinginfo")); BufferedWriter fwFastqWriter = new BufferedWriter(new FileWriter(od + "/fw.fastq")); BufferedWriter rvFastqWriter = new BufferedWriter(new FileWriter(od + "/rw.fastq"))) {

            Writer.writeMappingHeader(mappingWriter);

            for (var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()) {
                var gene = entry.getValue();
                var strand = gene.getStrand();

                geneSeq.setLength(0);
                geneSeq.append(seqExtractor.getSequence(gene.getSeqName(), gene.getStart(), gene.getStop()));

                for (var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()) {
                    var readCount = transcript.getReadCount();
                    transcriptSeq.setLength(0);
                    var exons = transcript.getTranscriptEntry().getPositions().values();
                    for (var exon : exons) {
                        transcriptSeq.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                    }

                    if (strand == '-') {
                        var revComp = GenomeSequenceExtractor.getReverseComplement(transcriptSeq.toString());
                        transcriptSeq.setLength(0);
                        transcriptSeq.append(revComp);
                    }

                    for (int i = 0; i < readCount; i++) {
                        int fragmentLength;

                        do{
                            fragmentLength = (int) Math.round(random.nextGaussian(frLength, sd));
                        } while (fragmentLength < readLength || transcriptSeq.length() < fragmentLength);

                        var binomialDistribution = binomialDistributionCache.computeIfAbsent(readLength, rl ->
                                new BinomialDistribution(rl, mutationRate / 100)
                        );

                        var fragmentStartPos = random.nextLong(0, transcriptSeq.length() - fragmentLength + 1);
                        var fragment = transcriptSeq.substring((int) fragmentStartPos, (int) (fragmentStartPos + fragmentLength));

                        var fwRead = fragment.substring(0, readLength);
                        var rvRead = GenomeSequenceExtractor.getReverseComplement(fragment.substring(fragmentLength - readLength, fragmentLength));
                        var mutatedPos = new ArrayList<Integer>();
                        var mutatedSeq = simulateRead(fwRead, binomialDistribution, mutatedPos);
                        var mutatedRevPos = new ArrayList<Integer>();
                        var mutatedRevSeq = simulateRead(rvRead, binomialDistribution, mutatedRevPos);

                        var fwStart = fragmentStartPos;
                        var fwEnd = fragmentStartPos + readLength;
                        var fwTranscriptVector = new long[]{fwStart, fwEnd};
                        var rvStart = fragmentStartPos + fragmentLength - readLength;
                        var rvEnd = fragmentStartPos + fragmentLength;
                        var rwTranscriptVector = new long[]{rvStart, rvEnd};

                        // map back to reality

                        var genomicRegions = getGenomicRegions(exons, fwStart, fwEnd, rvStart, rvEnd, strand);
                        var fwGenomicRegion = genomicRegions[0];
                        var rvGenomicRegion = genomicRegions[1];

                        var outputEntry = simulationOutputFactory.createObject()
                                .withReadId(readId.getAndIncrement())
                                .withChromosome(gene.getSeqName())
                                .withGeneId(gene.getGeneId())
                                .withTranscriptId(transcript.getTranscriptId())
                                .withTranscriptFwRegionVectors(fwTranscriptVector)
                                .withTranscriptRvRegionVectors(rwTranscriptVector)
                                .withGenomeFwRegionVectors(fwGenomicRegion)
                                .withGenomeRvRegionVectors(rvGenomicRegion)
                                .withFwMutationIdx(mutatedPos)
                                .withRvMutationIdx(mutatedRevPos)
                                .build();

                        Writer.writeMapping(mappingWriter, outputEntry);
                        Writer.writeFastq(fwFastqWriter, outputEntry.getReadId(), mutatedSeq);
                        Writer.writeFastq(rvFastqWriter, outputEntry.getReadId(), mutatedRevSeq);
                    }
                }
            }
        } catch (Exception e) {
            logger.error("Error while simulating reads", e);
        }
    }


    private static String simulateRead(String seq, BinomialDistribution bd, List<Integer> mutatedPositions) {

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

    private static List<Long[]>[] getGenomicRegions(
            Collection<FeatureRecord> exons,
            long fwStart, long fwEnd,
            long rvStart, long rvEnd,
            char strand) {

        var currentPos = 0L;
        var fwGenomicRegions = new ArrayList<Long[]>();
        var rvGenomicRegions = new ArrayList<Long[]>();

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

            if (exonTranscriptStart >= Math.max(fwEnd, rvEnd)) {
                break;
            }

            // Verarbeitung der fwGenomicRegion
            if (exonTranscriptEnd >= fwStart && exonTranscriptStart < fwEnd) {
                var overlapStart = Math.max(fwStart, exonTranscriptStart);
                var overlapEnd = Math.min(fwEnd, exonTranscriptEnd + 1);

                long genomicStart, genomicEnd;
                if (strand == '+') {
                    genomicStart = exonStart + (overlapStart - exonTranscriptStart);
                    genomicEnd = exonStart + (overlapEnd - exonTranscriptStart);
                } else {
                    genomicStart = exonEnd - (overlapEnd - exonTranscriptStart) + 1;
                    genomicEnd = exonEnd - (overlapStart - exonTranscriptStart) + 1;
                }

                fwGenomicRegions.add(new Long[]{genomicStart, genomicEnd});
            }

            // Verarbeitung der rvGenomicRegion
            if (exonTranscriptEnd >= rvStart && exonTranscriptStart < rvEnd) {
                var overlapStart = Math.max(rvStart, exonTranscriptStart);
                var overlapEnd = Math.min(rvEnd, exonTranscriptEnd + 1);

                long genomicStart, genomicEnd;
                if (strand == '+') {
                    genomicStart = exonStart + (overlapStart - exonTranscriptStart);
                    genomicEnd = exonStart + (overlapEnd - exonTranscriptStart);
                } else {
                    genomicStart = exonEnd - (overlapEnd - exonTranscriptStart) + 1;
                    genomicEnd = exonEnd - (overlapStart - exonTranscriptStart) + 1;
                }

                rvGenomicRegions.add(new Long[]{genomicStart, genomicEnd});
            }

            currentPos += exonLength;
        }

        if (strand == '-') {
            Collections.reverse(fwGenomicRegions);
            Collections.reverse(rvGenomicRegions);
        }

        return new List[]{fwGenomicRegions, rvGenomicRegions};
    }
}
