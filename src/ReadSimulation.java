import extensions.ExecutorServiceExtensions;
import model.FeatureRecord;
import model.FidxEntry;
import model.Genes;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import pooling.SimulationOutputFactory;
import utils.Constants;
import utils.GenomeSequenceExtractor;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class ReadSimulation {
    private static final StringBuilder transcriptSeq = new StringBuilder();
    private static final Logger logger = LoggerFactory.getLogger(ReadSimulation.class);
    private static final Object writeLock = new Object();
    private static final SplittableRandom random = new SplittableRandom();
    private static final ConcurrentHashMap<Long, SplittableRandom> randomPool = new ConcurrentHashMap<>();
    public static void simulate(Genes gtfData, String fasta, Map<String, FidxEntry> fidxData, int frLength, int sd, int readLength, double mutationRate, String od) throws FileNotFoundException, IOException {
        final var seqExtractor = new GenomeSequenceExtractor(new File(fasta), fidxData);
        final var geneSeq = new StringBuilder();
        final var simulationOutputFactory = new SimulationOutputFactory();
        final var readId = new AtomicInteger(0);
        final var qualityString = "I".repeat(readLength);
        lambda = Math.exp(-(readLength * (mutationRate / 100)));
        var numberOfThreads = Runtime.getRuntime().availableProcessors();
        var executor = Executors.newFixedThreadPool(numberOfThreads);
        var minZ = (readLength - frLength) / (double)sd;
        try (BufferedWriter mappingWriter = new BufferedWriter(new FileWriter(od + "/read.mappinginfo")); BufferedWriter fwFastqWriter = new BufferedWriter(new FileWriter(od + "/fw.fastq")); BufferedWriter rvFastqWriter = new BufferedWriter(new FileWriter(od + "/rw.fastq"))) {

            Writer.writeMappingHeader(mappingWriter);



            for (var entry : gtfData.getFeaturesByTranscriptByGene().entrySet()) {
                var geneId = entry.getKey();
                var gene = entry.getValue();
                var strand = gene.getStrand();
                var chromosome = gene.getSeqName();
                geneSeq.setLength(0);
                geneSeq.append(seqExtractor.getSequence(chromosome, gene.getStart(), gene.getStop()));

                for (var transcript : gene.getTranscriptMapArray()[Constants.EXON_INDEX].values()) {
                    var transcriptId = transcript.getTranscriptId();
                    var readCount = transcript.getReadCount();
                    transcriptSeq.setLength(0);
                    var exons = transcript.getTranscriptEntry().getPositions();
                    for (var exon : exons) {
                        transcriptSeq.append(geneSeq, exon.getStart() - gene.getStart(), exon.getStop() - gene.getStart() + 1);
                    }

                    var forwardSeq = transcriptSeq.toString();
                    String otherSeq;

                    if (strand == '-') {
                        otherSeq = forwardSeq;
                        transcriptSeq.setLength(0);
                        transcriptSeq.append(GenomeSequenceExtractor.getReverseComplement(forwardSeq));
                    } else {
                        otherSeq = GenomeSequenceExtractor.getReverseComplement(forwardSeq);
                    }

                    var transcriptLength = transcriptSeq.length();
                    var maxZ = (transcriptLength - frLength) / (double)sd;
                    var tasks = new ArrayList<Callable<Void>>();
                    for (int i = 0; i < readCount; i++) {
                        tasks.add(() -> {
                            var z = random.nextGaussian(); // Standard-Normalverteilung

                            var scaledZ = (z - minZ) / (maxZ - minZ);
                            scaledZ = Math.max(0, Math.min(1, scaledZ)); // Begrenzen auf [0, 1]

                            var fragmentLength = (int) (readLength + scaledZ * (transcriptLength - readLength));

                            int fragmentStartPos;
                            fragmentStartPos = random.nextInt(0, transcriptLength - fragmentLength + 1);

                            var fwRead = transcriptSeq.substring(fragmentStartPos, fragmentStartPos + readLength);
                            String rvRead;

                            if (strand == '-') {
                                rvRead = forwardSeq.substring(transcriptLength - fragmentStartPos - fragmentLength, transcriptLength - fragmentStartPos - fragmentLength + readLength);
                            } else {
                                rvRead = otherSeq.substring(transcriptLength - fragmentStartPos - fragmentLength, transcriptLength - fragmentStartPos - fragmentLength + readLength);
                            }

                            var mutatedPos = new TreeSet<Integer>();
                            var mutatedSeq = simulateRead(fwRead, mutatedPos);
                            var mutatedRevPos = new TreeSet<Integer>();
                            var mutatedRevSeq = simulateRead(rvRead, mutatedRevPos);

                            var fwStart = fragmentStartPos;
                            var fwEnd = fragmentStartPos + readLength;
                            var fwTranscriptVector = new int[]{fwStart, fwEnd};
                            var rvStart = fragmentStartPos + fragmentLength - readLength;
                            var rvEnd = fragmentStartPos + fragmentLength;
                            var rwTranscriptVector = new int[]{rvStart, rvEnd};

                            // Map back to reality
                            var genomicRegions = getGenomicRegions(exons, fwStart, fwEnd, rvStart, rvEnd, strand);
                            var fwGenomicRegion = genomicRegions[0];
                            var rvGenomicRegion = genomicRegions[1];

                            var outputEntry = simulationOutputFactory.createObject()
                                    .withReadId(readId.getAndIncrement())
                                    .withTranscriptFwRegionVectors(fwTranscriptVector)
                                    .withTranscriptRvRegionVectors(rwTranscriptVector)
                                    .withGenomeFwRegionVectors(fwGenomicRegion)
                                    .withGenomeRvRegionVectors(rvGenomicRegion)
                                    .withFwMutationIdx(mutatedPos)
                                    .withRvMutationIdx(mutatedRevPos)
                                    .build();

                            var currentReadId = outputEntry.getReadId();
                            synchronized (writeLock) {
                                Writer.writeMapping(mappingWriter, outputEntry, currentReadId, chromosome, geneId, transcriptId);
                                Writer.writeFastq(fwFastqWriter, currentReadId, mutatedSeq, qualityString);
                                Writer.writeFastq(rvFastqWriter, currentReadId, mutatedRevSeq, qualityString);
                            }

                            return null;
                        });
                    }

                    try {
                        executor.invokeAll(tasks);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                        logger.error("error", e);
                    }
                }
            }
        } catch (Exception e) {
            logger.error("Error while simulating reads", e);
        } finally {
            ExecutorServiceExtensions.shutdownExecutorService(executor);
        }
    }

    private static double lambda;

    public static int samplePoisson() {
        var k = 0;
        var p = 1.0;

        do {
            k++;
            p *= random.nextDouble();
        } while (p > lambda);

        return k - 1;
    }

    private static String simulateRead(String seq, Set<Integer> mutatedPositions) {

        var numMutations = samplePoisson();
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

    private static List<int[]>[] getGenomicRegions(
            TreeSet<FeatureRecord> exons,
            int fwStart, int fwEnd,
            int rvStart, int rvEnd,
            char strand) {

        var currentPos = 0;
        var fwGenomicRegions = new ArrayList<int[]>();
        var rvGenomicRegions = new ArrayList<int[]>();

        // Reverse iteration if needed
        Set<FeatureRecord> iterableExons = exons;
        if (strand == '-') {
            iterableExons = exons.descendingSet();
        }

        for (var exon : iterableExons) {
            var exonStart = exon.getStart();
            var exonEnd = exon.getStop();
            var exonLength = exonEnd - exonStart + 1;
            var exonTranscriptStart = currentPos;
            var exonTranscriptEnd = currentPos + exonLength - 1;
            var posExceededFwStart = false;

            if (exonTranscriptStart >= Math.max(fwEnd, rvEnd)) {
                break;
            }

            if (exonTranscriptStart >= fwEnd) {
                posExceededFwStart = true;
            }

            if (exonTranscriptEnd >= fwStart && !posExceededFwStart) {
                var overlapStart = Math.max(fwStart, exonTranscriptStart);
                var overlapEnd = Math.min(fwEnd, exonTranscriptEnd + 1);

                int genomicStart, genomicEnd;
                if (strand == '+') {
                    genomicStart = exonStart + (overlapStart - exonTranscriptStart);
                    genomicEnd = exonStart + (overlapEnd - exonTranscriptStart);
                } else {
                    genomicStart = exonEnd - (overlapEnd - exonTranscriptStart) + 1;
                    genomicEnd = exonEnd - (overlapStart - exonTranscriptStart) + 1;
                }

                fwGenomicRegions.add(new int[]{genomicStart, genomicEnd});
            }

            if (exonTranscriptStart >= rvEnd) {
                currentPos += exonLength;
                continue;
            }

            if (exonTranscriptEnd >= rvStart) {
                var overlapStart = Math.max(rvStart, exonTranscriptStart);
                var overlapEnd = Math.min(rvEnd, exonTranscriptEnd + 1);

                int genomicStart, genomicEnd;
                if (strand == '+') {
                    genomicStart = exonStart + (overlapStart - exonTranscriptStart);
                    genomicEnd = exonStart + (overlapEnd - exonTranscriptStart);
                } else {
                    genomicStart = exonEnd - (overlapEnd - exonTranscriptStart) + 1;
                    genomicEnd = exonEnd - (overlapStart - exonTranscriptStart) + 1;
                }

                rvGenomicRegions.add(new int[]{genomicStart, genomicEnd});
            }

            currentPos += exonLength;
        }

        // No need to reverse if already handled by descendingSet
        return new List[]{fwGenomicRegions, rvGenomicRegions};
    }
}
