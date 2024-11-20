import model.SimulationOutputEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Writer {
    private static final Logger logger = LoggerFactory.getLogger(Writer.class);

    public static void writeMappingHeader(BufferedWriter writer) {
        try {
            writer.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut\n");
        } catch (Exception e) {
            logger.error("Error while trying to write mapping header", e);
        }
    }

    private final static StringBuilder result = new StringBuilder();

    public static void writeMapping(BufferedWriter writer, SimulationOutputEntry simulationEntry, int readId, String chromosome, String geneId, String transcriptId) {
        try {
            result.setLength(0);

            result.append(readId).append("\t")
                    .append(chromosome).append("\t")
                    .append(geneId).append("\t")
                    .append(transcriptId).append("\t");

            var tFwRegVec = simulationEntry.getTranscriptFwRegionVectors();
            if (tFwRegVec.length != 2) {
                tFwRegVec = new int[2];
                logger.warn("You are a failure");
            }
            result.append(tFwRegVec[0]).append("-").append(tFwRegVec[1]).append("\t");

            // Verarbeitung der tRvRegVec
            var tRvRegVec = simulationEntry.getTranscriptRvRegionVectors();
            if (tRvRegVec.length != 2) {
                tRvRegVec = new int[2];
                logger.warn("You are a failure");
            }
            result.append(tRvRegVec[0]).append("-").append(tRvRegVec[1]).append("\t");

            var genomicFw = simulationEntry.getGenomeFwRegionVectors();
            var firstFwVec = genomicFw.getFirst();
            result.append(firstFwVec[0]).append("-").append(firstFwVec[1]);

            for (int i = 1; i < genomicFw.size(); i++) {
                var vec = genomicFw.get(i);
                result.append("|").append(vec[0]).append("-").append(vec[1]);
            }
            result.append("\t");

            var genomicRw = simulationEntry.getGenomeRvRegionVectors();
            if (genomicRw != null && !genomicRw.isEmpty()) {
                var firstRvVec = genomicRw.getFirst();
                result.append(firstRvVec[0]).append("-").append(firstRvVec[1]);
                for (int i = 1; i < genomicRw.size(); i++) {
                    var vec = genomicRw.get(i);
                    result.append("|").append(vec[0]).append("-").append(vec[1]);
                }
            }
            result.append("\t");

            var fwMutationPos = simulationEntry.getFwMutationIdx();
            if (fwMutationPos != null && !fwMutationPos.isEmpty()) {
                var fwIterator = fwMutationPos.iterator();
                result.append(fwIterator.next());
                while (fwIterator.hasNext()) {
                    result.append(",").append(fwIterator.next());
                }
            }
            result.append("\t");

            var rvMutationPos = simulationEntry.getRvMutationIdx();
            if (rvMutationPos != null && !rvMutationPos.isEmpty()) {
                var rvIterator = rvMutationPos.iterator();
                result.append(rvIterator.next());
                while (rvIterator.hasNext()) {
                    result.append(",").append(rvIterator.next());
                }
            }
            result.append("\n");

            writer.write(result.toString());

        } catch (IOException e) {
            logger.warn("Error while reading line", e);
        }
    }

    private static final Map<Integer, String> qualityCache = new HashMap<>();

    private static final StringBuilder sb = new StringBuilder();
    public static void writeFastq(BufferedWriter writer, int readId, String seq, String qualityString) {
        try {
            sb.setLength(0);
            var valOfReadId = String.valueOf(readId);

            sb.append("@").append(valOfReadId).append("\n")
                    .append(seq).append("\n+").append(valOfReadId).append("\n")
                    .append(qualityString).append("\n");

            writer.write(sb.toString());
        } catch (IOException e) {
            logger.error("Error while trying to write to file", e);
        }
    }

}
