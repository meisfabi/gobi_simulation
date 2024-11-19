import model.SimulationOutputEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
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

    public static void writeMapping(BufferedWriter writer, SimulationOutputEntry simulationEntry) {
        try {
            writer.write(String.format("%s\t%s\t%s\t%s\t",
                    simulationEntry.getReadId(),
                    simulationEntry.getChromosome(),
                    simulationEntry.getGeneId(),
                    simulationEntry.getTranscriptId()));

            var tFwRegVec = simulationEntry.getTranscriptFwRegionVectors();
            if (tFwRegVec.length != 2) {
                tFwRegVec = new long[2];
                logger.warn("You are a failure");
            }

            writer.write(String.format("%s-%s\t", tFwRegVec[0], tFwRegVec[1]));

            var tRvRegVec = simulationEntry.getTranscriptRvRegionVectors();
            if (tRvRegVec.length != 2) {
                tRvRegVec = new long[2];
                logger.warn("You are a failure");
            }
            writer.write(String.format("%s-%s\t", tRvRegVec[0], tRvRegVec[1]));
            var genomicFw = simulationEntry.getGenomeFwRegionVectors();

            var firstFwVec = genomicFw.getFirst();


            writer.write(String.format("%s-%s", firstFwVec[0], firstFwVec[1]));
            for (int i = 1; i < genomicFw.size(); i++) {
                var vec = genomicFw.get(i);
                writer.write(String.format("|%s-%s", vec[0], vec[1]));
            }
            writer.write("\t");

            var genomicRw = simulationEntry.getGenomeRvRegionVectors();
            if (!(genomicRw == null || genomicRw.isEmpty())) {
                var firstRvVec = genomicRw.getFirst();

                writer.write(String.format("%s-%s", firstRvVec[0], firstRvVec[1]));
                for (int i = 1; i < genomicRw.size(); i++) {
                    var vec = genomicRw.get(i);
                    writer.write(String.format("|%s-%s", vec[0], vec[1]));
                }
            }


            writer.write("\t");

            var fwMutationPos = simulationEntry.getFwMutationIdx();
            if (!fwMutationPos.isEmpty()) {
                writer.write(String.valueOf(fwMutationPos.getFirst()));
                for (int i = 1; i < fwMutationPos.size(); i++) {
                    writer.write(",");
                    writer.write(String.valueOf(fwMutationPos.get(i)));
                    ;
                }
            }

            writer.write("\t");

            var rvMutationPos = simulationEntry.getRvMutationIdx();
            if (rvMutationPos != null && !rvMutationPos.isEmpty()) {
                writer.write(String.valueOf(rvMutationPos.getFirst()));
                for (int i = 1; i < rvMutationPos.size(); i++) {
                    writer.write(",");
                    writer.write(String.valueOf(rvMutationPos.get(i)));
                    ;
                }
            }
            writer.write("\n");
        } catch (IOException e) {
            logger.warn("Error while reading line", e);
        }
    }

    private static final Map<Integer, String> qualityCache = new HashMap<>();

    private static String getQualityString(int length) {
        return qualityCache.computeIfAbsent(length, "I"::repeat);
    }

    public static void writeFastq(BufferedWriter writer, int readId, String seq) {
        try {
            var valOfReadId = String.valueOf(readId);
            writer.write("@" + valOfReadId + "\n" + seq + "\n+" + valOfReadId + "\n");
            writer.write(getQualityString(seq.length()) + "\n");
        } catch (IOException e) {
            logger.error("Error while trying to write to file", e);
        }
    }

}
