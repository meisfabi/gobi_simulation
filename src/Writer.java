import model.SimulationOutputEntry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Collection;

public class Writer {
    private static final Logger logger = LoggerFactory.getLogger(Writer.class);
    public static void writeMapping(Collection<SimulationOutputEntry> simulationEntries, String outputDir){
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir + "/read.mappinginfo", false))){
            var startTime = System.currentTimeMillis();
            writer.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut\n");
            for(var simulationEntry : simulationEntries) {
                writer.write(String.format("%s\t%s\t%s\t%s\t",
                        simulationEntry.getReadId(),
                        simulationEntry.getChromosome(),
                        simulationEntry.getGeneId(),
                        simulationEntry.getTranscriptId()));

                var tFwRegVec = simulationEntry.getTranscriptFwRegionVectors();
                if(tFwRegVec.length != 2){
                    tFwRegVec = new long[2];
                    logger.warn("You are a failure");
                }

                writer.write(String.format("%s-%s\t", tFwRegVec[0], tFwRegVec[1]));

                var tRvRegVec = simulationEntry.getTranscriptRvRegionVectors();
                if(tRvRegVec.length != 2){
                    tRvRegVec = new long[2];
                    logger.warn("You are a failure");
                }
                writer.write(String.format("%s-%s\t", tRvRegVec[0], tRvRegVec[1]));
                var genomicFw = simulationEntry.getGenomeFwRegionVectors();

                var firstFwVec = genomicFw.getFirst();


                writer.write(String.format("%s-%s", firstFwVec[0], firstFwVec[1]));
                for(int i = 1; i < genomicFw.size(); i++){
                    var vec = genomicFw.get(i);
                    writer.write(String.format("|%s-%s",vec[0], vec[1]));
                }
                writer.write("\t");

                var genomicRw = simulationEntry.getGenomeRvRegionVectors();
                if(!(genomicRw == null || genomicRw.isEmpty())){
                    var firstRvVec = genomicRw.getFirst();

                    writer.write(String.format("%s-%s", firstRvVec[0], firstRvVec[1]));
                    for(int i = 1; i < genomicRw.size(); i++){
                        var vec = genomicRw.get(i);
                        writer.write(String.format("|%s-%s",vec[0], vec[1]));
                    }
                }


                writer.write("\t");

                var fwMutationPos = simulationEntry.getFwMutationIdx();
                if(!fwMutationPos.isEmpty()){
                    writer.write(String.valueOf(fwMutationPos.getFirst()));
                    for(int i = 1; i < fwMutationPos.size(); i++){
                        writer.write(",");
                        writer.write(String.valueOf(fwMutationPos.get(i)));;
                    }
                }

                writer.write("\t");

                var rvMutationPos = simulationEntry.getRvMutationIdx();
                if(!rvMutationPos.isEmpty()){
                    writer.write(String.valueOf(rvMutationPos.getFirst()));
                    for(int i = 1; i < rvMutationPos.size(); i++){
                        writer.write(",");
                        writer.write(String.valueOf(rvMutationPos.get(i)));;
                    }
                }
                writer.write("\n");
            }
            logger.info(String.format("Time for writing to file: %s", (System.currentTimeMillis() - startTime) / 1000.0));

        } catch (Exception e){
            logger.error("Error while trying to write to file", e);
        }
    }

    public static void writeFastqFw(Collection<SimulationOutputEntry> simulationEntries, String output){
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(output, false))){
            var startTime = System.currentTimeMillis();
            for(var simulationEntry : simulationEntries) {
                writer.write("@");
                writer.write(String.valueOf(simulationEntry.getReadId()));
                writer.write("\n");
                writer.write(simulationEntry.getMutatedFwSeq());
                writer.write("\n");
                writer.write("+");
                writer.write(String.valueOf(simulationEntry.getReadId()));
                writer.write("\n");
                writer.write("I".repeat(simulationEntry.getMutatedFwSeq().length()));
                writer.write("\n");
            }
            logger.info(String.format("Time for writing to file: %s", (System.currentTimeMillis() - startTime) / 1000.0));

        } catch (Exception e){
            logger.error("Error while trying to write to file", e);
        }
    }

    public static void writeFastqRv(Collection<SimulationOutputEntry> simulationEntries, String output){
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(output, false))){
            var startTime = System.currentTimeMillis();
            for(var simulationEntry : simulationEntries) {
                writer.write("@");
                writer.write(String.valueOf(simulationEntry.getReadId()));
                writer.write("\n");
                writer.write(simulationEntry.getMutatedRvSeq());
                writer.write("\n");
                writer.write("+");
                writer.write(String.valueOf(simulationEntry.getReadId()));
                writer.write("\n");
                writer.write("I".repeat(simulationEntry.getMutatedRvSeq().length()));
                writer.write("\n");
            }
            logger.info(String.format("Time for writing to file: %s", (System.currentTimeMillis() - startTime) / 1000.0));

        } catch (Exception e){
            logger.error("Error while trying to write to file", e);
        }
    }
}
