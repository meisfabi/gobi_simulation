package parser;

import model.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import utils.Constants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;

public class GtfParser {
    private static final Logger logger = LoggerFactory.getLogger(GtfParser.class);

    private static int errorLines;

    // Map<Gene_Id, Map<Transcript_Id, TreeMap<StartPosition, FeatureRecord>>[cds or exon]>
    private static Genes parsedGTF;

    public static Genes parse(String inputPath, ReadCountData rcData) {
        errorLines = 0;
        logger.info("Starting to parse gtf file");
        parsedGTF = new Genes();
        Path path = Path.of(inputPath);

        try (Stream<String> lines = Files.lines(path, StandardCharsets.UTF_8)) {
            lines.parallel()
                    .filter(line -> !line.trim().isEmpty() && !line.startsWith("#"))
                    .forEach(line -> processLine(line.trim(), rcData));
        } catch (Exception e) {
            logger.error("Error while parsing gtf file", e);
        }

        logger.info("GTF-File parsed");
        if (errorLines > 0)
            logger.warn(String.format("%s could not be saved due to an error while parsing", errorLines));
        return parsedGTF;
    }

    private static void processLine(String line, ReadCountData rcData) {
        var splitLine = new String[9];
        var currentIdx = 0;
        var currentStart = 0;

        for (int i = 0; i < line.length(); i++) {
            var currentChar = line.charAt(i);
            if (currentChar == '\t') {
                splitLine[currentIdx++] = line.substring(currentStart, i);

                currentStart = i + 1;
                if (currentIdx == 3 && !splitLine[2].equals("exon") && !splitLine[2].equals("transcript") && !splitLine[2].equals("gene")) {
                    return;
                }

                if (currentIdx == 9) {
                    break;
                }
            }
        }

        if (currentIdx < 9) {
            splitLine[currentIdx] = line.substring(currentStart);
        }

        var featureRecord = new FeatureRecord();
        try {
            featureRecord.setStart(Integer.parseInt(splitLine[3]));
            featureRecord.setStop(Integer.parseInt(splitLine[4]));

            if (!splitLine[5].equals(".")) {
                featureRecord.setScore(Double.parseDouble(splitLine[5]));
            } else {
                featureRecord.setScore(-1.0);
            }

            featureRecord.setStrand(splitLine[6].charAt(0));

            if (!splitLine[7].equals(".")) {
                featureRecord.setFrame(Integer.parseInt(splitLine[7]));
            } else {
                featureRecord.setFrame(-1);
            }
            final var stringBuilder = new StringBuilder();
            var attributes = new ArrayList<String>();
            for (int i = 0; i < splitLine[8].length(); i++) {
                var currentChar = splitLine[8].charAt(i);
                if (currentChar == ';') {
                    attributes.add(stringBuilder.toString());
                    stringBuilder.setLength(0);
                } else {
                    stringBuilder.append(currentChar);
                }
            }
            stringBuilder.setLength(0);

            switch (splitLine[2]){
                case "exon":
                    read_exon(rcData, splitLine, attributes, stringBuilder, featureRecord);
                    break;
                case "gene":
                    read_gene(rcData, splitLine, attributes, stringBuilder);
                    break;
                case "transcript":
                    read_transcript(rcData, splitLine, attributes, stringBuilder);
                    break;
            }

        } catch (Exception e) {
            logger.error("Error while trying to parse line", e);
            errorLines++;
        }
    }

    private static void read_gene(ReadCountData rcData, String[] splitLine, ArrayList<String> attributes, StringBuilder stringBuilder){
        String key = null;
        var gene = new Gene();
        for (var attribute : attributes) {
            stringBuilder.setLength(0);
            for (int i = 0; i < attribute.length(); i++) {
                var currentChar = attribute.charAt(i);

                if (currentChar == '\"') continue;
                if (currentChar == ' ') {
                    key = stringBuilder.toString();
                    stringBuilder.setLength(0);
                } else {
                    stringBuilder.append(currentChar);
                }
            }

            var value = stringBuilder.toString();

            if (key == null || key.isBlank()) continue;

            switch (key) {
                case Constants.GENE_ID:
                    if(!rcData.getReadCountData().containsKey(value))
                        return;
                    gene.setGeneId(value);
                    break;
                case Constants.GENE_SOURCE:
                    gene.setGeneSource(value);
                    break;
                case Constants.GENE_BIOTYPE:
                    gene.setGeneBiotype(value);
                    break;
                case Constants.GENE_NAME:
                    gene.setGeneName(value);
                    break;
            }
        }
        var geneId = gene.getGeneId();
        Gene currentGene;


        currentGene = parsedGTF.getFeaturesByTranscriptByGene().get(geneId);
        if (currentGene == null) {
            currentGene = new Gene();
            currentGene.setGeneBiotype(gene.getGeneBiotype());
            currentGene.setGeneName(gene.getGeneName());
            currentGene.setGeneId(geneId);
            currentGene.setGeneSource(gene.getGeneSource());
            currentGene.setSeqName(splitLine[0]);
            currentGene.setStrand(splitLine[6].charAt(0));
        }

        currentGene.setStart(Integer.parseInt(splitLine[3]));
        currentGene.setStop(Integer.parseInt(splitLine[4]));

        var finalGene = currentGene;
        synchronized (parsedGTF){
            parsedGTF.getFeaturesByTranscriptByGene().computeIfAbsent(geneId, id -> finalGene);
        }
    }

    private static void read_exon(ReadCountData rcData, String[] splitLine, ArrayList<String> attributes, StringBuilder stringBuilder, FeatureRecord featureRecord){
        String key = null;
        var gene = new Gene();
        var transcript = new Transcript();

        for (var attribute : attributes) {
            stringBuilder.setLength(0);
            for (int i = 0; i < attribute.length(); i++) {
                var currentChar = attribute.charAt(i);

                if (currentChar == '\"') continue;
                if (currentChar == ' ') {
                    key = stringBuilder.toString();
                    stringBuilder.setLength(0);
                } else {
                    stringBuilder.append(currentChar);
                }
            }

            var value = stringBuilder.toString();

            if (key == null || key.isBlank()) continue;

            switch (key) {
                case Constants.GENE_ID:
                    if(!rcData.getReadCountData().containsKey(value))
                        return;
                    gene.setGeneId(value);
                    break;
                case Constants.TRANSCRIPT_ID:
                    var geneId = gene.getGeneId();
                    if(geneId != null && !rcData.getTranscriptsForGene(geneId).containsKey(value))
                        return;
                    transcript.setTranscriptId(value);
                    break;
                case Constants.EXON_NUMBER:
                    featureRecord.setExonNumber(value);
                    break;
                case Constants.GENE_SOURCE:
                    gene.setGeneSource(value);
                    break;
                case Constants.GENE_BIOTYPE:
                    gene.setGeneBiotype(value);
                    break;
                case Constants.TRANSCRIPT_NAME:
                    transcript.setTranscriptName(value);
                    break;
                case Constants.TRANSCRIPT_SOURCE:
                    transcript.setTranscriptSource(value);
                    break;
                case Constants.TAG:
                    featureRecord.setTag(value);
                    break;
                case Constants.CCDS_ID:
                    featureRecord.setCcdsId(value);
                    break;
                case Constants.PROTEIN_ID:
                    featureRecord.setProteinId(value);
                    break;
                case Constants.GENE_NAME:
                    gene.setGeneName(value);
                    break;
            }
        }
        var geneId = gene.getGeneId();
        Gene currentGene;


        currentGene = parsedGTF.getFeaturesByTranscriptByGene().get(geneId);

        if (currentGene == null) {
            currentGene = new Gene();
            currentGene.setGeneBiotype(gene.getGeneBiotype());
            currentGene.setGeneName(gene.getGeneName());
            currentGene.setGeneId(geneId);
            currentGene.setGeneSource(gene.getGeneSource());
            currentGene.setSeqName(splitLine[0]);
            currentGene.setStrand(splitLine[6].charAt(0));
            parsedGTF.getFeaturesByTranscriptByGene().put(geneId, currentGene);
        }

        var transcriptId = transcript.getTranscriptId();
        if (geneId.isEmpty() || transcriptId.isEmpty()) {
            logger.warn("Could not add GtfRecord because geneId or transcriptId was empty");
            return;
        }

        var transcriptMap = currentGene.getTranscriptMap();
        if (transcriptMap == null) {
            transcriptMap = new ConcurrentHashMap<>();
        }

        var transcriptEntry = transcriptMap
                .computeIfAbsent(transcriptId, id -> transcript)
                .getTranscriptEntry();

        transcriptEntry.addRecord(featureRecord);
    }

    private static void read_transcript(ReadCountData rcData, String[] splitLine, ArrayList<String> attributes, StringBuilder stringBuilder){
        String key = null;
        var gene = new Gene();
        var transcript = new Transcript();

        for (var attribute : attributes) {
            stringBuilder.setLength(0);
            for (int i = 0; i < attribute.length(); i++) {
                var currentChar = attribute.charAt(i);

                if (currentChar == '\"') continue;
                if (currentChar == ' ') {
                    key = stringBuilder.toString();
                    stringBuilder.setLength(0);
                } else {
                    stringBuilder.append(currentChar);
                }
            }

            var value = stringBuilder.toString();

            if (key == null || key.isBlank()) continue;

            switch (key) {
                case Constants.GENE_ID:
                    if(!rcData.getReadCountData().containsKey(value))
                        return;
                    gene.setGeneId(value);
                    break;
                case Constants.TRANSCRIPT_ID:
                    var geneId = gene.getGeneId();
                    if(geneId != null && !rcData.getTranscriptsForGene(geneId).containsKey(value))
                        return;
                    transcript.setTranscriptId(value);
                    break;
                case Constants.GENE_SOURCE:
                    gene.setGeneSource(value);
                    break;
                case Constants.GENE_BIOTYPE:
                    gene.setGeneBiotype(value);
                    break;
                case Constants.TRANSCRIPT_NAME:
                    transcript.setTranscriptName(value);
                    break;
                case Constants.TRANSCRIPT_SOURCE:
                    transcript.setTranscriptSource(value);
                    break;
                case Constants.GENE_NAME:
                    gene.setGeneName(value);
                    break;
            }
        }
        var geneId = gene.getGeneId();
        Gene currentGene;


        currentGene = parsedGTF.getFeaturesByTranscriptByGene().get(geneId);

        if (currentGene == null) {
            currentGene = new Gene();
            currentGene.setGeneBiotype(gene.getGeneBiotype());
            currentGene.setGeneName(gene.getGeneName());
            currentGene.setGeneId(geneId);
            currentGene.setGeneSource(gene.getGeneSource());
            currentGene.setSeqName(splitLine[0]);
            currentGene.setStrand(splitLine[6].charAt(0));
            parsedGTF.getFeaturesByTranscriptByGene().put(geneId, currentGene);
        }

        var transcriptId = transcript.getTranscriptId();
        if (geneId.isEmpty() || transcriptId.isEmpty()) {
            logger.warn("Could not add GtfRecord because geneId or transcriptId was empty");
            return;
        }

        transcript.setStart(Integer.parseInt(splitLine[3]));
        transcript.setStop(Integer.parseInt(splitLine[4]));

        var transcriptMap = currentGene.getTranscriptMap();
        if (transcriptMap == null) {
            transcriptMap = new ConcurrentHashMap<>();
        }

        transcript.setReadCount(rcData.getReadCountData().get(geneId).get(transcriptId).getReadCounts());

        transcriptMap
                .computeIfAbsent(transcriptId, id -> transcript)
                .getTranscriptEntry();    }

}

