package parser;

import model.FidxEntry;
import model.Genes;
import model.ReadCountData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class FidxParser {
    private static int errorLines;
    private static Logger logger = LoggerFactory.getLogger(FidxParser.class);
    private static List<FidxEntry> data;
    public static List<FidxEntry> parse(String inputPath) {
        errorLines = 0;
        logger.info("Starting to parse gtf file");
        data = new ArrayList<>();
        Path path = Path.of(inputPath);

        try (Stream<String> lines = Files.lines(path, StandardCharsets.UTF_8)) {
            lines//.parallel()
                    .filter(line -> !line.trim().isEmpty())
                    .forEach(line -> processLine(line.trim()));
        } catch (Exception e) {
            logger.error("Error while parsing gtf file", e);
        }

        logger.info("GTF-File parsed");
        if (errorLines > 0)
            logger.warn(String.format("%s could not be saved due to an error while parsing", errorLines));

        return data;
    }

    private static void processLine(String line){
        var splitLine = new String[5];
        var currentIdx = 0;
        var currentStart = 0;

        for (int i = 0; i < line.length(); i++) {
            var currentChar = line.charAt(i);
            if (currentChar == '\t') {
                splitLine[currentIdx++] = line.substring(currentStart, i);

                currentStart = i + 1;
            }
        }
        splitLine[currentIdx] = line.substring(currentStart);
        data.add(new FidxEntry(splitLine[0], Long.parseLong(splitLine[2]), Integer.parseInt(splitLine[3]),Integer.parseInt(splitLine[4])));
    }
}