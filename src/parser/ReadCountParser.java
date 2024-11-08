package parser;

import model.ReadCountData;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.stream.Stream;

public class ReadCountParser {
    private static Logger logger = LoggerFactory.getLogger(ReadCountParser.class);
    private static int errorLines;
    private static ReadCountData data;
    public static ReadCountData parse(String inputPath) {
        Path path = Path.of(inputPath);
        errorLines = 0;
        logger.info("Starting to parse rc file");
        data = new ReadCountData();

        try (Stream<String> lines = Files.lines(path, StandardCharsets.UTF_8)) {
            lines.parallel()
                    .filter(line -> !line.trim().isEmpty() && !line.startsWith("gene"))
                    .forEach(line -> processLine(line.trim()));
        } catch (Exception e) {
            logger.error("Error while parsing rc file", e);
        }

        logger.info("RC-File parsed");
        if (errorLines > 0)
            logger.warn(String.format("%s could not be saved due to an error while parsing", errorLines));

        return data;
    }

    private static void processLine(String line){
        var splitLine = new String[3];
        var currentIdx = 0;
        var currentStart = 0;

        for (int i = 0; i < line.length(); i++) {
            var currentChar = line.charAt(i);
            if (currentChar == '\t' || currentChar == '\n') {
                splitLine[currentIdx++] = line.substring(currentStart, i);

                currentStart = i + 1;
                if(currentIdx == 2)
                    break;
            }
        }
        splitLine[currentIdx] = line.substring(currentStart);

        data.getReadCountData()
                .computeIfAbsent(splitLine[0], id -> new ConcurrentHashMap<>())
                .computeIfAbsent(splitLine[1], id -> Integer.parseInt(splitLine[2]));
    }
}
