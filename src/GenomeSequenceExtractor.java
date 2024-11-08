import jdk.jshell.spi.ExecutionControl;
import model.FidxEntry;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class GenomeSequenceExtractor {
    private RandomAccessFile raf;
    Map<String, FidxEntry> fidxData;

    public GenomeSequenceExtractor(File fasta, Map<String, FidxEntry> fidxData) throws FileNotFoundException {
        raf = new RandomAccessFile(fasta, "r");
        this.fidxData = fidxData;
    }

    public String getSequence(String chr, int start, int end) throws IOException {
        var fidx = fidxData.get(chr);
        var lineLength = fidx.getLineLength();
        var lineLengthWithNewline = fidx.getLineLengthWithNewLine();
        var realStart = fidx.getStart() + start;
        var result = new StringBuilder();
        var newLineLength = lineLengthWithNewline - lineLength;

        var position = 0;
        raf.seek(realStart);
        var tempBuffer = new byte[lineLengthWithNewline];
        raf.readFully(tempBuffer);

        for (var b : tempBuffer) {
            var value = b & 0xFF;
            if (value < 32 || value == 127) {
                break;
            }
            position++;
        }

        raf.seek(realStart);
        tempBuffer = new byte[position];

        raf.readFully(tempBuffer);

        result.append(new String(tempBuffer));
        raf.skipBytes(newLineLength);

        var bytesToRead = end - start - position - (newLineLength - 1);
        var bytesRead = 0;

        while(bytesRead < bytesToRead){
            var currentBytesToRead = Math.min(lineLength, bytesToRead - bytesRead);
            if(tempBuffer.length != currentBytesToRead)
                tempBuffer = new byte[currentBytesToRead];

            raf.readFully(tempBuffer);
            result.append(new String(tempBuffer));
            raf.skipBytes(newLineLength);
            bytesRead += currentBytesToRead;
        }

        return result.toString();
    }


}
