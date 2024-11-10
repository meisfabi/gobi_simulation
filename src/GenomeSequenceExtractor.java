import model.FidxEntry;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.Map;

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
        var result = new StringBuilder();
        var newLineLength = lineLengthWithNewline - lineLength;
        // var realStart = fidx.getStart();
        var linesTillStart = (start) / lineLength;
        var totalLinebreaks = linesTillStart * newLineLength;

        var realStart = fidx.getStart() + start + totalLinebreaks - 1;

        var position = 0;
        raf.seek(realStart);
        var tempBuffer = new byte[Math.min(end - start + 1, lineLength)];
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

        var bytesToRead = end - start - position - (newLineLength - 1) + 1;
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

    private final Map<Character, Character> complementMap = Map.of(
    'A', 'T',
    'T', 'A',
    'G', 'C',
    'C', 'G'
    );
    public String getReverseComplement(String gene, int start, int end){
        var geneLength = gene.length();
        var seqLength = end - start;
        var subsequence = gene.substring((geneLength - start - seqLength), (geneLength - start));
        var reverseComplement = new StringBuilder();
        for (int i = subsequence.length() - 1; i >= 0; i--) {
            char base = subsequence.charAt(i);
            reverseComplement.append(complementMap.get(base));
        }
        return reverseComplement.toString();
    }


}
