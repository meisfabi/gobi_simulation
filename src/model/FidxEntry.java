package model;

public class FidxEntry {
    private String chromosome;
    private long start;
    private int lineLength;
    private int lineLengthWithNewLine;

    public FidxEntry(String chromosome, long start, int lineLength, int lineLengthWithNewLine) {
        this.chromosome = chromosome;
        this.start = start;
        this.lineLength = lineLength;
        this.lineLengthWithNewLine = lineLengthWithNewLine;
    }

    public long getStart() {
        return start;
    }

    public int getLineLength() {
        return lineLength;
    }

    public int getLineLengthWithNewLine() {
        return lineLengthWithNewLine;
    }

    public String getChromosome() {
        return chromosome;
    }
}
