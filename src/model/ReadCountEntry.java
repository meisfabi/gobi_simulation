package model;

import pooling.ObjectFactory;

public class ReadCountEntry {
    private String geneId;
    private String transcriptId;
    private int readCounts;

    public ReadCountEntry(String geneId, String transcriptId, int readCounts) {
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.readCounts = readCounts;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public int getReadCounts() {
        return readCounts;
    }
}
