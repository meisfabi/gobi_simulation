package model;

import augmentedTree.Interval;

public class Transcript implements Interval {

    private String transcriptId;
    private String transcriptName;
    private String transcriptSource;
    private TranscriptEntry transcriptEntry;
    private int start;
    private int stop;

    public TranscriptEntry getTranscriptEntry() {
        if(transcriptEntry == null)
            transcriptEntry = new TranscriptEntry();
        return transcriptEntry;
    }

    public void setTranscriptSource(String transcriptSource) {
        this.transcriptSource = transcriptSource;
    }

    public String getTranscriptSource() {
        return transcriptSource;
    }

    public void setTranscriptName(String transcriptName) {
        this.transcriptName = transcriptName;
    }

    public String getTranscriptName() {
        return transcriptName;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    @Override
    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }
}
