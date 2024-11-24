package model;

import augmentedTree.Interval;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class Gene implements Interval {
    private String geneId;
    private String geneName;
    private String geneSource;
    private String geneBiotype;
    private String seqName;
    private char strand;
    private int start;
    private int stop;

    private final Map<String, Transcript> transcriptMap = new ConcurrentHashMap<>();

    public Map<String, Transcript> getTranscriptMap() {
        return transcriptMap;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getGeneBiotype() {
        return geneBiotype;
    }

    public void setGeneBiotype(String geneBiotype) {
        this.geneBiotype = geneBiotype;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneSource() {
        return geneSource;
    }

    public void setGeneSource(String geneSource) {
        this.geneSource = geneSource;
    }
    public String getSeqName() {
        return seqName;
    }

    public void setSeqName(String seqName) {
        this.seqName = seqName;
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

    public char getStrand() {
        return strand;
    }

    public void setStrand(char strand) {
        this.strand = strand;
    }
}
